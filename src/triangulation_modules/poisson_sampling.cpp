#include "poisson_sampling.hpp"



/*
    Given a triangulated surface, this class returns a point cloud where the points are located on the surface
    and are separated from each other by at least a given distance of l_min. The theory on the technique can be 
    found in the paper: "Parallel Poisson disk sampling with spectrum analysis on surfaces" Bowers et al, SIGGRAPH ASIA 2010.
    DOI: 978-1-4503-0439-9.

    The class is a static class and cannot be istantiated
*/


//-------------------------------------------------------------------------------------------------
//The input cell is already triangulated and all the normals of its faces point outward
std::vector<oriented_point> poisson_sampling::compute_poisson_point_cloud(double const l_min, const cell_ptr input_cell) noexcept(false){

    //Make sure we didn't get non sense for the l_min
    assert(l_min > 0.);

    //Make sure there are enough points in the mesh
    assert(input_cell->get_node_lst().size() >= 4);

    const auto face_lst = input_cell->get_face_lst();

    //Make sure that the area of each face has been calculated
    assert(
        std::all_of(face_lst.begin(), face_lst.end(), [](const face& f) -> bool {return f.get_area() > 0.;})
    );

    //Perform a uniform sampling of the cell surface
    std::vector<oriented_point>  uniform_point_cloud = poisson_sampling::uniform_sampling(input_cell, l_min);

    //Get the dimensions of the cell
    const auto [min_x, min_y, min_z, max_x, max_y, max_z] = input_cell->get_aabb();

    //Create 2 uniform space partitionning grids (uspg)
    const double voxel_size = l_min / std::sqrt(3);
    uspg_4d<oriented_point> grid_1(min_x, min_y, min_z, max_x, max_y, max_z, l_min, uniform_point_cloud.size());
    uspg_4d<oriented_point> grid_2(min_x, min_y, min_z, max_x, max_y, max_z, l_min, uniform_point_cloud.size());

    //Place all the uniformly sampled point in a uspg
    for(oriented_point& point: uniform_point_cloud){grid_1.place_object(point, point.position_);}

    //Use the first and second grid to efficiently compute a poisson point cloud where the points
    //lye on the cell surface
    auto poisson_point_cloud = poisson_sampling::poisson_disk_sampling(grid_1, grid_2, l_min);

    return poisson_point_cloud;
}
//-------------------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------------------
//The input mesh should already be triangulated, this version of the poisson point cloud 
//is used during the cell division. It makes sure that the nodes of the mesh that are located at the 
//division interface will be included in the poisson point cloud 
std::vector<oriented_point> poisson_sampling::compute_poisson_point_cloud(
    double const l_min, 
    mesh& m,
    const unsigned division_point_ids_threshold,  // All the points with an id equal or greater than this value should be added to the sampling
    const unsigned division_face_ids_threshold,   //All the faces with an id equal or greater than this value will be sampled
    const vec3 division_plane_normal //The normal of the division plane
) noexcept(false){

    assert(division_point_ids_threshold < m.node_pos_lst.size() / 3);
    assert(division_face_ids_threshold < m.face_point_ids.size());
    assert(division_face_ids_threshold > 0);

    //Make sure we didn't get non sense for the l_min
    assert(l_min > 0.);

    //Make sure there are enough points in the mesh
    assert(static_cast<unsigned>(m.node_pos_lst.size() / 3) >= 4);

    //Make sure that the area of each face has been calculated
    assert(
        std::all_of(m.face_point_ids.begin(), m.face_point_ids.end(), [](const auto& f) -> bool {return f.size() == 3;})
    );

    //Perform a uniform sampling of the cell surface
    std::vector<oriented_point>  uniform_point_cloud = poisson_sampling::uniform_sampling(m, division_face_ids_threshold, l_min);

    //Get the dimensions of the surface to sample
    auto [min_x, min_y, min_z, max_x, max_y, max_z] = poisson_sampling::get_surface_aabb(m, division_face_ids_threshold);

    //Add some padding to the min_z and max_z otherwise they are both equal to 0
    min_z = min_z - (2. * l_min);
    max_z = max_z + (2. * l_min);


    //Create 2 uniform space partitionning grids (uspg)
    const double voxel_size = l_min / std::sqrt(3);
    uspg_4d<oriented_point> grid_1(min_x, min_y, min_z, max_x, max_y, max_z, l_min, uniform_point_cloud.size());
    uspg_4d<oriented_point> grid_2(min_x, min_y, min_z, max_x, max_y, max_z, l_min, uniform_point_cloud.size());

    //Place all the uniformly sampled point in a uspg
    for(oriented_point& point: uniform_point_cloud){grid_1.place_object(point, point.position_);}

    //Already include the nodes that should be included in the final poisson point cloud
    for(unsigned p_id = division_point_ids_threshold; p_id < m.node_pos_lst.size() / 3; p_id++){
        oriented_point new_point(vec3(m.node_pos_lst[p_id * 3], m.node_pos_lst[p_id * 3 + 1], m.node_pos_lst[p_id * 3 + 2]), division_plane_normal);
        new_point.created_by_poisson_sampling_ = false;
        grid_2.place_object(new_point, new_point.position_);
    }

    const auto grid_content = grid_2.get_grid_content();
    const size_t number_points = std::distance(grid_content.begin(), grid_content.end());

    //Use the first and second grid to efficiently compute a poisson point cloud where the points
    //lye on the cell surface
    auto poisson_point_cloud = poisson_sampling::poisson_disk_sampling(grid_1, grid_2, l_min);

    return poisson_point_cloud;
}
//-------------------------------------------------------------------------------------------------




//-------------------------------------------------------------------------------------------------
//Uniformaly sample the cell surface
std::vector<oriented_point> poisson_sampling::uniform_sampling(const cell_ptr input_cell, const double l_min) noexcept(false){

    //Calculate the area of a Poisson Disc
    const double disk_area = M_PI * std::pow(l_min, 2);

    //The list in which the newly sampled points will be stored
    std::vector<oriented_point> uniform_point_cloud;

    //Use random numbers to sample uniformly the cell surface
    std::minstd_rand gen(std::chrono::system_clock::now().time_since_epoch().count());
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    const auto face_lst = input_cell->get_face_lst();

    //Parallelize the random sampling of nodes, by allowing mutliple threads to write into the same vector
    #pragma omp declare reduction (merge : std::vector<oriented_point> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))

    //In parallel, perform the uniform random sampling
    #pragma omp parallel for reduction(merge: uniform_point_cloud)
    for(size_t i = 0; i < face_lst.size(); i++){

        const face& f = face_lst[i];

        //Get the area of the face
        const double f_area = f.get_area();

        //Get the points of the face
        const auto [n1_id, n2_id, n3_id] = f.get_node_ids();
        const node& n1 = input_cell->get_const_ref_node(n1_id);
        const node& n2 = input_cell->get_const_ref_node(n2_id);
        const node& n3 = input_cell->get_const_ref_node(n3_id);

        //Get the normal of the face
        const vec3 face_normal = f.get_normal();

        //Compute the number of points to insert to fully cover the triangle surface
        size_t nb_points_to_insert = static_cast<size_t>(f_area * 1000  / disk_area);

        for(size_t i = 0; i < nb_points_to_insert; i++){

            //Get 2 random numbers to generate random barycentric coodiates
            const double rn_1 = dist(gen), rn_2 = dist(gen);

            //Generate barycentric coordinates with this 2 random numbers
            const double root = std::sqrt(rn_1);
            assert(std::isfinite(root));

            const double u = 1. - root;
            const double v = rn_2 * root;
            const double w = 1. - (u + v);

            //Get a random point on the surface of the triangle
            const vec3 rn_point = (n1.pos() * u) + (n2.pos() * v) + (n3.pos() * w); 

            //Save the random point
            uniform_point_cloud.emplace_back(rn_point, face_normal);
        }
    }
    return uniform_point_cloud;
}
//-------------------------------------------------------------------------------------------------




//-------------------------------------------------------------------------------------------------
//Uniformaly sample the cell surface
std::vector<oriented_point> poisson_sampling::uniform_sampling(
    const mesh& m, 
    const unsigned division_face_ids_threshold, //All the faces with id >= to this value will be sampled
    const double l_min
) noexcept(false){
    assert(division_face_ids_threshold < m.face_point_ids.size());

    //Calculate the area of a Poisson Disc
    const double disk_area = M_PI * std::pow(l_min, 2);

    //The list in which the newly sampled points will be stored
    std::vector<oriented_point> uniform_point_cloud;

    //Use random numbers to sample uniformly the cell surface
    std::minstd_rand gen(std::chrono::system_clock::now().time_since_epoch().count());
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    //Loop over the faces of the mesh that should be triangulated
    for(unsigned f_id = division_face_ids_threshold; f_id < m.face_point_ids.size(); f_id++){
        assert(f_id < m.face_point_ids.size());
        const auto& f = m.face_point_ids[f_id];
        assert(f.size() == 3);

        //Get the points of the face
        const unsigned n_1_id = f[0];
        const unsigned n_2_id = f[1];
        const unsigned n_3_id = f[2];

        //Get the coordinates of the points
        assert(n_1_id < m.get_nb_nodes() && n_2_id < m.get_nb_nodes() && n_3_id < m.get_nb_nodes());
        const vec3 n1(m.node_pos_lst[n_1_id*3    ], m.node_pos_lst[n_1_id*3 + 1], m.node_pos_lst[n_1_id*3 + 2]);
        const vec3 n2(m.node_pos_lst[n_2_id*3    ], m.node_pos_lst[n_2_id*3 + 1], m.node_pos_lst[n_2_id*3 + 2]);
        const vec3 n3(m.node_pos_lst[n_3_id*3    ], m.node_pos_lst[n_3_id*3 + 1], m.node_pos_lst[n_3_id*3 + 2]);

        const vec3 face_normal = (n2 - n1).cross(n3 - n1);
        const double f_area = face_normal.norm() / 2.;
        assert(std::isfinite(f_area));

        //Compute the number of points to insert to fully cover the triangle surface
        size_t nb_points_to_insert = static_cast<size_t>(f_area * 1000  / disk_area);

        if(nb_points_to_insert > std::numeric_limits<unsigned>::max()){
            throw division_exception("The number of points to insert is too large");
        }

        for(size_t i = 0; i < nb_points_to_insert; i++){

            //Get 2 random numbers to generate random barycentric coodiates
            const double rn_1 = dist(gen), rn_2 = dist(gen);

            //Generate barycentric coordinates with this 2 random numbers
            const double root = std::sqrt(rn_1);
            assert(std::isfinite(root));

            const double u = 1. - root;
            const double v = rn_2 * root;
            const double w = 1. - (u + v);

            //Get a random point on the surface of the triangle
            const vec3 rn_point = (n1 * u) + (n2 * v) + (n3 * w); 

            //Save the random point
            uniform_point_cloud.emplace_back(rn_point, face_normal);
        }
    }
    
    return uniform_point_cloud;
}
//-------------------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------------------
//Carry out the Poisson Disk Sampling of the cell surface
std::vector<oriented_point> poisson_sampling::poisson_disk_sampling(
    const uspg_4d<oriented_point>& grid_1, 
    uspg_4d<oriented_point>& grid_2,
    const double l_min
    ) noexcept{

    const double l_min_squared = l_min * l_min;

    //Get the nb of voxels in each direction in the first grid
    const auto [nb_voxel_x, nb_voxel_y, nb_voxel_z] = grid_1.get_nb_voxels(); 

    //Loop over the voxels of the first grid
    for(unsigned voxel_x = 0; voxel_x < nb_voxel_x; voxel_x ++){
    for(unsigned voxel_y = 0; voxel_y < nb_voxel_y; voxel_y ++){
    for(unsigned voxel_z = 0; voxel_z < nb_voxel_z; voxel_z ++){

        //Get the nodes in that voxel
        const auto voxel_content_lst = grid_1.get_voxel_content(voxel_x, voxel_y, voxel_z);

        //Get all the nodes in the neighborhood of the corresponding voxel in the second grid
        const auto poisson_point_lst = grid_2.get_neighborhood(voxel_x, voxel_y, voxel_z);

        //Try a max of 30 times to insert a node of voxel_content_lst into the voxel (voxel_x, voxel_y, voxel_z) of grid_2
        const auto voxel_content_size = std::distance(voxel_content_lst.begin(), voxel_content_lst.end());
        for(size_t point_id = 0; point_id < voxel_content_size && point_id < 30; point_id++){

            //Get the the candidate point for insertion
            oriented_point candidate_point = *(std::next(voxel_content_lst.begin(), point_id));


            //Check that this candidate point is at least at a distance of l_min from the points already inserted
            //in grid_2
            bool insert_candidate_point = true;
            unsigned nb_tries = 0;
            //Tries 30 times to insert the point in the grid
            for(oriented_point poisson_point: poisson_point_lst){
                if((poisson_point.position_ - candidate_point.position_).squared_norm() < l_min_squared  || nb_tries >= 30){
                    insert_candidate_point = false;
                    break;
                } 
                nb_tries++;             
            }

            if(insert_candidate_point){
                //Plcace the point in the second grid and go to the next voxel of grid_1
                grid_2.place_object(candidate_point, candidate_point.position_);
                break;
            }
        }
    }}}

    //Get the poisson point cloud
    const auto poisson_point_cloud_lst = grid_2.get_grid_content();

    //Make a copy of the points
    std::vector<oriented_point> poisson_point_cloud_lst_copy;
    poisson_point_cloud_lst_copy.reserve(std::distance(poisson_point_cloud_lst.begin(), poisson_point_cloud_lst.end()));
    for(oriented_point point: poisson_point_cloud_lst){poisson_point_cloud_lst_copy.push_back(point);}
    
    //And return the copy
    return poisson_point_cloud_lst_copy;
}
//-------------------------------------------------------------------------------------------------



//-------------------------------------------------------------------------------------------------
//Utility function to get the AABB of the set of faces of a mesh
std::array<double, 6> poisson_sampling::get_surface_aabb(
        const mesh& m, 
        const std::vector<unsigned> face_to_sample_ids
    ) noexcept{


    double aabb_min_x = std::numeric_limits<double>::infinity();
    double aabb_min_y = std::numeric_limits<double>::infinity();
    double aabb_min_z = std::numeric_limits<double>::infinity();

    double aabb_max_x = -std::numeric_limits<double>::infinity();
    double aabb_max_y = -std::numeric_limits<double>::infinity();
    double aabb_max_z = -std::numeric_limits<double>::infinity();


    //Loop over the faces to sample
    for(const unsigned  face_id: face_to_sample_ids){

        //Get the face of the mesh 
        assert(face_id < m.face_point_ids.size());
        const auto f = m.face_point_ids[face_id];

        //Make sure the face is a triangle 
        assert(f.size() == 3);

        //Get the nodes of the face
        const unsigned n_1_id = f[0];
        const unsigned n_2_id = f[1];
        const unsigned n_3_id = f[2];

        //Make sure the nodes are in the mesh
        assert(n_1_id < m.node_pos_lst.size()/3);
        assert(n_2_id < m.node_pos_lst.size()/3);
        assert(n_3_id < m.node_pos_lst.size()/3);

        //Get the min and max coordinates of the nodes
        const double x_min_face = std::min(m.node_pos_lst[n_1_id*3    ], std::min(m.node_pos_lst[n_2_id*3    ], m.node_pos_lst[n_3_id*3    ]));
        const double x_max_face = std::max(m.node_pos_lst[n_1_id*3    ], std::max(m.node_pos_lst[n_2_id*3    ], m.node_pos_lst[n_3_id*3    ]));

        const double y_min_face = std::min(m.node_pos_lst[n_1_id*3 + 1], std::min(m.node_pos_lst[n_2_id*3 + 1], m.node_pos_lst[n_3_id*3 + 1]));
        const double y_max_face = std::max(m.node_pos_lst[n_1_id*3 + 1], std::max(m.node_pos_lst[n_2_id*3 + 1], m.node_pos_lst[n_3_id*3 + 1]));

        const double z_min_face = std::min(m.node_pos_lst[n_1_id*3 + 2], std::min(m.node_pos_lst[n_2_id*3 + 2], m.node_pos_lst[n_3_id*3 + 2]));
        const double z_max_face = std::max(m.node_pos_lst[n_1_id*3 + 2], std::max(m.node_pos_lst[n_2_id*3 + 2], m.node_pos_lst[n_3_id*3 + 2]));

        if(x_min_face < aabb_min_x){aabb_min_x = x_min_face;}
        if(y_min_face < aabb_min_y){aabb_min_y = y_min_face;}
        if(z_min_face < aabb_min_z){aabb_min_z = z_min_face;}

        if(x_max_face > aabb_max_x){aabb_max_x = x_max_face;}
        if(y_max_face > aabb_max_y){aabb_max_y = y_max_face;}
        if(z_max_face > aabb_max_z){aabb_max_z = z_max_face;}
    }

    

    return {aabb_min_x, aabb_min_y, aabb_min_z, aabb_max_x, aabb_max_y, aabb_max_z};
}
//-------------------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------------------
//Utility function to get the AABB of the set of faces of a mesh
std::array<double, 6> poisson_sampling::get_surface_aabb(
    const mesh& m, 
    const unsigned division_face_ids_threshold // All the faces with id >= division_face_ids_threshold will be included in the AABB
) noexcept{

    assert(division_face_ids_threshold < m.face_point_ids.size());

    double aabb_min_x = std::numeric_limits<double>::infinity();
    double aabb_min_y = std::numeric_limits<double>::infinity();
    double aabb_min_z = std::numeric_limits<double>::infinity();

    double aabb_max_x = -std::numeric_limits<double>::infinity();
    double aabb_max_y = -std::numeric_limits<double>::infinity();
    double aabb_max_z = -std::numeric_limits<double>::infinity();


    //Loop over the faces to sample
    for(unsigned face_id = division_face_ids_threshold; face_id < m.face_point_ids.size(); face_id++){

        //Get the face of the mesh 
        assert(face_id < m.face_point_ids.size());
        const auto f = m.face_point_ids[face_id];

        //Make sure the face is a triangle 
        assert(f.size() == 3);

        //Get the nodes of the face
        const unsigned n_1_id = f[0];
        const unsigned n_2_id = f[1];
        const unsigned n_3_id = f[2];

        //Make sure the nodes are in the mesh
        assert(n_1_id < m.node_pos_lst.size()/3);
        assert(n_2_id < m.node_pos_lst.size()/3);
        assert(n_3_id < m.node_pos_lst.size()/3);

        //Get the min and max coordinates of the nodes
        const double x_min_face = std::min(m.node_pos_lst[n_1_id*3    ], std::min(m.node_pos_lst[n_2_id*3    ], m.node_pos_lst[n_3_id*3    ]));
        const double x_max_face = std::max(m.node_pos_lst[n_1_id*3    ], std::max(m.node_pos_lst[n_2_id*3    ], m.node_pos_lst[n_3_id*3    ]));

        const double y_min_face = std::min(m.node_pos_lst[n_1_id*3 + 1], std::min(m.node_pos_lst[n_2_id*3 + 1], m.node_pos_lst[n_3_id*3 + 1]));
        const double y_max_face = std::max(m.node_pos_lst[n_1_id*3 + 1], std::max(m.node_pos_lst[n_2_id*3 + 1], m.node_pos_lst[n_3_id*3 + 1]));

        const double z_min_face = std::min(m.node_pos_lst[n_1_id*3 + 2], std::min(m.node_pos_lst[n_2_id*3 + 2], m.node_pos_lst[n_3_id*3 + 2]));
        const double z_max_face = std::max(m.node_pos_lst[n_1_id*3 + 2], std::max(m.node_pos_lst[n_2_id*3 + 2], m.node_pos_lst[n_3_id*3 + 2]));

        if(x_min_face < aabb_min_x){aabb_min_x = x_min_face;}
        if(y_min_face < aabb_min_y){aabb_min_y = y_min_face;}
        if(z_min_face < aabb_min_z){aabb_min_z = z_min_face;}

        if(x_max_face > aabb_max_x){aabb_max_x = x_max_face;}
        if(y_max_face > aabb_max_y){aabb_max_y = y_max_face;}
        if(z_max_face > aabb_max_z){aabb_max_z = z_max_face;}
    }

    return {aabb_min_x, aabb_min_y, aabb_min_z, aabb_max_x, aabb_max_y, aabb_max_z};
}


//-------------------------------------------------------------------------------------------------
