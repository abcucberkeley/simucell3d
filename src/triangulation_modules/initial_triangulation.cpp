#include "initial_triangulation.hpp"



/*
    This class has no constructor and all its methods are static. They can therefore
    be called without instantiating the class. 

    Given an arbitrary surface, this class returns a corresponding triangulated surface
    where the lengths of the edges are comprised in the range [l_min, l_max]. This class 
    does it by following the given procedure:

    1 - Coarsely triangulate the surface, by connecting the edges of each face to the face center
    2 - Sample the newly triangulated surface with the Poisson Disc Sampling algorithm
    3 - Reconstruct a finely triangulated surface with the Ball Pivoting Algorithm
    4 - Call the edge_manager to make sure that the lengths of the edges are in the prescribed range

    This class will throw an error if it receives a surface that is not closed

*/



//-----------------------------------------------------------------------------------------
mesh initial_triangulation::triangulate_surface(
    double const l_min, 
    double const l_max, 
    const mesh& surface, 
    const unsigned cell_id
) noexcept(false){

    //Makes a copy of the input mesh and work on the copy
    mesh surface_copy = surface;


    //Make sure the mesh does not have more points and faces than supported by the numbering system
    assert(surface_copy.get_nb_nodes() < std::numeric_limits<unsigned>::max());
    assert(surface_copy.get_nb_faces() < std::numeric_limits<unsigned>::max());

    //Check that the l_min and l_max make sense
    assert(l_min > 0.);
    assert(l_max > 0.);
    assert(l_max > l_min);

    //Make sure there are at least 4 nodes in the initial mesh
    if(!(surface_copy.get_nb_nodes() >= 4)){
        throw initial_triangulation_exception(std::string("There are less than 4 nodes in the input cell mesh: ") + std::to_string(surface_copy.get_nb_nodes()));
    } 

    //Make sure that all the faces have at least 3 nodes
    if(!(
            std::all_of(surface_copy.face_point_ids.begin(), surface_copy.face_point_ids.end(),
                [](const auto& f) -> bool {return f.size() >= 3;}
            )
        )
    ){
        throw initial_triangulation_exception("Not all the faces have at least 3 nodes");
    }


    //Triangulate the untriangulated mesh, acts in place on the mesh object
    coarse_triangulation(surface_copy);

    //Convert the mesh to a cell object. This allows to make sure that the face normals are oriented
    //correctly outward and that the cell is manifold (has no holes or self intersections)
    cell_ptr c1 = initial_triangulation::convert_mesh_to_cell(surface_copy);

    //Make sure the number of points sampled on the surface does not exceed the maximum number of points allowed 
    //by the numbering system
    const size_t estimated_nb_points =  static_cast<size_t>(c1->get_area() / (std::pow(l_max/4, 2) * M_PI));
    if(estimated_nb_points > std::numeric_limits<unsigned>::max()){
        throw initial_triangulation_exception(std::string("The estimated number of points sampled on the cell surface ("+std::to_string(estimated_nb_points)+") exceeds") +
            std::string(" the maximum number of points allowed by the numbering system ("+std::to_string(std::numeric_limits<unsigned>::max())+")"));
    }

    //Use the coarsely triangulated mesh to generate a Poisson point of the cell surface, where 
    //each point is separated by a distance of at least l_min from the rest
    const std::vector<oriented_point> poisson_point_cloud = initial_triangulation::generate_poisson_point_cloud(l_min, c1);

    //Get the bounding box of the cell
    const auto [min_x, min_y, min_z, max_x, max_y, max_z] = c1->get_aabb();

    //Use the ball pivoting algorithm to connect the points of the Poisson point cloud
    ball_pivoting_algorithm bpa(
        poisson_point_cloud, 
        l_min, 
        min_x, min_y, min_z, 
        max_x, max_y, max_z
    );

    //Run the ball pivoting algorithm
    bpa.run();



    //Get the ball pivoting triangulation
    std::vector<face> face_lst = bpa.get_face_lst();
    //std::vector<edge> edge_lst = bpa.get_edge_lst();
    std::vector<oriented_point> node_lst = bpa.get_node_lst();

    //Make a mesh object from the ball pivoting triangulation results
    mesh bp_triangulation;

    for(const auto& f: face_lst){
        const auto [n1_id, n2_id, n3_id] = f.get_node_ids();
        bp_triangulation.face_point_ids.push_back({n1_id, n2_id, n3_id});
    }

    for(const auto& n: node_lst){
        const auto [x, y, z] = n.position_.to_array();
        bp_triangulation.node_pos_lst.insert(bp_triangulation.node_pos_lst.end(), {x, y, z});
    }


    return bp_triangulation;

}
//-----------------------------------------------------------------------------------------



//-----------------------------------------------------------------------------------------
//Convert the mesh object into a cell. The cell class have a set of functions which automatically
//check the integrity of the input geometry and also orient correctly the face normals.
cell_ptr initial_triangulation::convert_mesh_to_cell(const mesh& surface) noexcept(false){ //throw mesh_integrity_exception


    cell_ptr c1 = std::make_shared<cell>(surface, 0);

    try{
        c1->initialize_cell_properties(false);
    }
    catch (const mesh_integrity_exception& e){

        //Rethrow the exception in another form, to make sure the program is killed
        throw initial_triangulation_exception(e.what());
    }

    return c1; 
}
//-----------------------------------------------------------------------------------------


//-----------------------------------------------------------------------------------------
//Use the coarsely triangulated mesh to generate a Poisson point of the cell surface, where 
//each point is separated by a distance of at least l_min from the rest
std::vector<oriented_point> initial_triangulation::generate_poisson_point_cloud(const double l_min, const cell_ptr c) noexcept{

    const std::vector<oriented_point> poisson_point_cloud = poisson_sampling::compute_poisson_point_cloud(l_min, c);

    return poisson_point_cloud;
}
//-----------------------------------------------------------------------------------------



//-----------------------------------------------------------------------------------------
//Trangulate the faces of the mesh that are not yet triangles by connecting each of their 
//edges to the center of the face
void initial_triangulation::coarse_triangulation(mesh& surface) noexcept{

    //Keep track of all the faces to delete
    std::vector<unsigned> face_to_delete_id_lst;
    std::vector<std::vector<unsigned>> face_to_add_lst;


    //Loop over the faces
    for(unsigned face_id = 0; face_id < surface.face_point_ids.size(); face_id++){
        std::vector<unsigned>& f =  surface.face_point_ids[face_id];
        assert(f.size() >= 3);

        //If the face is not alreadyy a triangle
        if(f.size() != 3){

            //This face will bee deleted later on
            face_to_delete_id_lst.push_back(face_id);

            //Compute the center of the face
            vec3 sum_node_pos;

            //Loop over the nodes of the face
            for(auto node_id: f){
                const auto [node_x, node_y, node_z] = surface.get_node_pos(node_id);

                //Sum their positions
                sum_node_pos.translate(node_x, node_y, node_z);
            }


            vec3 face_center = sum_node_pos / static_cast<double>(f.size());

            //Insert the centroid of the face in the mesh
            const auto face_center_pos = face_center.to_array();
            surface.node_pos_lst.insert(surface.node_pos_lst.end(), face_center_pos.begin(), face_center_pos.end());


            //Get the ids of this face_center in the mesh
            unsigned face_center_id = surface.get_nb_nodes() - 1;


            //Create the new faces
            for(unsigned i = f.size() -1, j = 0; j < f.size(); i=j++){
                face_to_add_lst.push_back({f[i], f[j], face_center_id});
            }
        }
    }
    
    //Add the newly created faces to the mesh 
    surface.face_point_ids.insert(surface.face_point_ids.end(), 
        std::make_move_iterator(face_to_add_lst.begin()), 
        std::make_move_iterator(face_to_add_lst.end())
    );

    //Remove all the faces that have been triangulated
    remove_index(surface.face_point_ids, face_to_delete_id_lst);


}
//-----------------------------------------------------------------------------------------





