#include "cell_divider.hpp"





//-----------------------------------------------------------------------------------------------
//If a cell is ready to be divided, this function will divide it into two daughter cells
void cell_divider::run(
    std::vector<cell_ptr>& cell_lst,
    const double l_min, //The min edge length
    const local_mesh_refiner& lmr, // Will be used to refine the meshes of the daughter cells
    unsigned& max_cell_id_, //Use this number to assign unique ids to the daughter cells
    bool verbose
) noexcept{


    //Keep track of the ids of the cells to delete 
    std::vector<unsigned> cells_to_delete_lst;

    //Loop over the cells
    #pragma omp parallel for
    for(size_t i = 0; i < cell_lst.size(); i++){

        bool is_ready_to_divide = cell_lst[i]->is_ready_to_divide();

        //If the cell is ready to be divided
        if(cell_lst[i]->is_ready_to_divide()){


            auto division_result  = divide_cell(cell_lst[i], l_min, lmr);

            //If the division was successful
            if(division_result.has_value()){

                #pragma omp critical
                {
                    auto [daughter_1, daughter_2] = division_result.value();

                    //Erase all the faces and nodes of the mother cell to avoid cyclic dependencies
                    cell_lst[i]->clear_data(); 

                    cells_to_delete_lst.push_back(i);

                    daughter_1->cell_id_ = max_cell_id_++;
                    daughter_2->cell_id_ = max_cell_id_++;                
                    
                    cell_lst.push_back(daughter_1);
                    cell_lst.push_back(daughter_2);
                } 
            }
        }
    }

    //Remove the mother cells at the end of the parallel loop 
    if(cells_to_delete_lst.size() > 0){
        std::sort(cells_to_delete_lst.begin(), cells_to_delete_lst.end(), std::less<unsigned>());
        remove_index(cell_lst, cells_to_delete_lst);

        for(size_t local_cell_id = 0; local_cell_id < cell_lst.size(); local_cell_id++){
            cell_lst[local_cell_id]->set_local_id(local_cell_id);
        }
    }




}
//-----------------------------------------------------------------------------------------------



//-----------------------------------------------------------------------------------------------
//Run the whole division pipeline, and cut the mother cell into 2 daughter cells
//The mother cell is not affected by the division operation, if anything goes wrong during the division
//operation the mother cell will continue to exist
std::optional<std::pair<cell_ptr, cell_ptr>>  cell_divider::divide_cell(
    cell_ptr c, 
    const double l_min, //The min edge length
    const local_mesh_refiner& lmr // Will be used to refine the meshes of the daughter cells
) noexcept{


    try{
        //Remove the unused nodes or faces in the cell
        c->rebase();
        
        //Get the origin and the normal of the division plane
        const vec3 centroid = c->compute_centroid();
        const vec3 division_plane_normal = c->get_cell_division_axis();

        //All the points that will be part of the division interface will have ids >= than this value
        const unsigned division_point_ids_threshold = c->get_node_lst().size();

        //Add the points where the plane intersects with the edges of the cell
        mesh m = cell_divider::add_intersection_points(c, centroid, division_plane_normal);

        //At this stage some faces are not triangles anymore and must be triangulated before the division can be performed
        cell_divider::divide_faces(m, division_point_ids_threshold);

        //All the faces that will be part of the division interface will have an id greater than this value
        const unsigned division_face_ids_threshold  = m.face_point_ids.size();

        //Form a polygon by connecting all the points where the plane intersects with the edges of the cell
        // and then divide this polygon in triangles
        const unsigned nb_division_face_points = m.node_pos_lst.size() / 3 - division_point_ids_threshold;
        std::vector<unsigned> division_face_point_ids(nb_division_face_points);
        std::iota(division_face_point_ids.begin(), division_face_point_ids.end(), division_point_ids_threshold);
        m.face_point_ids.push_back(division_face_point_ids);
        initial_triangulation::coarse_triangulation(m);
    
        //Map the points of the division face to the xy plane. This is done such that we can use the 2D
        //Delaunay triangulation and it also facilitates subsequent tests
        auto [translation, rotation] = cell_divider::map_points_to_xy_plane(m, division_point_ids_threshold, division_plane_normal);

        //Triangulate the the division interface in a way that respects the constraints on the edge lengths
        cell_divider::triangulate_division_interface(
            l_min, 
            m, 
            division_point_ids_threshold, 
            division_face_ids_threshold, 
            division_plane_normal
        );

        //Map the points of the division interface back to their original positions
        cell_divider::map_points_to_division_plane(m, division_point_ids_threshold, translation, rotation);

        //Split the mesh of the mother cell in 2 and create the daughter cells
        auto [daughter_cell_1, daughter_cell_2] = cell_divider::create_daughter_cells(
            c, 
            m, 
            division_point_ids_threshold, 
            division_face_ids_threshold, 
            division_plane_normal,
            centroid
        );




        //Refine the meshes of the daughter cells to make sure their edge lengths are in the range [l_min, l_max]
        lmr.refine_mesh(daughter_cell_1);
        lmr.refine_mesh(daughter_cell_2);

        //If this is not done there will be a sudden drop in the daughter cells pressures
        daughter_cell_1->target_volume_ = c->target_volume_ / 2;
        daughter_cell_2->target_volume_ = c->target_volume_ / 2;

        daughter_cell_1->rebase();
        daughter_cell_2->rebase();

        //Draw the growth rate and division volume from a normal distribution
        daughter_cell_1->initialize_random_properties();
        daughter_cell_2->initialize_random_properties();

        return std::make_pair(daughter_cell_1, daughter_cell_2);
    }

    catch(std::exception& e){return std::nullopt;}
}
//-----------------------------------------------------------------------------------------------






//-----------------------------------------------------------------------------------------------
//Finds where the division plane intersects with the edges of the cell and return a mesh of the 
//cell with the intersection points added.
mesh cell_divider::add_intersection_points(
    cell_ptr c, 
    const vec3& p,  //Position of a point on the plane
    const vec3& n   //Normal vector of the plane
)noexcept(false){

    assert(c->free_face_queue_.empty());
    assert(c->free_node_queue_.empty());


    //Get a mesh object of the cell, all the modifications are going to be done on this mesh object
    //and the mother cell will be left untouched, if anything goes wrong during the division operation
    //the mother cell will continue to exist
    mesh m = c->get_mesh();

    //Find one first edge that is intersected by the plane
    std::optional<edge> seed_edge_opt;
    for(const edge& e : c->get_edge_set()){

        //Get the nodes of the edge
        const vec3& n1_pos = c->get_node(e.n1()).pos();
        const vec3& n2_pos = c->get_node(e.n2()).pos();
        std::optional<vec3> intersection_point_opt = find_edge_plane_intersection(n1_pos, n2_pos, p, n);

        if(intersection_point_opt.has_value()){
            seed_edge_opt = e;
            const vec3 intersection_point = intersection_point_opt.value();
            m.node_pos_lst.insert(m.node_pos_lst.end(), {intersection_point.dx(), intersection_point.dy(), intersection_point.dz()});
            
            //Add the intersection point to the 2 faces that share the edge
            const unsigned f_1_id = e.f1();
            const unsigned f_2_id = e.f2();

            //Get the id of the intersection point
            const unsigned intersection_point_id = static_cast<double>(m.node_pos_lst.size()/3.) - 1;
     
            //Add the intersection point to the 2 faces
            add_point_to_face(m, f_1_id, e.n1(), e.n2(), intersection_point_id);
            add_point_to_face(m, f_2_id, e.n1(), e.n2(), intersection_point_id);
            
            break;
        }
    }
    if(!seed_edge_opt.has_value()){throw division_exception("The plane does not intersect with the cell");}
    edge edge_to_cut = seed_edge_opt.value();

    //We cut the cell starting from f2 and when the cut reaches f1 we stop
    const face stop_face   = c->get_face(edge_to_cut.f1());
    face face_to_cut       = c->get_face(edge_to_cut.f2());

    //This loop goes through all the edges that are intersected by the plane
    unsigned iteration = 0;
    while(face_to_cut != stop_face){

        //Get the ids of the nodes of the edge that has already been cut
        const unsigned n_a_id = edge_to_cut.n1();
        const unsigned n_b_id = edge_to_cut.n2();

        //Get the node opposite to the edge that has already been cut in the face
        const unsigned n_c_id = face_to_cut.get_opposite_node(edge_to_cut.n1(), edge_to_cut.n2());

        //Get the nodes a, b, and c
        const node& n_a = c->get_node(n_a_id);
        const node& n_b = c->get_node(n_b_id);
        const node& n_c = c->get_node(n_c_id);

        //The plane either cuts the edge ac or the edge bc
        std::optional<vec3> intersection_ac = find_edge_plane_intersection(n_a.pos(), n_c.pos(), p, n);
        std::optional<vec3> intersection_bc = find_edge_plane_intersection(n_b.pos(), n_c.pos(), p, n);

        assert((intersection_ac.has_value() && intersection_bc.has_value()) == false);
        vec3 intersection_point;
        //If the plane cuts the edge ac
        if(intersection_ac.has_value()){
            auto edge_to_cut_opt = c->get_edge(n_a_id, n_c_id);
            assert(edge_to_cut_opt.has_value());
            edge_to_cut = edge_to_cut_opt.value();
            intersection_point = intersection_ac.value();
        }

        //If the plane cuts the edge bc
        else if(intersection_bc.has_value()){
            auto edge_to_cut_opt = c->get_edge(n_b_id, n_c_id);
            assert(edge_to_cut_opt.has_value());
            edge_to_cut = edge_to_cut_opt.value();
            intersection_point = intersection_bc.value();
        }
        else{throw division_exception("The plane does not intersect with the cell");}

        //Add the inetsection point to the mesh object
        m.node_pos_lst.insert(m.node_pos_lst.end(), {intersection_point.dx(), intersection_point.dy(), intersection_point.dz()});

        //Add the intersection point to the 2 faces
        const unsigned intersection_point_id = static_cast<double>(m.node_pos_lst.size()/3.) - 1;
        add_point_to_face(m, edge_to_cut.f1(), edge_to_cut.n1(), edge_to_cut.n2(), intersection_point_id);
        add_point_to_face(m, edge_to_cut.f2(), edge_to_cut.n1(), edge_to_cut.n2(), intersection_point_id);

        //Get the next face which is cut by the plane
        const unsigned face_to_cut_id = (edge_to_cut.f1() == face_to_cut.get_local_id()) ? edge_to_cut.f2() : edge_to_cut.f1();
        face_to_cut = c->get_face(face_to_cut_id);

        if(iteration == c->get_edge_set().size() - 1){throw division_exception("Infinite loop stopped during cell division");}
        iteration++;
    }

    return m; 
}

//-----------------------------------------------------------------------------------------------




//-----------------------------------------------------------------------------------------------
//At this stage of the division process, some faces of the mesh have 5 nodes. This is because the plane
//has intersected with 2 of their edges. This function splits those faces in 3 triangular subfaces. 
//The argument "intersection_point_ids_threshold" corresponds to the id of the first intersection point
//that was added to the mesh. The other intersection points that were added to the mesh have ids greater
void cell_divider::divide_faces(mesh& m, const unsigned intersection_point_ids_threshold) noexcept(false){
    assert(intersection_point_ids_threshold > 0 && intersection_point_ids_threshold < m.node_pos_lst.size()/3);
    assert(m.face_point_ids.size() > 0);

    //Keep track of the faces to delete
    std::vector<unsigned> faces_to_delete;
    std::vector<std::vector<unsigned>> faces_to_add;

    //Loop through all the faces of the mesh
    for(unsigned k = 0; k < m.face_point_ids.size(); k++){
        auto& f = m.face_point_ids[k];
        assert(f.size() == 3 || f.size() == 5);

        //If the face has 5 nodes, it means that the plane has intersected with 2 of its edges
        if(f.size() == 5){
            faces_to_delete.push_back(k);

            //Get the ids of the intersection points
            auto it1 = std::find_if(f.begin(), f.end(), [intersection_point_ids_threshold](const unsigned& id){return id >= intersection_point_ids_threshold;});
            assert(it1 != f.end());

            auto it2 = std::find_if(it1 + 1, f.end(), [intersection_point_ids_threshold](const unsigned& id){return id >= intersection_point_ids_threshold;});
            assert(it2 != f.end());
            assert(*it1 != *it2);

            //Get the local ids of the 2 intesection pointspoints in the face
            const unsigned local_point_1_id = std::distance(f.begin(), it1);
            const unsigned local_point_2_id = std::distance(f.begin(), it2);
            assert(local_point_2_id - local_point_1_id == 3 || local_point_2_id - local_point_1_id == 2);

            //Based on where the face has been divided there are 2 ways to create the new triangular subfaces
            if(local_point_2_id - local_point_1_id  == 3){

                //Add the new triangular subfaces
                faces_to_add.insert(faces_to_add.end(), {f[local_point_1_id],  f[local_point_2_id],              f[(local_point_2_id + 1) % 5]});
                faces_to_add.insert(faces_to_add.end(), {f[local_point_1_id],  f[local_point_1_id + 1],          f[local_point_1_id + 2]});
                faces_to_add.insert(faces_to_add.end(), {f[local_point_1_id],  f[(local_point_1_id + 2) % 5],    f[local_point_2_id]});

            }
            else{
                faces_to_add.insert(faces_to_add.end(), {f[local_point_1_id],            f[local_point_1_id + 1],    f[local_point_2_id]});
                faces_to_add.insert(faces_to_add.end(), {f[(local_point_2_id + 2) % 5],  f[local_point_1_id],        f[local_point_2_id]});
                faces_to_add.insert(faces_to_add.end(), {f[(local_point_2_id + 2) % 5],  f[local_point_2_id],        f[(local_point_2_id + 1) % 5]});
            }
        }
    }

    
    //Delete the faces that have been divided
    std::sort(faces_to_delete.begin(), faces_to_delete.end(), std::less<unsigned>());
    remove_index(m.face_point_ids, faces_to_delete);

    //Add the faces that have been created 
    std::move(faces_to_add.begin(), faces_to_add.end(), std::back_inserter(m.face_point_ids));
}
//-----------------------------------------------------------------------------------------------






//-----------------------------------------------------------------------------------------------
//Find the intersection point between an edge (line segment) and a plane
//If there is no intersection returns an empty optional
std::optional<vec3> cell_divider::find_edge_plane_intersection(
    const vec3& e1, //Position of the 1st point of the edge
    const vec3& e2, //Position of the 2nd point of the edge

    const vec3& p,  //Position of a point on the plane
    const vec3& n   //Normal vector of the plane
) noexcept(false){


    //Find where a plane cuts an edge, returns nullptr if it doesn't cut it
    // Compute the t value for the directed line ab intersecting the plane
    const double dot1 = n.dot(p  - e1);
    const double dot2 = n.dot(e2 - e1);

    double t;
    if(dot2 == 0.0){return std::nullopt;} //The edge is colinear with the plane
    else{t = dot1 / dot2;}
    if(t < 0.0 || t > 1.0){return std::nullopt;}
    
    vec3 intersection_point = e1 + (e2 - e1) * t;
    return intersection_point;
}
//-----------------------------------------------------------------------------------------------







//-----------------------------------------------------------------------------------------------
//Based on the side of a face wrt the plane, return true or false based on the side
bool cell_divider::face_side_wrt_plane(
    const face& f,
    const cell_ptr c,

    const vec3& p,  //Position of a point on the plane
    const vec3& n   //Normal vector of the plane
    
) noexcept{

    //Get the position of the 3 nodes of the face
    auto [n1_id, n2_id, n3_id] = f.get_node_ids();

    const vec3& p1 = c->get_node(n1_id).pos();
    const vec3& p2 = c->get_node(n2_id).pos();
    const vec3& p3 = c->get_node(n3_id).pos();

    //Compute the centroid of the face
    const vec3 centroid = (p1 + p2 + p3) / 3.0;

    return (centroid - p).dot(n) > 0;
}


//-----------------------------------------------------------------------------------------------



//-----------------------------------------------------------------------------------------------
//Add the id of the point p to the given face and make sure that the point p is inserted between the nodes n_a and n_b
void cell_divider::add_point_to_face(
    mesh& m, 
    const unsigned face_id, 
    const unsigned n_a_id, 
    const unsigned n_b_id,
    const unsigned n_p_id) noexcept(false){

    assert(face_id < m.face_point_ids.size());

    //Get the face
    auto& f = m.face_point_ids[face_id];
    assert(f.size() >= 3);

    //Loop over the edges of the face
    for(unsigned i = f.size() -1, j = 0; j < f.size(); i = j++){
        //If the edge is the one between n_a and n_b
        if((f[i] == n_a_id && f[j] == n_b_id) || (f[i] == n_b_id && f[j] == n_a_id)){
            //Insert the point between the two nodes
            f.insert(f.begin() + j, n_p_id);
            return;
        }
    }
    throw division_exception("The point could not be inserted between the two nodes");
}
//-----------------------------------------------------------------------------------------------





//-----------------------------------------------------------------------------------------------
//As it names indicates, this function maps the points in intersection_points to the xy plane
//Returns the translation and the rotation that were applied to the points
std ::pair<vec3, mat33> cell_divider::map_points_to_xy_plane(
    mesh& m, 
    const unsigned division_point_ids_threshold,
    const vec3& division_plane_normal 
) noexcept{

    assert(division_point_ids_threshold < m.node_pos_lst.size() / 3);

    //Compute the number of points that are on the division interface
    const double nb_points_division_interface = (static_cast<double>(m.node_pos_lst.size()) / 3.) - static_cast<double>(division_point_ids_threshold);

    //First compute the centroid of the intersection points
    vec3 centroid(0.0, 0.0, 0.0);
    for(unsigned p_id = division_point_ids_threshold; p_id < m.node_pos_lst.size() / 3; p_id++){
        centroid.translate(m.node_pos_lst[p_id *3    ], m.node_pos_lst[p_id *3 + 1], m.node_pos_lst[p_id *3 + 2]);
    }
    const vec3 translation = centroid / (-1 * nb_points_division_interface);

    //Move all the points so that the centroid is at the origin
    for(unsigned p_id = division_point_ids_threshold; p_id < m.node_pos_lst.size() / 3; p_id++){

        m.node_pos_lst[p_id *3    ] += translation.dx();
        m.node_pos_lst[p_id *3 + 1] += translation.dy();
        m.node_pos_lst[p_id *3 + 2] += translation.dz();
    }

    //Compute the rotation matrix that maps the plane to the xy plane
    const vec3 xy_plane_normal(0., 0., 1.);
    mat33 rotation_matrix;

    //If the division plane is already the xy plane, then the rotation matrix is the identity matrix
    if(xy_plane_normal.dot(division_plane_normal) == 1.0){rotation_matrix = mat33::identity();}

    //Otherwise we need to create a rotation matrix from a quaternion
    else{
        const vec3 a = division_plane_normal.cross(xy_plane_normal);
        const double w = 1.0 + division_plane_normal.dot(xy_plane_normal);

        //Use a quaternion to make sure that the rotation matrix is orthogonal
        quaternion q(w, a.dx(), a.dy(), a.dz());
        q = q.normalize();
        rotation_matrix = q.to_matrix();
    }


    //Rotate all the intesections points
    for(unsigned p_id = division_point_ids_threshold; p_id < m.node_pos_lst.size() / 3; p_id++){
        assert(p_id * 3 + 2 < m.node_pos_lst.size());

        vec3 pos_after_rotation = rotation_matrix.dot(vec3(
            m.node_pos_lst[p_id *3    ],
            m.node_pos_lst[p_id *3 + 1],
            m.node_pos_lst[p_id *3 + 2]
        ));

        m.node_pos_lst[p_id *3    ] = pos_after_rotation.dx();
        m.node_pos_lst[p_id *3 + 1] = pos_after_rotation.dy();
        m.node_pos_lst[p_id *3 + 2] = 0.;
    }
    


    return std::make_pair(translation, rotation_matrix);
}
//-----------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------
//Maps back the points of the division face to their original plane
void cell_divider::map_points_to_division_plane(
    mesh& m, 
    const unsigned division_point_ids_threshold,
    const vec3& translation,  
    const mat33& rotation_matrix
) noexcept{
    
    assert(division_point_ids_threshold < m.node_pos_lst.size() / 3);

    //Map back the nodes of the division face to their original positions
    const mat33 rotation_inv = rotation_matrix.transpose();

    for(unsigned p_id  = division_point_ids_threshold; p_id < m.node_pos_lst.size() / 3; p_id++){

        //Get the position of the point as a vector
        vec3 p(m.node_pos_lst[p_id * 3], m.node_pos_lst[p_id * 3 + 1], m.node_pos_lst[p_id * 3 + 2]);

        //Position back this point to the correct position of the division interface
        p = rotation_inv.dot(p) - translation;
        m.node_pos_lst[p_id * 3    ] = p.dx();
        m.node_pos_lst[p_id * 3 + 1] = p.dy();
        m.node_pos_lst[p_id * 3 + 2] = p.dz();
    }    
   
}
//-----------------------------------------------------------------------------------------------



//-----------------------------------------------------------------------------------------------
//Triangulate the interface between the 2 daughter cells in a way that respects the constraints on the edge lengths.
void cell_divider::triangulate_division_interface(
    const double l_min,                             //The minimum edge length
    mesh& m,                                        //The mesh containing all the faces and nodes
    const unsigned division_point_ids_threshold,    //All the nodes that belong to the division interface have an id greater than this threshold
    const unsigned division_face_ids_threshold,     //All the faces that belong to the division interface have an id greater than this threshold
    const vec3& division_face_normal               //The normal vector of the division interface
) noexcept(false){


    assert(l_min > 0.0);
    assert(m.node_pos_lst.size() % 3 == 0 && m.node_pos_lst.size() > 0);
    assert(m.face_point_ids.size() > 0);
    assert(division_point_ids_threshold < m.node_pos_lst.size() / 3);
    assert(division_face_ids_threshold < m.face_point_ids.size());
    
    //Compute the number of points  that are part of the division interface at this stage
    unsigned nb_points_division_interface = static_cast<unsigned>(m.node_pos_lst.size() / 3.) - division_point_ids_threshold;

    //Get a list of the coordinates of all the points located on the division interface
    std::vector<vec3> division_face_points_coords;

    //The last point is the centroid of the division interface and should not be added
    division_face_points_coords.reserve(nb_points_division_interface - 1); 

    //The last point at this stage is the centroid of the division interaface and should not be added
    for(unsigned p_id = division_point_ids_threshold; p_id < m.node_pos_lst.size() / 3 - 1; p_id++){
        division_face_points_coords.emplace_back(m.node_pos_lst[p_id * 3    ], m.node_pos_lst[p_id * 3 + 1], 0.);
    }

    //Use the poisson point cloud algorithm to triangulate the division interface
    std::vector<oriented_point> poisson_point_lst = poisson_sampling::compute_poisson_point_cloud(
        l_min, 
        m, 
        division_point_ids_threshold, 
        division_face_ids_threshold, 
        division_face_normal
    );
    
    //Store the points that have been created by the poisson sampling in the mesh
    for(const auto& poisson_point: poisson_point_lst){

        if(poisson_point.created_by_poisson_sampling_){
            //Add the point to the mesh
            m.node_pos_lst.insert(m.node_pos_lst.end(), {poisson_point.position_.dx(), poisson_point.position_.dy(), 0.});
        }
    }

    //Update the number of points  that are part of the division interface at this stage
    nb_points_division_interface = static_cast<unsigned>(m.node_pos_lst.size() / 3.) - division_point_ids_threshold;
  
    //Remove all the triangles of the division interface, they will be replaced by 
    //triangles respecting the constraints on the edge lengths
    m.face_point_ids.erase(m.face_point_ids.begin() + division_face_ids_threshold, m.face_point_ids.end());

    //Store in 2D the coordinates of the points of the division interface
    std::vector<double> division_interface_points_2d_coords;
    division_interface_points_2d_coords.reserve(nb_points_division_interface * 2);

    //Convert to 2D the positions of all the nodes located on the division interface
    for(unsigned p_id = division_point_ids_threshold; p_id < m.node_pos_lst.size() / 3; p_id++){
        division_interface_points_2d_coords.push_back(m.node_pos_lst[p_id * 3    ]);
        division_interface_points_2d_coords.push_back(m.node_pos_lst[p_id * 3 + 1]);
    }
    
    //Run the 2D delaunay triangulation
    try{
        delaunator::Delaunator d(division_interface_points_2d_coords);

        //Saves the faces created by the Delaunay algorithm into the mesh
        for(std::size_t i = 0; i < d.triangles.size(); i+=3) {
            assert(d.triangles[i + 2] < d.coords.size() / 2);

            //The Delaunay triangulation can create some triangles outside of the division interface. 
            //The following code checks if the triangle is inside the division interface
            const vec3 p0(d.coords[d.triangles[i    ] * 2], d.coords[d.triangles[i    ] * 2 + 1], 0.);
            const vec3 p1(d.coords[d.triangles[i + 1] * 2], d.coords[d.triangles[i + 1] * 2 + 1], 0.);
            const vec3 p2(d.coords[d.triangles[i + 2] * 2], d.coords[d.triangles[i + 2] * 2 + 1], 0.);

            //Compute the centroid of the triangle
            const vec3 face_centroid = (p0 + p1 + p2) / 3.;

            //If the triangle is in the division interface, then add it to the mesh
            if(cell_divider::point_is_in_polygon(division_face_points_coords, face_centroid)){

                assert(d.triangles[i]     < nb_points_division_interface);
                assert(d.triangles[i + 1] < nb_points_division_interface);
                assert(d.triangles[i + 2] < nb_points_division_interface);  

                //Get the ids of the 3 nodes in the mesh
                const unsigned n1_id =  division_point_ids_threshold + d.triangles[i    ];
                const unsigned n2_id =  division_point_ids_threshold + d.triangles[i + 1];
                const unsigned n3_id =  division_point_ids_threshold + d.triangles[i + 2];

                //Add the face to the mesh
                m.face_point_ids.insert(m.face_point_ids.end(),{n1_id, n2_id, n3_id});
            }
        }

    } catch (const std::exception& e) {
        throw division_exception("The Delaunay algorithm failed to triangulate the division interface.");
    }
    //catch (const std::runtime_error& error) {
    //    throw division_exception("The Delaunay algorithm failed to triangulate the division interface.");
    //}


}
//-----------------------------------------------------------------------------------------------


//-----------------------------------------------------------------------------------------------
//Check if a point is inside a polygon
bool cell_divider::point_is_in_polygon(
    const std::vector<vec3>& polygon, 
    const vec3& point
    ) noexcept{

    bool c = false;
    const double x = point.dx();
    const double y = point.dy();
    unsigned i, j;
    const unsigned nvert =  polygon.size();

    for (i = 0, j = nvert-1; i < nvert; j = i++) {
        const double vert_i_x = polygon[i].dx();
        const double vert_i_y = polygon[i].dy();
        const double vert_j_x = polygon[j].dx();
        const double vert_j_y = polygon[j].dy();

        //The point is part of the division face
        if((vert_i_x == x && vert_i_y == y)){return true;}
        if (((vert_i_y>y) != (vert_j_y>y)) &&
            (x < (vert_j_x-vert_i_x) * (y-vert_i_y) / (vert_j_y-vert_i_y) + vert_i_x)){
            c = !c;
        }
    }
  return c;
}
//-----------------------------------------------------------------------------------------------





//-----------------------------------------------------------------------------------------------




//-----------------------------------------------------------------------------------------------
//Split the mesh of the mother cell into 2 meshes, one for each daughter cell. And create the daughter cells
std::pair<cell_ptr, cell_ptr> cell_divider::create_daughter_cells(
    cell_ptr mother_cell, 
    mesh& mother_cell_mesh,                             //The mesh containing all the faces and nodes
    const unsigned division_point_ids_threshold,        //All the nodes that belong to the division interface have an id greater than this threshold
    const unsigned division_face_ids_threshold,         //All the faces that belong to the division interface have an id greater than this threshold
    const vec3& division_plane_normal,                  //The normal vector of the division interface
    const vec3& division_plane_origin                   //A point located on the division plane
) noexcept(false){


    assert(mother_cell->get_cell_type()!= nullptr);
    assert(division_point_ids_threshold < mother_cell_mesh.node_pos_lst.size() / 3);
    assert(division_face_ids_threshold < mother_cell_mesh.face_point_ids.size());

    //Create a mesh of the mother cell which does not contain the division interface
    mesh mother_cell_mesh_no_division_interface;
    mother_cell_mesh_no_division_interface.node_pos_lst = mother_cell_mesh.node_pos_lst;
    std::copy(mother_cell_mesh.face_point_ids.begin(), mother_cell_mesh.face_point_ids.begin() + division_face_ids_threshold, std::back_inserter(mother_cell_mesh_no_division_interface.face_point_ids));

    
    //Create the two daughter cells
    cell_ptr daughter_cell_1 = mother_cell->get_cell_same_type(mother_cell_mesh_no_division_interface);
    cell_ptr daughter_cell_2 = mother_cell->get_cell_same_type(mother_cell_mesh_no_division_interface);

    //Get a reference of the daughter cells node_lst
    auto& mother_cell_node_lst = mother_cell->node_lst_;
    auto& daughter_cell_1_node_lst = daughter_cell_1->node_lst_;
    auto& daughter_cell_2_node_lst = daughter_cell_2->node_lst_;

    //Map the attributes of the nodes of the mother cell to the nodes of the daughter cell
    for(unsigned i = 0; i < mother_cell_node_lst.size(); i++){
        daughter_cell_1_node_lst[i] = mother_cell_node_lst[i];
        daughter_cell_2_node_lst[i] = mother_cell_node_lst[i];
    }

    //For each daughter cell, we remove all the faces that are on a given side of the division plane
    std::vector<unsigned> faces_to_remove_1;
    std::vector<unsigned> faces_to_remove_2;

    for(unsigned f_id = 0; f_id < daughter_cell_1->face_lst_.size(); f_id++){
        //Get the side of the face with respect to the division plane
        const face& f =  daughter_cell_1->face_lst_[f_id];

        //If the face 1 is not deleted in the daughter cell 1, then it is deleted in the daughter cell 2
        if(cell_divider::face_side_wrt_plane(f, daughter_cell_1, division_plane_origin, division_plane_normal)){
            faces_to_remove_1.push_back(f_id);
        }
        else{
            faces_to_remove_2.push_back(f_id);
        }
    }
    assert(faces_to_remove_1.size() != 0 && faces_to_remove_2.size() != 0);
    
    //Remove the faces that are on the wrong side of the division plane
    remove_index(daughter_cell_1->face_lst_, faces_to_remove_1);
    remove_index(daughter_cell_2->face_lst_, faces_to_remove_2);

    //We add all the faces of the division interface
    for(unsigned f_id = division_face_ids_threshold; f_id < mother_cell_mesh.face_point_ids.size(); f_id++){
        
        //Transform the mesh face into a cell face
        auto& f_mesh = mother_cell_mesh.face_point_ids[f_id];
        assert(f_mesh.size() == 3); // Make sure the face is a triangle

        //Create the cell face
        face f1(f_mesh[0], f_mesh[2], f_mesh[1], f_id);
        face f2(f_mesh[0], f_mesh[1], f_mesh[2], f_id);


        //Add the face to the daughter cell 1
        daughter_cell_1->face_lst_.push_back(f1);

        //Add the face to the daughter cell 2
        daughter_cell_2->face_lst_.push_back(f2);
    }

    //Initialize the properties of the 2 daughter cells
    daughter_cell_1->initialize_cell_properties();  
    daughter_cell_2->initialize_cell_properties(); 

    //Check that the daughter cells have the same type as the mother cell
    assert(daughter_cell_1->get_cell_type()->global_type_id_ == mother_cell->get_cell_type()->global_type_id_);
    assert(daughter_cell_2->get_cell_type()->global_type_id_ == mother_cell->get_cell_type()->global_type_id_);

    return std::make_pair(daughter_cell_1, daughter_cell_2);
}
//-----------------------------------------------------------------------------------------------


