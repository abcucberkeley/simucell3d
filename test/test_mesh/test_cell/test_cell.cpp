#include "cell_tester.hpp"

/*

Test the method of the cell class

*/




//---------------------------------------------------------------------------------------------------------
//Create a dummy cell with the shape of a cube where all faces are triangles
cell_ptr generate_dummy_cell(){

    std::vector<double> cell_node_pos_lst{
        0,0,0,  1,0,0,  1,0,1,
        0,0,1,  0,1,0,  1,1,0,
        0,1,1,  1,1,1
    };



    std::vector<unsigned> cell_face_connectivity{
        0, 1, 3, //0
        2, 3, 1, //1
        0, 4, 1, //2
        5, 1, 4, //3
        0, 3, 4, //4
        6, 4, 3, //5
        1, 5, 2, //6
        7, 2, 5, //7
        5, 4, 7, //8
        6, 7, 4, //9
        3, 2, 6, //10
        7, 6, 2  //11
    };

    cell_ptr c0 = std::make_shared<cell>(cell_node_pos_lst, cell_face_connectivity, 0);
    return c0;
}
//---------------------------------------------------------------------------------------------------------




//---------------------------------------------------------------------------------------------------------
//Test all the constructors of the cell class
int cell_tester::deleted_constructors_test(){


    //Test the default constructor
    bool t1 = !std::is_default_constructible<cell>::value;

    //Test the copy constructor
    bool t2 = std::is_copy_constructible<cell>::value;

    //Test the move constructor
    bool t3 = !std::is_move_constructible<cell>::value;

    return !(t1 && t2 && t3);
}



int cell_tester::trivial_default_constructor_1_test(){

    //Test the trivial default constructor with node and face objects
    //-----------------------------------------------------------------------
    node n1(1), n2(2), n3(3);
    face f1(n1, n2, n3, 0), f2(n1, n2, n3, 1);

    std::vector<node> node_lst = {n1, n2, n3};
    std::vector<face> face_lst = {f1, f2};

    //Generate the cell
    cell_ptr c1 = std::make_shared<cell>(node_lst, face_lst, 0);

    //Make sure the node_lst and face_lst have the correct size
    bool t1 = c1->get_node_lst().size() == node_lst.size();
    bool t2 = c1->get_face_lst().size() == face_lst.size();

    //Get the cell node and face lst
    const auto& cell_node_lst = c1->get_node_lst();
    const auto& cell_face_lst = c1->get_face_lst();


    //Check the coordinates of the nodes
    bool t3 = true; 
    for(size_t i = 0; i < cell_node_lst.size(); i++){

        const auto& n = cell_node_lst[i];
        if(n.pos().dx() != node_lst[i].pos().dx()) t3 == false;
        if(n.pos().dy() != node_lst[i].pos().dy()) t3 == false;
        if(n.pos().dz() != node_lst[i].pos().dz()) t3 == false;
        if(t3 == false) break;
    }

    //Check that each face is connected with the correct nodes
    const auto& f1_cell = cell_face_lst[0];
    const auto& f2_cell = cell_face_lst[1];

    //Extract the node ids of each faces
    const auto& [n_id1, n_id2, n_id3] = f1_cell.get_node_ids();
    const auto& [n_id4, n_id5, n_id6] = f2_cell.get_node_ids();

    //Check that they are correct
    bool t4 = (n_id1 == 1 && n_id2 == 2 && n_id3 == 3);
    bool t5 = (n_id4 == 1 && n_id5 == 2 && n_id6 == 3);

    //Check the cell id
    bool t6 = c1->get_id() == 0;
    //-----------------------------------------------------------------------

    return !(t1 && t2 && t3 && t4 && t5 && t6);
}
//---------------------------------------------------------------------------------------------------------



//---------------------------------------------------------------------------------------------------------
//Check if the instantiation of a cell with its node positions and face node ids works
int cell_tester::trivial_default_constructor_2_test(){

    //Generate some node positions
    //                               ----n1----  ----n2---- ----n3----         
    std::vector<double> node_pos_lst{0., 0., 0., 1., 0., 0., 0., 1., 0.};


    //Generate some faces
    //                                 ---f1--- ---f2---
    std::vector<unsigned> face_node_lst{0, 1, 2, 0, 1, 2};
    
    //Generate the cell
    cell_ptr c1 = std::make_shared<cell>(node_pos_lst, face_node_lst, 0);


    //Make sure the node_lst and face_lst have the correct size
    bool t1 = c1->get_node_lst().size() == (size_t) node_pos_lst.size()  / 3;
    bool t2 = c1->get_face_lst().size() == (size_t) face_node_lst.size() / 3;

    //Get the cell node and face lst
    const auto& cell_node_lst = c1->get_node_lst();
    const auto& cell_face_lst = c1->get_face_lst();


    //Check the coordinates of the nodes
    bool t3 = true; 
    for(size_t i = 0; i < cell_node_lst.size(); i++){

        const auto& n = cell_node_lst[i];

        if(n.pos().dx() != node_pos_lst[i * 3 + 0]) t3 == false;
        if(n.pos().dy() != node_pos_lst[i * 3 + 1]) t3 == false;
        if(n.pos().dz() != node_pos_lst[i * 3 + 2]) t3 == false;

        if(t3 == false) break;
    }

    //Check that each face is connected with the correct nodes
    bool t4 = true; 
    for(size_t i = 0; i < cell_face_lst.size(); i++){

        //Extract the face 
        const auto& f = cell_face_lst[i];

        //Extract its node ids
        const auto& [n1, n2, n3] = f.get_node_ids();

        //Check that each face has the correct ids
        if(n1 != node_pos_lst[i * 3 + 0]) t4 == false;
        if(n2 != node_pos_lst[i * 3 + 1]) t4 == false;
        if(n3 != node_pos_lst[i * 3 + 2]) t4 == false;
        if(t4 == false) break;
    }

    //Check the id of the cell
    bool t5 = c1->get_id() == 0;


    return !(t1 && t2 && t3 && t4 && t5);
}
//---------------------------------------------------------------------------------------------------------


//---------------------------------------------------------------------------------------------------------
//Test if the cell can be instantiated from a mesh object
int cell_tester::instantiation_from_mesh_test(){


    //Create a mesh object with dummy faces and nodes
    mesh m;
    m.node_pos_lst = {        
        0,0,0,  1,0,0,  1,0,1,
        0,0,1,  0,1,0,  1,1,0,
        0,1,1,  1,1,1
    };

    m.face_point_ids = {        
        {0, 1, 3},
        {2, 3, 1},
        {0, 4, 1},
        {5, 1, 4},
        {3, 0, 4},
        {6, 4, 3},
        {5, 1, 2},
        {7, 2, 5},
        {5, 7, 4},
        {6, 7, 4},
        {3, 2, 6},
        {7, 6, 2} 
    };

    //Create a cell object from the mesh object
    cell_ptr c = std::make_shared<cell>(m, 0);

    //Get the cell node_lst and face_lst
    const auto cell_node_coord = c->get_flat_node_coord_lst();
    
    //Create a 2D vector of node indices based on the faces of the cell
    std::vector<unsigned> face_node_ids_lst;
    face_node_ids_lst.reserve(c->face_lst_.size());

    for(const face& f: c->face_lst_){
        const auto face_node_ids = f.get_node_ids();
        face_node_ids_lst.insert(face_node_ids_lst.end(), face_node_ids.begin(), face_node_ids.end());
    }

    const auto correct_face_id_lst = flatten(m.face_point_ids);

    //Compare the nodes and faces of the mesh and the cell
    bool t1 = std::equal(correct_face_id_lst.begin(), correct_face_id_lst.end(), face_node_ids_lst.begin());
    bool t2 = std::equal(cell_node_coord.begin(), cell_node_coord.end(), m.node_pos_lst.begin());

    return !(t1 && t2);
}
//---------------------------------------------------------------------------------------------------------


//---------------------------------------------------------------------------------------------------------
int cell_tester::all_nodes_are_used_test(){

    //Generate a cell where all the nodes are used
    std::vector<double> node_pos_lst_1{0., 0., 0., 1., 0., 0., 0., 1., 0.};

    //Generate some faces
    //                                 ---f1--- ---f2---
    std::vector<unsigned> face_node_lst_1{0, 1, 2, 0, 1, 2};
    
    //Generate the cell
    cell_ptr c1 = std::make_shared<cell>(node_pos_lst_1, face_node_lst_1, 0);
    c1->remove_unused_nodes();

    bool t1 = c1->free_node_queue_.size() == 0;

    //Generate a cell where one node is not used
    std::vector<double> node_pos_lst_2{0., 0., 0., 1., 0., 0., 0., 1., 0., 0., 1., 0.};

    //Generate some faces
    //                                 ---f1--- ---f2---
    std::vector<unsigned> face_node_lst_2{0, 2, 3, 0, 2, 3};
    
    //Generate the cell
    cell_ptr c2 = std::make_shared<cell>(node_pos_lst_2, face_node_lst_2, 0);
    c2->remove_unused_nodes();
    bool t2 = c2->free_node_queue_.size() == 1;

    std::cout << t1 << " " << t2 << std::endl;
    return !(t1 && t2);
}
//---------------------------------------------------------------------------------------------------------





//---------------------------------------------------------------------------------------------------------
//Test the generation of the edge_set
int cell_tester::edge_set_generation_test(){

    //Generate some arbitrary nodes and faces
    node n0(0), n1(1), n2(2), n3(3);
    face f1(n0, n1, n2, 0), f2(n0, n1, n3, 1);

    std::vector<node> node_lst = {n1, n2, n3};
    std::vector<face> face_lst = {f1, f2};

    //Generate the cell
    cell_ptr c1 = std::make_shared<cell>(node_lst, face_lst, 0);


    //Generate the edges of the cell
    c1->generate_edge_set();

    //Get the cell edge_Set
    const edge_set& cell_edge_set = c1->get_edge_set();

    //Check that the size is correct
    bool t1 = cell_edge_set.size() == 5; 

    //Find some of the edges generated by the cell
    auto it1  = std::find(cell_edge_set.begin(), cell_edge_set.end(), edge(0, 1));
    auto it2  = std::find(cell_edge_set.begin(), cell_edge_set.end(), edge(1, 2));
    auto it3  = std::find(cell_edge_set.begin(), cell_edge_set.end(), edge(1, 3));

    //Check that those edges exist
    bool t2 = it1 != cell_edge_set.end();
    bool t3 = it2 != cell_edge_set.end();
    bool t4 = it3 != cell_edge_set.end();


    //Check that the edge are correctly manifold or not
    bool t5 = false, t6 = false, t7 = false, t8 = false;
    
    if(t2) t5 = it1->is_manifold() == true;
    if(t3) t6 = it2->is_manifold() == false;
    if(t4) t7 = it3->is_manifold() == false;

    //Now check that the edge between the node 0 and 1 is correctly created
    if(t5){
        const auto [n_id1, n_id2] = it1->get_node_ids();
        t8 = (n_id1 == 0 && n_id2 == 1);
    }

    return !(t1 && t2 && t3 && t4 && t5 && t6 && t7 && t8);

}
//---------------------------------------------------------------------------------------------------------




//---------------------------------------------------------------------------------------------------------
//Test the rebase method of the cell class, on the node side. Removes 2 nodes from the 
//cell and check that the integrity of the cell mesh is conserved as expected.
int cell_tester::rebase_node_test(){

    //Generate a node with default coordinates 
    //The node 0 and 4 are not used by any faces and can be deleted
    node n0(0), n1(1), n2(2), n3(3), n4(4), n5(5), n6(6), n7(7);
    face f0(n1, n2, n3, 0), f1(n5, n6, n7, 1);
    std::vector<node> node_lst{n0, n1, n2, n3, n4, n5, n6, n7};
    std::vector<face> face_lst{f0, f1};

    //Create a dummy cell
    cell_ptr c1 = std::make_shared<cell>(node_lst, face_lst, 1);
    bool t1 = c1->node_lst_.size() == 8;


    //Reset the first and third node of the cell
    c1->node_lst_[0].reset();
    c1->node_lst_[4].reset();

    //Also add these two nodes to the queue of free nodes
    c1->add_free_node(0);
    c1->add_free_node(4);

    bool t2 = c1->free_node_queue_.size() == 2;

    //Call the function rebase to discard this nodes from the cell mesh
    c1->rebase();

    bool t3 = c1->free_node_queue_.size() == 0;

    //Check the ids of the remaining nodes
    bool t4 =   c1->node_lst_[0].get_local_id() == 0 && 
                c1->node_lst_[1].get_local_id() == 1 && 
                c1->node_lst_[2].get_local_id() == 2;

    //Check that the node ids in the faces have correctly been updated
    face f3= c1->face_lst_[0];
    face f4= c1->face_lst_[1];

    bool t5 =   f3.get_node_ids()[0] == 0 && 
                f3.get_node_ids()[1] == 1 && 
                f3.get_node_ids()[2] == 2;

    bool t6 =   f4.get_node_ids()[0] == 3 && 
                f4.get_node_ids()[1] == 4 && 
                f4.get_node_ids()[2] == 5;

    bool t7 = c1->node_lst_.size() == 6;

    return !(t1 && t2 && t3 && t4 && t5 && t6 && t7);
}
//---------------------------------------------------------------------------------------------------------



//---------------------------------------------------------------------------------------------------------
//Test the rebase method of the cell class, on the face side. Remove 1 face from the 
//cell and check that the integrity of the cell mesh is conserved as expected.
int cell_tester::rebase_face_test(){

    //Generate a node with default coordinates 
    //The node 0 and 4 are not used by any faces and can be deleted
    node n0(0), n1(1), n2(2), n3(3), n4(4), n5(5);
    face f0(n0, n1, n2, 0), f1(n3, n4, n5, 1);
    std::vector<node> node_lst{n0, n1, n2, n3, n4, n5};
    std::vector<face> face_lst{f0, f1};

    //Create a dummy cell
    cell_ptr c1 = std::make_shared<cell>(node_lst, face_lst, 1);
    bool t1 = c1->node_lst_.size() == 6;
    bool t2 = c1->face_lst_.size() == 2;


    //Reset the second face
    c1->face_lst_[1].reset();

    //Indicate to the cell that this face is free to be reused
    c1->add_free_face(1);

    bool t3 = c1->free_face_queue_.size() == 1;

    //Call the function rebase to discard this nodes from the cell mesh
    c1->rebase();

    bool t4 = c1->free_face_queue_.size() == 0;
    bool t5 = c1->face_lst_.size() == 1;

    //Check that the first face has been unaffected
    face f3= c1->face_lst_[0];
    bool t6 =   f3.get_node_ids()[0] == 0 && 
                f3.get_node_ids()[1] == 1 && 
                f3.get_node_ids()[2] == 2;


    return !(t1 && t2 && t3 && t4 && t5 && t6);
}
//---------------------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------------------
int cell_tester::get_flat_node_coord_lst_test(){

//The node 0 and 4 are not used by any faces and can be deleted
    node n0(1., 2., 3., 0), n1(3., 4., 5., 1), n2(6., 7., 8., 2);
    std::vector<node> node_lst{n0, n1, n2};
    std::vector<face> face_lst{};

    //Create a dummy cell
    cell_ptr c1 = std::make_shared<cell>(node_lst, face_lst, 1);

    std::vector<double> node_coord_lst = c1->get_flat_node_coord_lst();
    std::vector<double> correct_node_coord_lst{1., 2., 3., 3., 4., 5., 6., 7., 8.};
    bool t1 = std::equal(correct_node_coord_lst.begin(), correct_node_coord_lst.end(), node_coord_lst.begin());
    bool t2 = node_coord_lst.size() == 9;
    return !(t1 && t2);
}
//---------------------------------------------------------------------------------------------------------



//---------------------------------------------------------------------------------------------------------
int cell_tester::get_face_normal_test(){

    //Generate a node with default coordinates 
    //The node 0 and 4 are not used by any faces and can be deleted
    node n1(0., 0., 0., 0);
    node n2(1., 0., 0., 1);
    node n3(0., 1., 0., 2);

    //Create a face_lst and node_lst for the new cell
    face f0(n1, n2, n3, 0);
    std::vector<node> node_lst{n1, n2, n3};
    std::vector<face> face_lst{f0};

    //Create a dummy cell
    cell_ptr c1 = std::make_shared<cell>(node_lst, face_lst, 1);

    //Get the normal of the first face
    c1->update_face_normal_and_area(0);
    vec3 normal = c1->get_face_normal(0);
    const auto [n_x, n_y, n_z] = normal.to_array();

    bool t1 = almost_equal(n_x, 0.) && almost_equal(n_y, 0.) && almost_equal(n_z, 1.);
    return !t1; 
}
//---------------------------------------------------------------------------------------------------------




//---------------------------------------------------------------------------------------------------------
int cell_tester::check_face_normal_orientation_test(){

    //Generate a dmmy cell
    std::vector<double> cell_node_pos{
        0,0,0,  1,0,0,  1,0,1,
        0,0,1,  0,1,0,  1,1,0,
        0,1,1,  1,1,1
    };

    //Not all the normals are pointing outward
    std::vector<unsigned> cell_face_connectivity{
        0, 1, 3,
        2, 3, 1,
        0, 4, 1,
        5, 1, 4,
        3, 0, 4,
        6, 4, 3,
        5, 1, 2,
        7, 2, 5,
        5, 7, 4,
        6, 7, 4,
        3, 2, 6,
        7, 6, 2 
    };

    cell_ptr c0 = std::make_shared<cell>(cell_node_pos, cell_face_connectivity, 0);

    c0->generate_edge_set();
    c0->check_face_normal_orientation();
    c0->check_face_normal_orientation();

    c0->update_all_face_normals_and_areas();

    c0->area_ = c0->compute_area();

    //Since the cell geometry is convex
    vec3 cell_centroid = c0->compute_centroid();
    bool t1 = true;

    //Loop over the faces of the cell
    for(unsigned face_id = 0; face_id < c0->get_face_lst().size(); face_id++){

        //Get the normal of the face
        vec3 face_normal = c0->get_face_normal(face_id);

        //Get the nodes of the face
        const auto [n1_id, n2_id, n3_id] = c0->get_face_lst()[face_id].get_node_ids();
        assert(n1_id < c0->get_node_lst().size());
        const node& n1 = c0->get_node_lst()[n1_id];

        //Get a vector from the cell_centroid to the face
        const vec3 v = n1.pos_ - cell_centroid;

        if(v.dot(face_normal) < 0){
            t1 = false; 
        }
    }

    return !(t1);
}
//---------------------------------------------------------------------------------------------------------



//---------------------------------------------------------------------------------------------------------

int cell_tester::compute_centroid_test(){

     std::vector<double> cell_node_pos_1{
        0,0,0,  1,0,0,  1,0,1,
        0,0,1,  0,1,0,  1,1,0,
        0,1,1,  1,1,1,  1,1,6
    };


    std::vector<unsigned> cell_face_connectivity{
        0, 1, 3,
        2, 3, 1,
        0, 4, 1,
        5, 1, 4,
        0, 3, 4,
        6, 4, 3,
        1, 5, 2,
        7, 2, 5,
        5, 4, 7,
        6, 7, 4,
        3, 2, 6,
        7, 6, 2 
    };


    cell_ptr c0 = std::make_shared<cell>(cell_node_pos_1, cell_face_connectivity, 0);


    c0->initialize_cell_properties(false);


    //Since the cell geometry is convex
    vec3 cell_centroid = c0->compute_centroid();

    cell_centroid.print();

    bool t1 = cell_centroid == vec3(0.5, 0.5, 0.5);

    return !(t1);
}
//---------------------------------------------------------------------------------------------------------




//---------------------------------------------------------------------------------------------------------
int cell_tester::compute_volume_test(){
     
     std::vector<double> cell_node_pos_1{
        0,0,0,  1,0,0,  1,0,1,
        0,0,1,  0,1,0,  1,1,0,
        0,1,1,  1,1,1
    };


    std::vector<unsigned> cell_face_connectivity{
        0, 1, 3,
        2, 3, 1,
        0, 4, 1,
        5, 1, 4,
        0, 3, 4,
        6, 4, 3,
        1, 5, 2,
        7, 2, 5,
        5, 4, 7,
        6, 7, 4,
        3, 2, 6,
        7, 6, 2 
    };

    cell_ptr c0 = std::make_shared<cell>(cell_node_pos_1, cell_face_connectivity, 0);
    double cell_volume = c0->compute_volume();

    return !(cell_volume == 1.);
}
//---------------------------------------------------------------------------------------------------------


//---------------------------------------------------------------------------------------------------------
//Detects if there are holes in the input surface
int cell_tester::is_manifold_test(){
     //Generate a dmmy cell
    std::vector<double> cell_node_pos{
        0,0,0,  1,0,0,  1,0,1,
        0,0,1,  0,1,0,  1,1,0,
        0,1,1,  1,1,1
    };

    //Not all the normals are pointing outward
    std::vector<unsigned> cell_face_connectivity{
        0, 1, 3,
        2, 3, 1,
        0, 4, 1,
        5, 1, 4,
        3, 0, 4,
        6, 4, 3,
        5, 1, 2,
        7, 2, 5,
        5, 7, 4,
        6, 7, 4,
        3, 2, 6,
        7, 6, 2 
    };


    //This cell has no holes
    cell_ptr c0 = std::make_shared<cell>(cell_node_pos, cell_face_connectivity, 0);
    c0->generate_edge_set();
    bool t1 = c0->is_manifold() == true;



    //Not all the normals are pointing outward
    std::vector<unsigned> cell_face_connectivity_2{
        0, 1, 3,
        2, 3, 1,
        5, 1, 4,
        3, 0, 4,
        6, 4, 3,
        5, 1, 2,
        7, 2, 5,
        5, 7, 4,
        6, 7, 4,
        3, 2, 6,
        7, 6, 2 
    };


    //This cell has a hole in its surface
    cell_ptr c1 = std::make_shared<cell>(cell_node_pos, cell_face_connectivity_2, 0);
    c1->generate_edge_set();
    
    bool t2 = c1->is_manifold() == false;

    return !(t1 && t2);
}
//---------------------------------------------------------------------------------------------------------



//---------------------------------------------------------------------------------------------------------
int cell_tester::update_face_area_test(){
    //Get face area 


    //Generate some node positions
    //                               ----n1----  ----n2---- ----n3----         
    std::vector<double> node_pos_lst{0., 0., 0., 1., 0., 0., 0., 1., 0.};


    //Generate some faces
    //                                 ---f1--- 
    std::vector<unsigned> face_node_lst{0, 1, 2};
    
    //Generate the cell
    cell_ptr c1 = std::make_shared<cell>(node_pos_lst, face_node_lst, 0);


    //Update the area of the face
    c1->update_all_face_normals_and_areas();

    //Get the area of the cell face
    const double area = c1->get_face_lst()[0].get_area();
    return !(almost_equal(area, 0.5));
}
//---------------------------------------------------------------------------------------------------------


//---------------------------------------------------------------------------------------------------------
int cell_tester::get_aabb_test(){
    
    //Generate a dummy cell
    std::vector<double> node_pos_lst{
        0,0,0,  1,0,0,  1,0,1,
        0,0,1,  0,1,0,  1,1,0,
        0,1,1,  1,1,1
    };

    //Not all the normals are pointing outward
    std::vector<unsigned> face_node_lst{
        0, 1, 3,
        2, 3, 1,
        0, 4, 1,
        5, 1, 4,
        3, 0, 4,
        6, 4, 3,
        5, 1, 2,
        7, 2, 5,
        5, 7, 4,
        6, 7, 4,
        3, 2, 6,
        7, 6, 2 
    };

    //Generate the cell
    cell_ptr c1 = std::make_shared<cell>(node_pos_lst, face_node_lst, 0);

    //Get the face AABB
    std::array<double, 6> cell_aabb = c1->get_aabb();
    std::array<double, 6> correct_aabb{0., 0., 0., 1., 1., 1.};
    bool t1 = std::equal(cell_aabb.begin(), cell_aabb.end(), correct_aabb.begin());

    return !(t1);
}
//---------------------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------------------
int cell_tester::compute_area_test(){
    
    //Generate a dummy cell
    std::vector<double> node_pos_lst{
        0,0,0,  1,0,0,  1,0,1,
        0,0,1,  0,1,0,  1,1,0,
        0,1,1,  1,1,1,  2,2,2,
        3,3,3,  4,4,4
    };

    //Not all the normals are pointing outward
    std::vector<unsigned> face_node_lst{
        0, 1, 3, 
        2, 3, 1,
        0, 4, 1,
        5, 1, 4,
        3, 0, 4,
        6, 4, 3,
        5, 1, 2,
        7, 2, 5,
        5, 7, 4,
        6, 7, 4,
        3, 2, 6,
        7, 6, 2, 
        8, 9, 10
    };

    //Generate the cell
    cell_ptr c1 = std::make_shared<cell>(node_pos_lst, face_node_lst, 0);


    //Get the last face of the cell and reset it
    c1->face_lst_.back().reset();

    //Remove a face from the cell to spice things up
    c1->add_free_face(c1->face_lst_.back().get_local_id());

    //Compute the area of the cell (which is a cube)
    c1->update_all_face_normals_and_areas();

    const double cell_area = c1->compute_area();
    return !(almost_equal(cell_area, 6.));
}
//---------------------------------------------------------------------------------------------------------


//---------------------------------------------------------------------------------------------------------
//Unit test of the get_edge_test function
int cell_tester::get_edge_test(){

    //Get a cell with the shape of a cube where all the faces are triangulated
    cell_ptr c = generate_dummy_cell();

    //Generate the edges of the cell
    c->generate_edge_set();

    //Get an existing egde
    auto edge_opt = c->get_edge(0, 1);

    bool t1 = edge_opt.has_value();

    //Get a non-existing edge
    edge_opt = c->get_edge(0, 2);

    bool t2 = !edge_opt.has_value();


    return !(t1 && t2);
}
//---------------------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------------------
int cell_tester::delete_face_test(){

    cell_ptr c = generate_dummy_cell();

    //Generate the edges of the cell
    c->generate_edge_set();

    //Get the face 0
    face& f = c->face_lst_[0];

    bool t1 = f.is_used() == true;

    //Get the edges of the face 0
    const auto [n1_id, n2_id, n3_id] = f.get_node_ids();

    //Delete the face         node_ids.push_back(node_id == e.n1() ? e.n2() : e.n1());
    c->delete_face(f.get_local_id());

    //Make sure the face has been marked as non used
    bool t2 = f.is_used() == false;

    //Make sure the spot 0 of the face lst has been marked as free
    bool t3 = c->free_face_queue_.size() == 1;
    bool t4 = (t3) ? c->free_face_queue_.front() == 0 : false;

    //Make sure the edges of the face 0 has been updated of the face deletion
    auto edge_opt_1 = c->get_edge(n1_id, n2_id);
    auto edge_opt_2 = c->get_edge(n2_id, n3_id);
    auto edge_opt_3 = c->get_edge(n3_id, n1_id);
    if (!edge_opt_1.has_value() || !edge_opt_2.has_value() || !edge_opt_3.has_value()) return 1;

    //get_edge returns a copy of the edge and not a reference
    edge edge_1 = edge_opt_1.value();
    edge edge_2 = edge_opt_2.value();
    edge edge_3 = edge_opt_3.value();

    bool t5 = edge_1.is_manifold() == false;
    bool t6 = edge_2.is_manifold() == false;
    bool t7 = edge_3.is_manifold() == false;

    //Now delete a second face of the mesh
    c->delete_face(2);

    //The edge between node 0 and 1 should be have been deleted
    auto edge_opt_4 = c->get_edge(0, 1);

    bool t8 = !edge_opt_4.has_value();

    std::cout << "t1: " << t1 << std::endl;
    std::cout << "t2: " << t2 << std::endl;
    std::cout << "t3: " << t3 << std::endl;
    std::cout << "t4: " << t4 << std::endl;
    std::cout << "t5: " << t5 << std::endl;
    std::cout << "t6: " << t6 << std::endl;
    std::cout << "t7: " << t7 << std::endl;
    std::cout << "t8: " << t8 << std::endl;


    return !(t1 && t2 && t3 && t4 && t5 && t6 && t7 && t8);
    
}
//---------------------------------------------------------------------------------------------------------


//---------------------------------------------------------------------------------------------------------
int cell_tester::add_face_test_1(){
    
    cell_ptr c = generate_dummy_cell();

    c->generate_edge_set();

    std::cout << "face_lst size: " << c->face_lst_.size() << std::endl;


    //Get the face 1
    face& f = c->face_lst_[1];

 
    //Get the node ids of the face 1
    const auto [n1_id, n2_id, n3_id] = f.get_node_ids();


    //Try to create a face with the same node ids, we should get an error
    try{
        c->create_face(n1_id, n2_id, n3_id);

        std::cout << "no error" << std::endl;
        return 1;
    }
    catch (const std::exception& e){}

    face& f2 = c->face_lst_[1];

    //Delete the face 1
    c->delete_face(f2);

    //Create a new face with the same node ids, this time we should not get an error
    c->create_face(n1_id, n2_id, n3_id);

    //Make sure the face 1 has been marked as used
    bool t1 = c->face_lst_[1].is_used() == true;
    bool t2 = c->free_face_queue_.size() == 0;
    bool t3 = c->face_lst_[1].n1_id_ == n1_id;
    bool t4 = c->face_lst_[1].n2_id_ == n2_id;
    bool t5 = c->face_lst_[1].n3_id_ == n3_id;
    bool t6 = c->face_lst_[1].get_local_id() == 1;

    //Make sure the new face 0 has been added to its edges
    auto edge_opt_1 = c->get_edge(n1_id, n2_id);
    auto edge_opt_2 = c->get_edge(n2_id, n3_id);
    auto edge_opt_3 = c->get_edge(n3_id, n1_id);
    if (!edge_opt_1.has_value() || !edge_opt_2.has_value() || !edge_opt_3.has_value()) return 1;

    edge edge_1 = edge_opt_1.value();
    edge edge_2 = edge_opt_2.value();
    edge edge_3 = edge_opt_3.value();

    if(!edge_1.is_manifold() || !edge_2.is_manifold() || !edge_3.is_manifold()) return 1;

    bool t7 = (edge_1.f1() == 1) || (edge_1.f2() == 1);
    bool t8 = (edge_2.f1() == 1) || (edge_2.f2() == 1);
    bool t9 = (edge_3.f1() == 1) || (edge_3.f2() == 1);

    std::cout << "t1: " << t1 << std::endl;
    std::cout << "t2: " << t2 << std::endl;
    std::cout << "t3: " << t3 << std::endl;
    std::cout << "t4: " << t4 << std::endl;
    std::cout << "t5: " << t5 << std::endl;
    std::cout << "t6: " << t6 << std::endl;
    std::cout << "t7: " << t7 << std::endl;
    std::cout << "t8: " << t8 << std::endl;
    std::cout << "t9: " << t9 << std::endl;
    return !(t1 && t2 && t3 && t4 && t5 && t6 && t7 && t8 && t9);
}
//---------------------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------------------
int cell_tester::add_face_test_2(){

    cell_ptr c = generate_dummy_cell();
    c->generate_edge_set();

    //Create 3 new nodes
    node n1(0., 0., 0., 0);
    node n2(1., 0., 0., 0);
    node n3(0., 1., 0., 0);

    c->add_node(n1);
    node& n1_ = c->node_lst_.back();
    
    c->add_node(n2);
    node& n2_ = c->node_lst_.back();

    c->add_node(n3);
    node& n3_ = c->node_lst_.back();

    //Create a new face with the new nodes
    c->create_face(n1_.get_local_id(), n2_.get_local_id(), n3_.get_local_id());

    //Make sure the created face is at the end of the face list
    face& f = c->face_lst_.back();
    bool t1 = f.n1_id_ == n1_.get_local_id();
    bool t2 = f.n2_id_ == n2_.get_local_id();
    bool t3 = f.n3_id_ == n3_.get_local_id();

    //Make sure the index of the face is correct
    bool t4 = f.get_local_id() == c->face_lst_.size() - 1;

    //Make sure new edges have been created
    auto edge_opt_1 = c->get_edge(n1_.get_local_id(), n2_.get_local_id());
    auto edge_opt_2 = c->get_edge(n2_.get_local_id(), n3_.get_local_id());
    auto edge_opt_3 = c->get_edge(n3_.get_local_id(), n1_.get_local_id());
    if (!edge_opt_1.has_value() || !edge_opt_2.has_value() || !edge_opt_3.has_value()) return 1;

    edge edge_1 = edge_opt_1.value();
    edge edge_2 = edge_opt_2.value();
    edge edge_3 = edge_opt_3.value();
    bool t5 = !edge_1.is_manifold() && !edge_2.is_manifold() && !edge_3.is_manifold();


    bool t6 = (edge_1.f1() == f.get_local_id()) || (edge_1.f2() == f.get_local_id());
    bool t7 = (edge_2.f1() == f.get_local_id()) || (edge_2.f2() == f.get_local_id());
    bool t8 = (edge_3.f1() == f.get_local_id()) || (edge_3.f2() == f.get_local_id());

    std::cout << "t1: " << t1 << std::endl;
    std::cout << "t2: " << t2 << std::endl;
    std::cout << "t3: " << t3 << std::endl;
    std::cout << "t4: " << t4 << std::endl;
    std::cout << "t5: " << t5 << std::endl;
    std::cout << "t6: " << t6 << std::endl;
    std::cout << "t7: " << t7 << std::endl;
    std::cout << "t8: " << t8 << std::endl;

    return !(t1 && t2 && t3 && t4 && t5 && t6 && t7 && t8);
    
}
//---------------------------------------------------------------------------------------------------------


//---------------------------------------------------------------------------------------------------------
int cell_tester::add_node_test(){

    cell_ptr c = generate_dummy_cell();

    //Create a new node
    node n(0., 0., 0., 0);

    c->add_node(n);

    //Make sure the created node is at the end of the node list
    node& n_ = c->node_lst_.back();
    bool t1 = n_.get_local_id() == c->node_lst_.size() - 1;
    bool t2 = n_.pos() == n.pos();

    //Remove a node from the node_lst_
    c->delete_node(1);
    node n2(1., 1., 1., 1);
    c->add_node(n2);
    bool t3 = c->node_lst_[1].pos() == n2.pos();

    std::cout << "t1: " << t1 << std::endl;
    std::cout << "t2: " << t2 << std::endl;
    std::cout << "t3: " << t3 << std::endl;

    return !(t1 && t2 && t3);
}
//---------------------------------------------------------------------------------------------------------


//---------------------------------------------------------------------------------------------------------
int cell_tester::delete_node_test(){

    cell_ptr c = generate_dummy_cell();

    //Delete the first node of the cell
    node& n = c->node_lst_[0];

    bool t1 = n.is_used() == true;

    c->delete_node(0);

    bool t2 = n.is_used() == false;
    bool t3 = (c->free_node_queue_.size() == 1) ? c->free_node_queue_.front() == 0 : false;

    std::cout << "t1: " << t1 << std::endl;
    std::cout << "t2: " << t2 << std::endl;
    std::cout << "t3: " << t3 << std::endl;

    return !(t1 && t2 && t3);
}
//---------------------------------------------------------------------------------------------------------


//---------------------------------------------------------------------------------------------------------

int cell_tester::replace_node_test(){

    cell_ptr c = generate_dummy_cell();
    c->initialize_cell_properties();

    //Create a new node
    const unsigned new_node_id = c->create_node(c->get_node_lst()[0].pos());
    const double cell_area_before = c->get_area();
    const double cell_volume_before = c->get_volume();

    auto edge_opt = c->get_edge(0, 1);
    if (!edge_opt.has_value()) return 1;
    edge& e = const_cast<edge&>(edge_opt.value());

    //Replace the first node of the cell
    c->replace_node(e, 0, new_node_id);

    //Make sure that no edges still contain the node 0 in the edge set
    bool t1 = std::all_of(c->edge_set_.begin(), c->edge_set_.end(), [](const edge& e){return e.n1() != 0 && e.n2() != 0;});

    //Make sure that no face in the cell use the node 0
    bool t2 = std::all_of(c->face_lst_.begin(), c->face_lst_.end(), [](const face& f){
        auto [n1, n2, n3] = f.get_node_ids();
        return n1 != 0 && n2 != 0 && n3 != 0;
    });

    //Make sure the edge 0 has been marked as unused
    bool t3 = !c->node_lst_[0].is_used();

    //Make sure the new node has been marked as used
    bool t4 = c->node_lst_[new_node_id].is_used();

    //Make sure the cell properties have been unchanged since the positions of the nodes have not changed
    const double cell_area_after =   c->compute_area();
    const double cell_volume_after = c->compute_volume();

    bool t5 = almost_equal(cell_area_before, cell_area_after);
    bool t6 = almost_equal(cell_volume_before, cell_volume_after);

    std::cout << "t1 : " << t1 << std::endl;
    std::cout << "t2 : " << t2 << std::endl;
    std::cout << "t3 : " << t3 << std::endl;
    std::cout << "t4 : " << t4 << std::endl;
    std::cout << "t5 : " << t5 << std::endl;
    std::cout << "t6 : " << t6 << std::endl;

    return !(t1 && t2 && t3 && t4 && t5 && t6);



}
//---------------------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------------------
int cell_tester::copy_constructor_test(){

    //Check that the default copy constructor of c++ does the job as expected
    cell_ptr c1 = generate_dummy_cell();
    c1->initialize_cell_properties();

    //Set some of the cell atributes
    c1->free_face_queue_.push_back(0);
    c1->free_face_queue_.push_back(1);

    //Test the copy assignment operator
    cell_ptr c2 = std::make_shared<cell>(*c1);

    bool t1 = c1->get_id() == c2->get_id();
    bool t2 = std::equal(c1->get_node_lst().begin(), c1->get_node_lst().end(), c2->get_node_lst().begin());
    bool t3 = std::equal(c1->get_face_lst().begin(), c1->get_face_lst().end(), c2->get_face_lst().begin());
    bool t4 = std::equal(c1->get_edge_set().begin(), c1->get_edge_set().end(), c2->get_edge_set().begin());
    bool t5 = std::equal(c1->free_face_queue_.begin(), c1->free_face_queue_.end(), c2->free_face_queue_.begin());
    bool t6 = std::equal(c1->free_node_queue_.begin(), c1->free_node_queue_.end(), c2->free_node_queue_.begin());

    bool t7 = c1->get_area() == c2->get_area();
    bool t8 = c1->get_volume() == c2->get_volume();

    std::cout << "t1 : " << t1 << std::endl;
    std::cout << "t2 : " << t2 << std::endl;
    std::cout << "t3 : " << t3 << std::endl;
    std::cout << "t4 : " << t4 << std::endl;
    std::cout << "t5 : " << t5 << std::endl;
    std::cout << "t6 : " << t6 << std::endl;
    std::cout << "t7 : " << t7 << std::endl;
    std::cout << "t8 : " << t8 << std::endl;

    return !(t1 && t2 && t3 && t4 && t5 && t6 && t7 && t8);
}
//---------------------------------------------------------------------------------------------------------



//---------------------------------------------------------------------------------------------------------
int cell_tester::get_node_faces_test() {

    cell_ptr c = generate_dummy_cell();
    c->generate_edge_set();


    //Get the edge 0,1 
    auto edge_opt = c->get_edge(0, 1);
    if (!edge_opt.has_value()){

        std::cout << "Edge 0,1 not found" << std::endl;
        return 1;
    }
    edge& e = const_cast<edge&>(edge_opt.value());

    //Get all the faces of the node 0
    std::vector<unsigned> node_faces_0 = c->get_connected_nodes(0, e);
    std::sort(node_faces_0.begin(), node_faces_0.end());
    bool t1 = std::equal(node_faces_0.begin(), node_faces_0.end(), std::vector<unsigned>{1, 3, 4}.begin());


    //Get all the faces of the node 1
    std::vector<unsigned> node_faces_1 = c->get_connected_nodes(1, e);
    std::sort(node_faces_1.begin(), node_faces_1.end());
    bool t2 = std::equal(node_faces_1.begin(), node_faces_1.end(), std::vector<unsigned>{0, 2, 3, 4, 5}.begin());

    std::cout << "node_nodes_0 : ";
    print_container(node_faces_0);

    std::cout << "node_nodes_1 : ";
    print_container(node_faces_1);


    std::cout << t1 << std::endl;
    std::cout << t2 << std::endl;

    return !(t1 && t2);
}




//---------------------------------------------------------------------------------------------------------
int cell_tester::get_cell_longest_axis_test() {

    
    std::vector<double> cell_node_pos_lst{
        0,0,0,  1,0,0,  1,0,2,
        0,0,2,  0,1,0,  1,1,0,
        0,1,2,  1,1,2
    };



    std::vector<unsigned> cell_face_connectivity{
        0, 1, 3, //0
        2, 3, 1, //1
        0, 4, 1, //2
        5, 1, 4, //3
        0, 3, 4, //4
        6, 4, 3, //5
        1, 5, 2, //6
        7, 2, 5, //7
        5, 4, 7, //8
        6, 7, 4, //9
        3, 2, 6, //10
        7, 6, 2  //11
    };

    cell_ptr c0 = std::make_shared<cell>(cell_node_pos_lst, cell_face_connectivity, 0);

    //print the cell 
    mesh_writer::write_cell_data_file("./get_cell_longest_axis_test_mesh.vtk",  {c0});

    c0->update_all_face_normals_and_areas();
    c0->area_ = c0->compute_area();

    //Get the longest axis of the celll
    const vec3 longest_axis = c0->get_cell_longest_axis();
    longest_axis.print();

    return !(longest_axis == vec3(0., 0., 1.) || longest_axis == vec3(0., 0., -1.));
}
//---------------------------------------------------------------------------------------------------------



//---------------------------------------------------------------------------------------------------------
int cell_tester::get_angle_gradient_test(){

    const vec3 i(1., 5., 1.);
    const vec3 j(2., 1., 2.);
    const vec3 k(1., 2., 2.);


    //Those are the ground truth values directly calculated in mathematica
    const double grad_theta_ix = 0.08375315127160103;
    const double grad_theta_iy =-0.0636523949664169;
    const double grad_theta_iz =-0.006700252101728042;

    const double grad_theta_jx = 0.2177581933061626;
    const double grad_theta_jy = 0.033501260508640544;
    const double grad_theta_jz =-0.08375315127160103;

    const double grad_theta_kx =-0.30151134457776363;
    const double grad_theta_ky = 0.030151134457776358;
    const double grad_theta_kz = 0.09045340337332908;

    //Calculate the gradient of the angle
    const auto [grad_theta_i, grad_theta_j, grad_theta_k] = cell::get_angle_gradient(i, j, k);

    bool t1 = almost_equal(grad_theta_i.dx(), grad_theta_ix);
    bool t2 = almost_equal(grad_theta_i.dy(), grad_theta_iy);
    bool t3 = almost_equal(grad_theta_i.dz(), grad_theta_iz);

    bool t4 = almost_equal(grad_theta_j.dx(), grad_theta_jx);
    bool t5 = almost_equal(grad_theta_j.dy(), grad_theta_jy);
    bool t6 = almost_equal(grad_theta_j.dz(), grad_theta_jz);

    bool t7 = almost_equal(grad_theta_k.dx(), grad_theta_kx);
    bool t8 = almost_equal(grad_theta_k.dy(), grad_theta_ky);
    bool t9 = almost_equal(grad_theta_k.dz(), grad_theta_kz);

    std::cout<< t1 << std::endl;
    std::cout<< t2 << std::endl;
    std::cout<< t3 << std::endl;
    std::cout<< t4 << std::endl;
    std::cout<< t5 << std::endl;
    std::cout<< t6 << std::endl;
    std::cout<< t7 << std::endl;
    std::cout<< t8 << std::endl;
    std::cout<< t9 << std::endl;

    return !(t1 && t2 && t3 && t4 && t5 && t6 && t7 && t8 && t9);
}
//---------------------------------------------------------------------------------------------------------











//---------------------------------------------------------------------------------------------------------
// The main function
int main (int argc, char** argv){

    //Check that the command line input is correctly formatted
    assert(argc == 2); 

    //Get the name of the test to run
    std::string test_name = argv[1];
    
    //Run the selected test
    cell_tester tester;


    if (test_name == "deleted_constructors_test")               return tester.deleted_constructors_test();
    if (test_name == "trivial_default_constructor_1_test")      return tester.trivial_default_constructor_1_test();
    if (test_name == "trivial_default_constructor_2_test")      return tester.trivial_default_constructor_2_test();
    if (test_name == "edge_set_generation_test")                return tester.edge_set_generation_test();
    if (test_name == "rebase_node_test")                        return tester.rebase_node_test();
    if (test_name == "rebase_face_test")                        return tester.rebase_face_test();
    if (test_name == "get_flat_node_coord_lst_test")                 return tester.get_flat_node_coord_lst_test();
    if (test_name == "get_face_normal_test")                    return tester.get_face_normal_test();
    if (test_name == "check_face_normal_orientation_test")      return tester.check_face_normal_orientation_test();
    if (test_name == "compute_volume_test")                     return tester.compute_volume_test();
    if (test_name == "is_manifold_test")                        return tester.is_manifold_test();
    if (test_name == "compute_centroid_test")                   return tester.compute_centroid_test();
    if (test_name == "update_face_area_test")                   return tester.update_face_area_test();
    if (test_name == "get_aabb_test")                           return tester.get_aabb_test();
    if (test_name == "compute_area_test")                       return tester.compute_area_test();
    if (test_name == "instantiation_from_mesh_test")            return tester.instantiation_from_mesh_test();
    if (test_name == "all_nodes_are_used_test")                 return tester.all_nodes_are_used_test();
    if (test_name == "get_edge_test")                           return tester.get_edge_test();
    if (test_name == "delete_face_test")                        return tester.delete_face_test();
    if (test_name == "add_face_test_1")                         return tester.add_face_test_1();
    if (test_name == "add_face_test_2")                         return tester.add_face_test_2();
    if (test_name == "add_node_test")                           return tester.add_node_test();
    if (test_name == "delete_node_test")                        return tester.delete_node_test();
    if (test_name == "replace_node_test")                       return tester.replace_node_test();
    if (test_name == "copy_constructor_test")                   return tester.copy_constructor_test();
    if (test_name == "get_node_faces_test")                     return tester.get_node_faces_test();
    if (test_name == "get_cell_longest_axis_test")              return tester.get_cell_longest_axis_test();
    if (test_name == "get_angle_gradient_test")                 return tester.get_angle_gradient_test();




    std::cout << "TEST NAME :" << test_name << " DOES NOT EXIST" << std::endl;
    return 1;

}
//---------------------------------------------------------------------------------------------------------
