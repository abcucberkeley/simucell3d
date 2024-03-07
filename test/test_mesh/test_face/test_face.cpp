#include "face_tester.hpp"


/*

Test the method of the node class

*/


typedef std::shared_ptr<cell> cell_ptr;



//---------------------------------------------------------------------------------------------------------
//Test the default instantiation of a face object
int face_tester::default_instantiation_test(){

    //Check that the default instantiation is not possible
    return std::is_default_constructible<face>::value;
}
//---------------------------------------------------------------------------------------------------------



//---------------------------------------------------------------------------------------------------------
//Test the default constructors of the face class
int face_tester::default_constructors_test(){

    //Check the default copy constructor

    //Create a dummy face
    face f1(1, 2, 3, 0);

    //This copy has the same cell_ptr but is not registered in the cell as one of the faces
    face f2(f1);
    auto [n_id1, n_id2, n_id3 ] = f2.get_node_ids();
    bool t2 = n_id1 == 1 &&  n_id2 == 2 && n_id3 == 3;
    bool t3 = f1.get_local_id() == f2.get_local_id();

    //Do the same test with the copy assignment function
    face f3 =  f1;
    auto [n_id4, n_id5, n_id6 ] = f3.get_node_ids();
    bool t4 = n_id4 == 1 &&  n_id5 == 2 && n_id6 == 3;
    bool t5 = f1.get_local_id() == f3.get_local_id();

    //Test the move assignment operator
    face f4 = face(1, 2, 3, 0);
    auto [n_id7, n_id8, n_id9 ] = f4.get_node_ids();
    bool t6 = n_id7 == 1 &&  n_id8 == 2 && n_id9 == 3;
    bool t7 = 0 == f4.get_local_id();

    return !(t2 && t3 && t4 && t5 && t6 && t7);
}
//---------------------------------------------------------------------------------------------------------



//---------------------------------------------------------------------------------------------------------
//Test the constructor of the face class
int face_tester::copy_constructor_test(){

    //Create arbitrary nodes
    node n1(1), n2(2), n3(3);

    //Create a face
    face f1(n1, n2, n3, 0);

    //Check the id of the face
    bool t1 = f1.get_local_id() == 0;

    //Chekc that t1 is connected to the correct nodes
    const auto& [f1_n1_id, f1_n2_id, f1_n3_id] = f1.get_node_ids();
    bool t2 = (f1_n1_id == 1 && f1_n2_id == 2 && f1_n3_id == 3);

    return !(t1 && t2);
}
//---------------------------------------------------------------------------------------------------------



//---------------------------------------------------------------------------------------------------------
//Test the constructor of where the ids of the nodes are directly given to the face 
int face_tester::node_id_constructor_test(){

    //Create a face
    face f1(1, 2, 3, 0);

    //Check the id of the face
    bool t1 = f1.get_local_id() == 0;

    //Chekc that t1 is connected to the coorect nodes
    const auto& [f1_n1_id, f1_n2_id, f1_n3_id] = f1.get_node_ids();
    bool t2 = (f1_n1_id == 1 && f1_n2_id == 2 && f1_n3_id == 3);


    return !(t1 && t2);
}
//---------------------------------------------------------------------------------------------------------



//---------------------------------------------------------------------------------------------------------
//Test that the id of a face can correctly be reset 
int face_tester::reset_local_id_test(){

    //Create a face
    face f1(1, 2, 3, 0);

    //Reset the id of the face
    f1.set_local_id(1);

    return !(f1.get_local_id() == 1);
}
//---------------------------------------------------------------------------------------------------------



//---------------------------------------------------------------------------------------------------------

//Check if the ids of the nodes of a face can correctly be updated
int face_tester::update_node_ids_test(){

    //Create arbitrary nodes
    node n1(1), n2(2), n3(3);

    //Create a face
    face f1(n1, n2, n3, 0);

    //Create a map with the old <--> new node ids correspondence
    std::map<unsigned, unsigned> node_id_correspondence;
    node_id_correspondence.insert({1, 10});
    node_id_correspondence.insert({2, 20});
    node_id_correspondence.insert({3, 30});

    //Update the ids of each node
    f1.update_node_ids(node_id_correspondence);

    //Chekc that t1 is connected to the correct nodes
    const auto& [f1_n1_id, f1_n2_id, f1_n3_id] = f1.get_node_ids();
    return !(f1_n1_id == 10 && f1_n2_id == 20 && f1_n3_id == 30);
}
//---------------------------------------------------------------------------------------------------------




//---------------------------------------------------------------------------------------------------------
//Check if the face is used by nodes or not
int face_tester::is_used_test(){

    //Create arbitrary nodes
    node n1(1), n2(2), n3(3);

    //Create a face
    face f1(n1, n2, n3, 0);
    
    bool t1 =  f1.is_used() ==  true;

    //Remove the nodes from the face
    f1.reset();

    bool t2 = f1.is_used() == false;
    
    return !(t1 && t2);

}
//---------------------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------------------
int face_tester::swap_nodes_test(){
    face f1(0, 1, 2, 0);
    f1.swap_nodes();
    const auto [n1, n2, n3] =  f1.get_node_ids();
    return !(n1 == 0, n2 == 2, n3 == 1);
}
//---------------------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------------------
int face_tester::get_opposite_node_test(){
    face f1(0, 1, 2, 0);
    bool t1 = f1.get_opposite_node(0, 1) == 2;
    bool t2 = f1.get_opposite_node(1, 2) == 0;
    bool t3 = f1.get_opposite_node(2, 0) == 1;

    return !(t1 && t2 && t3);
}
//---------------------------------------------------------------------------------------------------------



//---------------------------------------------------------------------------------------------------------
int face_tester::replace_node_test(){
    face f1(0, 1, 2, 0);

    f1.replace_node(0, 3);
    f1.replace_node(1, 4);

    const auto [n1, n2, n3] =  f1.get_node_ids();
    return !(n1 == 3, n2 == 4, n3 == 2);
}
//---------------------------------------------------------------------------------------------------------



//---------------------------------------------------------------------------------------------------------
int face_tester::face_class_size_test(){
    face f1(0, 1, 2, 0);

    std::cout << "Size of face class: " << sizeof(f1) << std::endl;

    return 0;
}
//---------------------------------------------------------------------------------------------------------





//---------------------------------------------------------------------------------------------------------
// The main function
int main (int argc, char** argv){

    //Check that the command line input is correctly formatted
    assert(argc == 2); 

    //Get the name of the test to run
    std::string test_name = argv[1];

    face_tester tester;
    
    //Run the selected test
    if (test_name == "default_constructors_test")                   return tester.default_constructors_test();
    if (test_name == "copy_constructor_test")                       return tester.copy_constructor_test();
    if (test_name == "node_id_constructor_test")                    return tester.node_id_constructor_test();
    if (test_name == "reset_local_id_test")                         return tester.reset_local_id_test();
    if (test_name == "update_node_ids_test")                        return tester.update_node_ids_test();
    if (test_name == "is_used_test")                                return tester.is_used_test();
    if (test_name == "swap_nodes_test")                             return tester.swap_nodes_test();
    if (test_name == "get_opposite_node_test")                      return tester.get_opposite_node_test();
    if (test_name == "replace_node_test")                           return tester.replace_node_test();
    if (test_name == "face_class_size_test")                        return tester.face_class_size_test();


    

    std::cout << "TEST NAME :" << test_name << " DOES NOT EXIST" << std::endl;
    return 1;
}
//---------------------------------------------------------------------------------------------------------
