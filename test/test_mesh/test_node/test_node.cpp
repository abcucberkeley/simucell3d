#include "node_tester.hpp"


/*

Test the method of the node class

*/


typedef std::shared_ptr<cell> cell_ptr;

//---------------------------------------------------------------------------------------------------------
//Test the default constructors of the node class
int node_tester::default_constructors_test(){


    //Check the default copy constructor

    //Create a dummy face
    node n1(1., 2., 3., 0);

    //This copy has the same cell_ptr but is not registered in the cell as one of the faces
    node n2(n1);
    auto [coord_1, coord_2, coord_3 ] = n2.pos().to_array();
    bool t2 = coord_1 == 1. &&  coord_2 == 2. && coord_3 == 3.;
    bool t3 = n1.get_local_id() == n2.get_local_id();

    //Do the same test with the copy assignment function
    node n3 =  n1;
    auto [coord_4, coord_5, coord_6 ] = n3.pos().to_array();
    bool t4 = coord_4 == 1. &&  coord_5 == 2. && coord_6 == 3.;
    bool t5 = n1.get_local_id() == n3.get_local_id();
    return !(t2 && t3 && t4 && t5);
}

//---------------------------------------------------------------------------------------------------------
//Test the constructor of the node class
int node_tester::trivial_instantiation_test(){

    //Generate a node with default coordinates 
    node n1(1);

    //Check the id of the node 
    bool t1 = n1.get_local_id() == 1;

    //Check the position of the node
    bool t2 = (n1.pos().dx() == 0. && n1.pos().dy() == 0. && n1.pos().dz() == 0.);

    //Create a node with arbitrary coordinates
    node n2(1., 2. , 3., 2);

    //Check the id of the node 
    bool t3 = n2.get_local_id() == 2;

    //Check the position of the node
    bool t4 = (n2.pos().dx() == 1. && n2.pos().dy() == 2. && n2.pos().dz() == 3.);

    return !(t1 && t2 && t3 && t4);
}
//---------------------------------------------------------------------------------------------------------



//---------------------------------------------------------------------------------------------------------
//Test the reset of a node
int node_tester::node_reset_test(){

    //Generate a node with default coordinates 
    node n0(0), n1(1), n2(2), n3(3);
    std::vector<node> node_lst = {n0, n1, n2, n3};
    
    face f1(n0, n1, n2, 0);
    std::vector<face> face_lst{f1};

    //Create a dummy cell
    cell_ptr c1 = std::make_shared<cell>(node_lst, face_lst, 1);

    //Reset the first and third node of the cell
    c1->node_lst_[0].reset();
    c1->node_lst_[3].reset();

    //Check that these two nodes are not used anymore by the cel
    bool t1 = c1->node_lst_[0].is_used_ == false;
    bool t2 = c1->node_lst_[3].is_used_ == false;

    return !(t1 && t2);
}
//---------------------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------------------
int node_tester::reset_local_id_test(){

    //Generate a node with default coordinates 
    node n1(1);

    //Reset its id
    n1.set_local_id(0);

    return !(n1.get_local_id() == 0);
}

//---------------------------------------------------------------------------------------------------------


//---------------------------------------------------------------------------------------------------------
int node_tester::set_is_used_test(){

    //Generate a node with default coordinates 
    node n1(1);

    bool t1 = n1.is_used() == true;

    n1.set_is_used(false);

    bool t2 = n1.is_used() == false;

    return !(t1 && t2);
}

//---------------------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------------------
int node_tester::operator_minus_test(){

    node n1(0., 0., 0., 1);
    node n2(1., 1., 1., 2);
    vec3 v = n2 - n1;

    return !(v == vec3(1, 1, 1));
}

//---------------------------------------------------------------------------------------------------------




//---------------------------------------------------------------------------------------------------------
// The main function
int main (int argc, char** argv){

    //Check that the command line input is correctly formatted
    assert(argc == 2); 

    //Get the name of the test to run
    std::string test_name = argv[1];

    node_tester tester;
    
    //Run the selected test
    if (test_name == "default_constructors_test")         return tester.default_constructors_test();
    if (test_name == "trivial_instantiation_test")        return tester.trivial_instantiation_test();
    if (test_name == "node_reset_test")                   return tester.node_reset_test();
    if (test_name == "reset_local_id_test")                     return tester.reset_local_id_test();
    if (test_name == "set_is_used_test")                  return tester.set_is_used_test();
    if (test_name == "operator_minus_test")               return tester.operator_minus_test();




    std::cout << "TEST NAME :" << test_name << " DOES NOT EXIST" << std::endl;
    return 1;
}
//---------------------------------------------------------------------------------------------------------

