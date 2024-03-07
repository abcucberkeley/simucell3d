#include <cassert>
#include <string>
#include <type_traits>
#include <unordered_set>
#include <algorithm>


#include "edge.hpp"

/*

Test the method of the edge class

*/


typedef std::shared_ptr<cell> cell_ptr;

//---------------------------------------------------------------------------------------------------------
//Test the default constructors of the edge class
int default_constructors_test(){

    //Check that the default instantiation is not possible
    bool t1 =  std::is_default_constructible<edge>::value;

    //Test the copy constructor
    edge e1(1, 2, 3, 4);
    edge e2(e1);

    const auto [n_id1, n_id2] = e2.get_node_ids();
    bool t2 = n_id1 == 1 && n_id2 == 2; 

    const auto [f_id1, f_id2] = e2.get_face_ids();
    bool t3 = f_id1 == 3 && f_id2 == 4; 

    //Test the copy assignment operator
    edge e3 = e1;
    const auto [n_id3, n_id4] = e3.get_node_ids();
    bool t4 = n_id3 == 1 && n_id4 == 2; 

    const auto [f_id3, f_id4] = e3.get_face_ids();
    bool t5 = f_id3 == 3 && f_id4 == 4; 

    return !(t1 && t2 && t3 && t4 && t5);
}
//---------------------------------------------------------------------------------------------------------



//---------------------------------------------------------------------------------------------------------
//Test the constructor of the face class
int trivial_constructor_test(){
   
    //Create arbitrary nodes
    node n1(1), n2(2), n3(3);

    //Create arbitrary faces
    face f1(n1, n2, n3, 0), f2(n1, n2, n3, 1);

    //Create an edge with the copy constructor
    edge e1(n2, n1);

    //Check the ids of the nodes
    const auto [id1, id2] = e1.get_node_ids();
    bool t1 = (id1 == 1 && id2 == 2);

    //Create an edge with the node ids
    edge e2(4, 3);
    const auto [id3, id4] = e2.get_node_ids();
    bool t2 = (id3 == 3 && id4 == 4);

    //Create an edge with the node and face ids
    edge e3(4, 3, 6, 7);
    const auto [id5, id6] = e3.get_node_ids();
    bool t3 = (id5 == 3 && id6 == 4);

    const auto [id7, id8] = e3.get_face_ids();
    bool t4 = (id7 == 6 && id8 == 7);

    return !(t1 && t2 && t3 && t4);

}
//---------------------------------------------------------------------------------------------------------



//---------------------------------------------------------------------------------------------------------




//---------------------------------------------------------------------------------------------------------
//Check the edge hasher
int edge_hash_test(){

    //Create arbitraty edges
    edge e1(1, 2), e2(2, 1), e3(1, 4);

    edge_hasher hasher;
    int hash_1 = hasher(e1);
    int hash_2 = hasher(e2);
    int hash_3 = hasher(e3);

    bool t1 = hash_1 == hash_2;
    bool t2 = hash_1 != hash_3;
    
    return !(t1 && t2);
}

//---------------------------------------------------------------------------------------------------------


//---------------------------------------------------------------------------------------------------------
//Check the edge comparator
int edge_comparator_test(){


    //Create arbitraty edges
    edge e1(1, 2), e2(2, 1), e3(1, 4);
    edge_comparator comparator;

    bool eq1 = comparator(e1, e1);
    bool eq2 = comparator(e1, e2);
    bool eq3 = comparator(e1, e3);

    return !(eq1 && eq2 && !eq3);




}
//---------------------------------------------------------------------------------------------------------






//---------------------------------------------------------------------------------------------------------
//Try to create an unordered_set of edges
int edge_unordered_set_test(){


    //Create arbitraty edges
    edge e1(1, 2), e2(2, 1), e3(1, 4);

    //Create an edge set
    std::unordered_set<edge, edge_hasher,  edge_comparator>  edge_set;
    
    //Add edges to the set
    edge_set.insert(e1);
    edge_set.insert(e2);
    edge_set.insert(e3);

    //Test the sizde of the set.
    bool t1 = edge_set.size() == 2;

    //Check that all the edges declared are in it
    bool t2 = std::find(edge_set.begin(), edge_set.end(), e1 ) !=  edge_set.end();
    bool t3 = std::find(edge_set.begin(), edge_set.end(), e2 ) !=  edge_set.end();
    bool t4 = std::find(edge_set.begin(), edge_set.end(), e3 ) !=  edge_set.end();


    //Check that only 2 elements are in the set
    return !(t1 && t2 && t3 && t4);
}
//---------------------------------------------------------------------------------------------------------



//---------------------------------------------------------------------------------------------------------
//Test the is manifold function
int is_manifold_test(){

    edge e1(1,2,3,4);
    edge e2(1, 2);
    edge e3(4, 5);
    e3.add_face(6);


    bool t1 = e1.is_manifold();
    bool t2 = !(e2.is_manifold());
    bool t3 = !(e3.is_manifold());

    return !(t1 && t2 && t3);

}
//---------------------------------------------------------------------------------------------------------
int swap_face_ids_test(){
    edge e1(1,2,3,4);

    const auto f1_id_before = e1.f1();
    const auto f2_id_before = e1.f2();

    bool t1 = f1_id_before == 3 && f2_id_before == 4;

    e1.swap_face_ids();

    const auto f1_id_after = e1.f1();
    const auto f2_id_after = e1.f2();

    bool t2 = f1_id_after == 4 && f2_id_after == 3;

    return !(t1 && t2);


}


//---------------------------------------------------------------------------------------------------------
int add_face_test(){

    edge e1(1,2);

    e1.add_face(3);
    bool t1 = e1.f1() == 3;

    e1.add_face(4);
    bool t2 = e1.f2() == 4;

    bool t3 = false;
    try {e1.add_face(6);}
    catch(const mesh_integrity_exception& e){t3 = true;}

    return !(t1 && t2 && t3);

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
    if (test_name == "default_constructors_test")                   return default_constructors_test();
    if (test_name == "trivial_constructor_test")                    return trivial_constructor_test();
    if (test_name == "edge_hash_test")                              return edge_hash_test();
    if (test_name == "edge_comparator_test")                        return edge_comparator_test();
    if (test_name == "edge_unordered_set_test")                     return edge_unordered_set_test();
    if (test_name == "add_face_test")                               return add_face_test();
    if (test_name == "is_manifold_test")                            return is_manifold_test();
    if (test_name == "swap_face_ids_test")                          return swap_face_ids_test();

    std::cout << "TEST NAME :" << test_name << " DOES NOT EXIST" << std::endl;
    return 1;
}
//---------------------------------------------------------------------------------------------------------
