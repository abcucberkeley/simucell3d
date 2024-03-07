#include <cassert>
#include <string>
#include <iostream>

#include "simulation_initializer.hpp"

//---------------------------------------------------------------------------------------------------------
int simulation_initializer_constructor_test(){
    
    simulation_initializer sim_init(std::string(PROJECT_SOURCE_DIR) + "/test/test_io/test_simulation_initializer/test_parameter_file.xml");

    //Get the list of cells
    std::vector<cell_ptr> cell_lst = sim_init.get_cell_lst();

    bool t1 = cell_lst.size() == 2;

    if(!t1){
        std::cout << "Number of cells: " << cell_lst.size() << std::endl;
        return 1;
    }

    std::cout << "Test 2" << std::endl;

    //Make sure the cells are of the correct type ptr 
    bool t2 = cell_lst[0]->get_cell_type() != nullptr;
    bool t3 = cell_lst[1]->get_cell_type() != nullptr;


    std::cout << "Test 1: " << t1 << std::endl;
    std::cout << "Test 2: " << t2 << std::endl;
    std::cout << "Test 3: " << t3 << std::endl;

    if(!t2 || !t3){
        std::cout << "cell type == nullptr" << std::endl;
        return 1;
    }


    bool t4 = cell_lst[0]->get_cell_type()->global_type_id_ == 0;
    bool t5 = cell_lst[1]->get_cell_type()->global_type_id_ == 1;



    //Make sure they have been cast to the correct type
    bool t6 = dynamic_cast<epithelial_cell*>(cell_lst[0].get()) != nullptr;
    bool t7 = dynamic_cast<ecm_cell*>(cell_lst[1].get()) != nullptr;


    //Check the id of the cell
    bool t8 = cell_lst[0]->get_id() == 0;
    bool t9 = cell_lst[1]->get_id() == 1;


    //Test that the faces know to which cell they belong
    bool t10  = true;
    for(cell_ptr c: cell_lst){
        if(!std::all_of(c->get_face_lst().begin(), c->get_face_lst().end(), [c](const face& f){return f.get_owner_cell() == c;})){
            t10 = false;
            break;
        }
    }


    std::cout << "Test 4: " << t4 << std::endl;
    std::cout << "Test 5: " << t5 << std::endl;
    std::cout << "Test 6: " << t6 << std::endl;
    std::cout << "Test 7: " << t7 << std::endl;
    std::cout << "Test 8: " << t8 << std::endl;
    std::cout << "Test 9: " << t9 << std::endl;
    std::cout << "Test 10: " << t10 << std::endl;



    return !(t1 && t2 && t3 && t4 && t5 && t6 && t7 && t8 && t9 && t10);
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
    if (test_name == "simulation_initializer_constructor_test")             return simulation_initializer_constructor_test();



    std::cout << "TEST NAME :" << test_name << " DOES NOT EXIST" << std::endl;
    return 1;
}
//---------------------------------------------------------------------------------------------------------



