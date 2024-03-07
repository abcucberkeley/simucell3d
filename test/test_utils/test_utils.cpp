#include <cassert>
#include <string>
#include <iostream>
#include <vector>

#include "utils.hpp"

int remove_index_test(){

    //Create a vector with dummy double values
    std::vector<double> test_vec{0.0, 0.1, 0.2, 0.3, 0.4, 0.5};

    //The indices of the elements to dellte in the test vector
    std::vector<unsigned> index_lst{1, 2, 4};

    remove_index(test_vec, index_lst);

    std::vector<double> correct_result{0.0, 0.3, 0.5};

    bool t1 = std::equal(correct_result.begin(), correct_result.end(), test_vec.begin());
    
    return !(t1);
}


int flatten_vector_test(){
    std::vector<std::vector<int>> test_vec{
        {1,2,3},
        {4, 5},
        {6},
        {7,8,9}
    };


    std::vector<int> flattened_vec = flatten(test_vec);

    std::vector<double> correct_result{1, 2, 3, 4, 5, 6, 7, 8, 9};

    bool t1 = std::equal(correct_result.begin(), correct_result.end(), flattened_vec.begin());
    
    return !(t1);
}




//---------------------------------------------------------------------------------------------------------
// The main function
int main (int argc, char** argv){

    //Check that the command line input is correctly formatted
    assert(argc == 2); 

    //Get the name of the test to run
    std::string test_name = argv[1];
    
    if (test_name == "remove_index_test")     return remove_index_test();
    if (test_name == "flatten_vector_test")   return flatten_vector_test();

    std::cout << "TEST NAME :" << test_name << " DOES NOT EXIST" << std::endl;
    return 1;

}
//---------------------------------------------------------------------------------------------------------
