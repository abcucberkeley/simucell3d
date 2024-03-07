#ifndef DEF_TEST_PARAMETER_READER
#define DEF_TEST_PARAMETER_READER


#include <cassert>
#include <string>
#include <iostream>


#include "parameter_reader.hpp"



class test_parameter_reader
{   

    public:
        test_parameter_reader() = default;

        int read_cell_type_parameters_test() const;
        int read_face_type_parameters_test() const;
        int read_biomechanical_parameters_test() const;
    


};

#endif