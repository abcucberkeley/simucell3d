#ifndef DEF_TESTER_CONTACT_MODEL_ABSTRACT
#define DEF_TESTER_CONTACT_MODEL_ABSTRACT

#include <cassert>
#include <string>
#include <vector>
#include <algorithm>



#include "custom_structures.hpp"
#include "contact_model_abstract.hpp"

#include "vec3.hpp"



class tester_contact_model_abstract{   
    public:
        int compute_node_triangle_distance_test();
         
};

#endif