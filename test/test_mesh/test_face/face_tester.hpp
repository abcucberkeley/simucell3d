#ifndef DEF_FACE_TESTER
#define DEF_FACE_TESTER

#include <cassert>
#include <string>
#include <type_traits>
#include <vector>
#include <algorithm>
#include <functional>
#include <map>


#include "cell.hpp"
#include "vec3.hpp"


class face_tester
{   
    public:

        int default_instantiation_test();
        int default_constructors_test();
        int copy_constructor_test();
        int node_id_constructor_test();
        int reset_local_id_test();
        int update_node_ids_test();
        int is_used_test();
        int swap_nodes_test();
        int get_opposite_node_test();
        int replace_node_test();
        int face_class_size_test();
 
};

#endif