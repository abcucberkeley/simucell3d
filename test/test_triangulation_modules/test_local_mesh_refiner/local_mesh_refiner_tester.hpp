#ifndef DEF_LOCAL_MESH_REFINER_TESTER
#define DEF_LOCAL_MESH_REFINER_TESTER

#include <cassert>
#include <string>
#include <type_traits>
#include <vector>
#include <algorithm>
#include <functional>

#include "local_mesh_refiner.hpp"
#include "mesh_writer.hpp"

#include "mesh_reader.hpp"



class local_mesh_refiner_tester
{   
    public:

        int can_be_merged_test();
        //int prepare_merger_test_1();
        //int prepare_merger_test_2();
        int edge_swap_test();
        int split_edge_test();
        int merge_edge_test();
        int merge_edge_test_2();
        int get_triangle_score_test();

};

#endif