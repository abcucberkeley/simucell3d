#ifndef DEF_CELL_TESTER
#define DEF_CELL_TESTER

#include <cassert>
#include <string>
#include <type_traits>
#include <vector>
#include <algorithm>

#include "global_configuration.hpp"

#include "cell.hpp"
#include "epithelial_cell.hpp"
#include "vec3.hpp"
#include "mesh_writer.hpp"

typedef std::shared_ptr<cell> cell_ptr;


class cell_tester
{   
    public:

        int deleted_constructors_test();
        int trivial_default_constructor_1_test();
        int trivial_default_constructor_2_test();
        int instantiation_from_mesh_test();
        int connect_cell_to_nodes_and_faces();
        int edge_set_generation_test();
        int rebase_node_test();
        int rebase_face_test();
        int get_flat_node_coord_lst_test();
        int get_face_normal_test();
        int check_face_normal_orientation_test();
        int compute_centroid_test();
        int compute_volume_test();      
        int is_manifold_test();
        int update_face_area_test();
        int get_aabb_test();
        int compute_area_test();
        int all_nodes_are_used_test();
        int get_edge_test();
        int delete_face_test();
        int add_face_test_1();
        int add_face_test_2();
        int add_node_test();
        int delete_node_test();
        int replace_node_test();
        int copy_constructor_test();
        int get_node_faces_test();
        int get_cell_longest_axis_test();
        int get_angle_gradient_test();

};

#endif