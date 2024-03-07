cmake_minimum_required(VERSION 3.0)
set(CMAKE_CXX_STANDARD 17)


#Compile the files containing the tests to run
add_executable(test_cell
    test_cell.cpp
)

target_link_libraries(test_cell PUBLIC  src 
                                        mesh
)

target_include_directories(test_cell PUBLIC  
                                        "${PROJECT_SOURCE_DIR}/include/"
                                        "${PROJECT_SOURCE_DIR}/include/math_modules"
                                        "${PROJECT_SOURCE_DIR}/include/mesh"
                                        "${PROJECT_SOURCE_DIR}/test/test_mesh/test_cell"                                        

                                        
)



#Run the test
add_test(NAME mesh_cell_deleted_constructors_test           COMMAND test_cell deleted_constructors_test)
add_test(NAME mesh_cell_trivial_default_constructor_1_test  COMMAND test_cell trivial_default_constructor_1_test)
add_test(NAME mesh_cell_trivial_default_constructor_2_test  COMMAND test_cell trivial_default_constructor_2_test)
add_test(NAME mesh_cell_edge_set_generation_test            COMMAND test_cell edge_set_generation_test)
add_test(NAME mesh_cell_rebase_node_test                    COMMAND test_cell rebase_node_test)
add_test(NAME mesh_cell_rebase_face_test                    COMMAND test_cell rebase_face_test)
add_test(NAME mesh_cell_get_flat_node_coord_lst_test             COMMAND test_cell get_flat_node_coord_lst_test)
add_test(NAME mesh_cell_get_face_normal_test                COMMAND test_cell get_face_normal_test)
add_test(NAME mesh_cell_check_face_normal_orientation_test  COMMAND test_cell check_face_normal_orientation_test)
add_test(NAME mesh_cell_compute_centroid_test               COMMAND test_cell compute_centroid_test)
add_test(NAME mesh_cell_compute_volume_test                 COMMAND test_cell compute_volume_test)
add_test(NAME mesh_cell_is_manifold_test                    COMMAND test_cell is_manifold_test)
add_test(NAME mesh_cell_update_face_area_test               COMMAND test_cell update_face_area_test)
add_test(NAME mesh_cell_get_aabb_test                       COMMAND test_cell get_aabb_test)
add_test(NAME mesh_cell_compute_area_test                   COMMAND test_cell compute_area_test)
add_test(NAME mesh_cell_instantiation_from_mesh_test        COMMAND test_cell instantiation_from_mesh_test)
add_test(NAME mesh_cell_all_nodes_are_used_test             COMMAND test_cell all_nodes_are_used_test)
add_test(NAME mesh_cell_get_edge_test                       COMMAND test_cell get_edge_test)
add_test(NAME mesh_cell_delete_face_test                    COMMAND test_cell delete_face_test)
add_test(NAME mesh_cell_add_face_test_1                     COMMAND test_cell add_face_test_1)
add_test(NAME mesh_cell_add_face_test_2                     COMMAND test_cell add_face_test_2)
add_test(NAME mesh_cell_add_node_test                       COMMAND test_cell add_node_test)
add_test(NAME mesh_cell_delete_node_test                    COMMAND test_cell delete_node_test)
add_test(NAME mesh_cell_replace_node_test                   COMMAND test_cell replace_node_test)
add_test(NAME mesh_cell_copy_constructor_test               COMMAND test_cell copy_constructor_test)
add_test(NAME mesh_cell_get_node_faces_test                 COMMAND test_cell get_node_faces_test)
add_test(NAME mesh_cell_get_cell_longest_axis_test          COMMAND test_cell get_cell_longest_axis_test)
add_test(NAME mesh_cell_get_angle_gradient_test             COMMAND test_cell get_angle_gradient_test)




























