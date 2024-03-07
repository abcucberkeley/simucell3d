cmake_minimum_required(VERSION 3.0)
set(CMAKE_CXX_STANDARD 17)


#Compile the files containing the tests to run
add_executable(test_face
    test_face.cpp
)

target_link_libraries(test_face PUBLIC  src 
                                        mesh
)

target_include_directories(test_face PUBLIC     "${PROJECT_SOURCE_DIR}/include/"
                                                "${PROJECT_SOURCE_DIR}/include/math_modules"
                                                "${PROJECT_SOURCE_DIR}/include/mesh"
                                                "${PROJECT_SOURCE_DIR}/test/test_mesh/test_face"                                        
)

#Run the test
add_test(NAME mesh_face_default_constructors_test       COMMAND test_face default_constructors_test)
add_test(NAME mesh_face_node_id_constructor             COMMAND test_face node_id_constructor_test)
add_test(NAME mesh_face_copy_constructor                COMMAND test_face copy_constructor_test)
add_test(NAME mesh_face_reset_local_id_test                   COMMAND test_face reset_local_id_test)
add_test(NAME mesh_face_update_node_ids_test            COMMAND test_face update_node_ids_test)
add_test(NAME mesh_face_is_used_test                    COMMAND test_face is_used_test)
add_test(NAME mesh_face_swap_nodes_test                 COMMAND test_face swap_nodes_test)
add_test(NAME mesh_face_get_opposite_node_test          COMMAND test_face get_opposite_node_test)
add_test(NAME mesh_face_replace_node_test               COMMAND test_face replace_node_test)








