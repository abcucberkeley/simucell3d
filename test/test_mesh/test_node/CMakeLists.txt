cmake_minimum_required(VERSION 3.0)
set(CMAKE_CXX_STANDARD 17)


#Compile the files containing the tests to run
add_executable(test_node
    test_node.cpp
)

target_link_libraries(test_node PUBLIC  src 
                                        mesh
)

target_include_directories(test_node PUBLIC     "${PROJECT_SOURCE_DIR}/include/"
                                                "${PROJECT_SOURCE_DIR}/include/math_modules"
                                                "${PROJECT_SOURCE_DIR}/include/mesh"
                                                "${PROJECT_SOURCE_DIR}/test/test_mesh/test_node"                                        

                                        
)



#Run the test
add_test(NAME mesh_node_default_constructors_test           COMMAND test_node default_constructors_test)
add_test(NAME mesh_node_trivial_instantiation_reference     COMMAND test_node trivial_instantiation_test)
add_test(NAME mesh_node_reset_test                          COMMAND test_node node_reset_test)
add_test(NAME mesh_node_reset_local_id_test                       COMMAND test_node reset_local_id_test)
add_test(NAME mesh_node_set_is_used_test                    COMMAND test_node set_is_used_test)
add_test(NAME mesh_node_operator_minus_test                 COMMAND test_node operator_minus_test)






