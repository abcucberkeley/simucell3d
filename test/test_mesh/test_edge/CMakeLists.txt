cmake_minimum_required(VERSION 3.0)
set(CMAKE_CXX_STANDARD 17)


#Compile the files containing the tests to run
add_executable(test_edge
    test_edge.cpp
)

target_link_libraries(test_edge PUBLIC  src 
                                        mesh
)

target_include_directories(test_edge PUBLIC  "${PROJECT_SOURCE_DIR}/include/"
                                             "${PROJECT_SOURCE_DIR}/include/math_modules"
                                             "${PROJECT_SOURCE_DIR}/include/mesh"
                                        
)



#Run the test
add_test(NAME mesh_edge_default_constructors_test       COMMAND test_edge default_constructors_test)
add_test(NAME mesh_edge_trivial_constructor             COMMAND test_edge trivial_constructor_test)
add_test(NAME mesh_edge_hash_test                       COMMAND test_edge edge_hash_test)
add_test(NAME mesh_edge_unordered_set_test              COMMAND test_edge edge_unordered_set_test)
add_test(NAME mesh_edge_add_face_test                   COMMAND test_edge add_face_test)
add_test(NAME mesh_edge_is_manifold_test                COMMAND test_edge is_manifold_test)
add_test(NAME mesh_edge_swap_face_ids_test              COMMAND test_edge swap_face_ids_test)
















