cmake_minimum_required(VERSION 3.0)
set(CMAKE_CXX_STANDARD 17)


#Compile the files containing the tests to run
add_executable(test_local_mesh_refiner
    test_local_mesh_refiner.cpp
)

target_link_libraries(test_local_mesh_refiner PUBLIC  
    src 
    mesh
    triangulation_modules
    io
)

target_include_directories(test_local_mesh_refiner PUBLIC     
    "${PROJECT_SOURCE_DIR}/include/"
    "${PROJECT_SOURCE_DIR}/include/triangulation_modules"
    "${PROJECT_SOURCE_DIR}/include/math_modules"
    "${PROJECT_SOURCE_DIR}/include/mesh"
    "${PROJECT_SOURCE_DIR}/test/test_triangulation_modules/test_local_mesh_refiner"                                                                           
    "${PROJECT_SOURCE_DIR}/include/io"
)



#Run the test
add_test(NAME triangulation_modules_local_mesh_refiner_can_be_merged      COMMAND test_local_mesh_refiner can_be_merged_test)
#add_test(NAME triangulation_modules_local_mesh_refiner_prepare_merger_test_1    COMMAND test_local_mesh_refiner prepare_merger_test_1)
#add_test(NAME triangulation_modules_local_mesh_refiner_prepare_merger_test_2    COMMAND test_local_mesh_refiner prepare_merger_test_2)
add_test(NAME triangulation_modules_local_mesh_refiner_edge_swap_test           COMMAND test_local_mesh_refiner edge_swap_test)
add_test(NAME triangulation_modules_local_mesh_refiner_split_edge_test          COMMAND test_local_mesh_refiner split_edge_test)
add_test(NAME triangulation_modules_local_mesh_refiner_merge_edge_test          COMMAND test_local_mesh_refiner merge_edge_test)
add_test(NAME triangulation_modules_local_mesh_refiner_get_triangle_score_test  COMMAND test_local_mesh_refiner get_triangle_score_test)


