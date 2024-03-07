cmake_minimum_required(VERSION 3.0)
set(CMAKE_CXX_STANDARD 17)


#Compile the files containing the tests to run
add_executable(test_bpa test_bpa.cpp)


target_link_libraries(test_bpa PUBLIC          
    src 
    mesh
    triangulation_modules
    uspg
)


#Indicate to main the path to project source directory
add_definitions(-DPROJECT_SOURCE_DIR="${PROJECT_SOURCE_DIR}")

target_include_directories(test_bpa PUBLIC  
    "${PROJECT_SOURCE_DIR}/include/"
    "${PROJECT_SOURCE_DIR}/include/math_modules"     
    "${PROJECT_SOURCE_DIR}/include/uspg"
    "${PROJECT_SOURCE_DIR}/include/mesh"
    "${PROJECT_SOURCE_DIR}/include/triangulation_modules"                                                                    
)

#Run the test
add_test(NAME triangulation_modules_bpa_points_coherently_oriented_test          COMMAND test_bpa points_coherently_oriented_test)
add_test(NAME triangulation_modules_bpa_compute_oriented_normal_test             COMMAND test_bpa compute_oriented_normal_test)
add_test(NAME triangulation_modules_bpa_node_is_compatible_with_edge_test        COMMAND test_bpa node_is_compatible_with_edge_test)
add_test(NAME triangulation_modules_bpa_compute_angle_between_circumspheres_test COMMAND test_bpa compute_angle_between_circumspheres_test)
add_test(NAME triangulation_modules_bpa_compute_circum_sphere_center_test        COMMAND test_bpa compute_circum_sphere_center_test)
add_test(NAME triangulation_modules_bpa_circum_sphere_is_empty_test              COMMAND test_bpa circum_sphere_is_empty_test)
add_test(NAME triangulation_modules_bpa_get_edge_test                            COMMAND test_bpa get_edge_test)
add_test(NAME triangulation_modules_bpa_node_is_inner_vertex_test                COMMAND test_bpa node_is_inner_vertex_test)
#add_test(NAME triangulation_modules_bpa_fill_surface_holes_test                  COMMAND test_bpa fill_surface_holes_test)





































