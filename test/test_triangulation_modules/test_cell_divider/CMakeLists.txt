cmake_minimum_required(VERSION 3.0)
set(CMAKE_CXX_STANDARD 17)


#Compile the files containing the tests to run
add_executable(test_cell_divider test_cell_divider.cpp)

#Indicate to main the path to project source directory
add_definitions(-DPROJECT_SOURCE_DIR="${PROJECT_SOURCE_DIR}")

target_link_libraries(test_cell_divider PUBLIC     src 
                                                            io
                                                            mesh
                                                            triangulation_modules
                                                            uspg
)

target_include_directories(test_cell_divider PUBLIC  
    "${PROJECT_SOURCE_DIR}/include/"
    "${PROJECT_SOURCE_DIR}/include/io"
    "${PROJECT_SOURCE_DIR}/include/math_modules"     
    "${PROJECT_SOURCE_DIR}/include/uspg"
    "${PROJECT_SOURCE_DIR}/include/mesh"
    "${PROJECT_SOURCE_DIR}/include/mesh/cell_types"
    "${PROJECT_SOURCE_DIR}/include/triangulation_modules"                                                                    
)

#Run the test
add_test(NAME triangulation_modules_cell_divider_add_intersection_point_test            COMMAND test_cell_divider add_intersection_point_test)
add_test(NAME triangulation_modules_cell_divider_divide_faces_test                      COMMAND test_cell_divider divide_faces_test)
add_test(NAME triangulation_modules_cell_divider_add_division_interface_test            COMMAND test_cell_divider add_division_interface_test)
add_test(NAME triangulation_modules_cell_divider_map_points_to_xy_plane_test            COMMAND test_cell_divider map_points_to_xy_plane_test)
add_test(NAME triangulation_modules_cell_divider_triangulate_division_interface_test    COMMAND test_cell_divider triangulate_division_interface_test)
add_test(NAME triangulation_modules_cell_divider_create_daughter_cell_meshes_test       COMMAND test_cell_divider create_daughter_cell_meshes_test)




































