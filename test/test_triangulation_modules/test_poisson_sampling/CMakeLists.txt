cmake_minimum_required(VERSION 3.0)
set(CMAKE_CXX_STANDARD 17)


#Compile the files containing the tests to run
add_executable(test_poisson_sampling test_poisson_sampling.cpp)

target_link_libraries(test_poisson_sampling PUBLIC          src 
                                                            mesh
                                                            triangulation_modules
                                                            uspg
)

target_include_directories(test_poisson_sampling PUBLIC  
    "${PROJECT_SOURCE_DIR}/include/"
    "${PROJECT_SOURCE_DIR}/include/math_modules"     
    "${PROJECT_SOURCE_DIR}/include/uspg"
    "${PROJECT_SOURCE_DIR}/include/mesh"
    "${PROJECT_SOURCE_DIR}/include/triangulation_modules"                                                                    
)

#Run the test
add_test(NAME triangulation_modules_poisson_sampling_point_min_distance_test     COMMAND test_poisson_sampling point_min_distance_test)
add_test(NAME triangulation_modules_poisson_sampling_point_on_surface_test       COMMAND test_poisson_sampling point_on_surface_test)































