cmake_minimum_required(VERSION 3.0)
set(CMAKE_CXX_STANDARD 17)


#Compile the files containing the tests to run
add_executable(test_initial_triangulation test_initial_triangulation.cpp)

target_link_libraries(test_initial_triangulation PUBLIC     src 
                                                            io
                                                            mesh
                                                            triangulation_modules
                                                            uspg
)

target_include_directories(test_initial_triangulation PUBLIC  
    "${PROJECT_SOURCE_DIR}/include/"
    "${PROJECT_SOURCE_DIR}/include/io"
    "${PROJECT_SOURCE_DIR}/include/math_modules"     
    "${PROJECT_SOURCE_DIR}/include/uspg"
    "${PROJECT_SOURCE_DIR}/include/mesh"


    "${PROJECT_SOURCE_DIR}/include/triangulation_modules"                                                                    
)

#Run the test
add_test(NAME triangulation_modules_initial_triangulation_coarse_triangulation_test COMMAND test_initial_triangulation coarse_triangulation_test)



























