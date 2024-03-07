

cmake_minimum_required(VERSION 3.0)
set(CMAKE_CXX_STANDARD 17)


#Compile the files containing the tests to run
add_executable(test_uspg_4d test_uspg_4d.cpp)

target_link_libraries(test_uspg_4d PUBLIC uspg)


target_include_directories(test_uspg_4d PUBLIC  
                                        "${PROJECT_SOURCE_DIR}/test/test_uspg/test_uspg_4d/"  
                                        "${PROJECT_SOURCE_DIR}/include/**"
                                        "${PROJECT_SOURCE_DIR}/include/uspg"
                                        "${PROJECT_SOURCE_DIR}/include/math_modules"
                                      
)



#Run the test
add_test(NAME uspg_4d_update_dimensions_test       COMMAND test_uspg_4d update_dimensions_test)
add_test(NAME uspg_4d_get_neighborhood_test             COMMAND test_uspg_4d get_neighborhood_test)
add_test(NAME uspg_4d_place_object_test                 COMMAND test_uspg_4d place_object_test)
add_test(NAME uspg_4d_get_grid_content_test             COMMAND test_uspg_4d get_grid_content_test
)





