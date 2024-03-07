cmake_minimum_required(VERSION 3.0)

#Use CPP 17
set(CMAKE_CXX_STANDARD 17)

#Compile the files containing the tests to run
add_executable(vec3_test
    vec3_test.cpp
)

target_link_libraries(vec3_test PUBLIC src)

#Add all the header files required to run the root file of the project
target_include_directories(vec3_test PUBLIC "${PROJECT_SOURCE_DIR}/include/")
target_include_directories(vec3_test PUBLIC "${PROJECT_SOURCE_DIR}/include/math_modules")

#Run the test
add_test(NAME math_modules_vec3_test_constructors           COMMAND vec3_test test_constructors)
add_test(NAME math_modules_vec3_coordinate_getter           COMMAND vec3_test coordinate_getter_test)
add_test(NAME math_modules_vec3_sum_lvalue                  COMMAND vec3_test sum_lvalue_test)
add_test(NAME math_modules_vec3_sum_rvalue                  COMMAND vec3_test sum_rvalue_test)
add_test(NAME math_modules_vec3_substraction_lvalue         COMMAND vec3_test substraction_lvalue_test)
add_test(NAME math_modules_vec3_substraction_rvalue         COMMAND vec3_test substraction_rvalue_test)
add_test(NAME math_modules_vec3_multiplication              COMMAND vec3_test multiplication_test)
add_test(NAME math_modules_vec3_division                    COMMAND vec3_test division_test)
add_test(NAME math_modules_vec3_norm                        COMMAND vec3_test norm_test)
add_test(NAME math_modules_vec3_normalize                   COMMAND vec3_test normalize_test)
add_test(NAME math_modules_vec3_dot_product                 COMMAND vec3_test dot_product_test)
add_test(NAME math_modules_vec3_equality                    COMMAND vec3_test equality_test)
add_test(NAME math_modules_vec3_cross_product               COMMAND vec3_test cross_product_test)
add_test(NAME math_modules_vec3_reset                       COMMAND vec3_test reset_test)
add_test(NAME math_modules_vec3_translate                   COMMAND vec3_test translate_test)
add_test(NAME math_modules_vec3_get_angle_with_test         COMMAND vec3_test get_angle_with_test)
add_test(NAME math_modules_vec3_rotate_around_axis_test     COMMAND vec3_test rotate_around_axis_test)
add_test(NAME math_modules_vec3_cap_test                    COMMAND vec3_test cap_test)









