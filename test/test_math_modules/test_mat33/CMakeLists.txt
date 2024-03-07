cmake_minimum_required(VERSION 3.0)

#Use CPP 17
set(CMAKE_CXX_STANDARD 17)

#Compile the files containing the tests to run
add_executable(mat33_test
    mat33_test.cpp
)

target_link_libraries(mat33_test PUBLIC src)

#Add all the header files required to run the root file of the project
target_include_directories(mat33_test PUBLIC "${PROJECT_SOURCE_DIR}/include/")
target_include_directories(mat33_test PUBLIC "${PROJECT_SOURCE_DIR}/include/math_modules")

#Run the test
add_test(NAME math_modules_mat33_test_constructors          COMMAND mat33_test test_constructors)
add_test(NAME math_modules_mat33_get_row_col_test           COMMAND mat33_test get_row_col_test)
add_test(NAME math_modules_mat33_matrix_dot_test            COMMAND mat33_test matrix_dot_test)
add_test(NAME math_modules_mat33_vector_dot_test            COMMAND mat33_test vector_dot_test)
add_test(NAME math_modules_mat33_transpose_test             COMMAND mat33_test transpose_test)
add_test(NAME math_modules_mat33_determinant_test           COMMAND mat33_test determinant_test)
add_test(NAME math_modules_mat33_inverse_test               COMMAND mat33_test inverse_test)
add_test(NAME math_modules_mat33_operator_minus_test        COMMAND mat33_test operator_minus_test)
add_test(NAME math_modules_mat33_operator_plus_test         COMMAND mat33_test operator_plus_test)
add_test(NAME math_modules_mat33_operator_multiply_test     COMMAND mat33_test operator_multiply_test)
add_test(NAME math_modules_mat33_identity_test              COMMAND mat33_test identity_test)
























