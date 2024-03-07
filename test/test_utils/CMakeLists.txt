cmake_minimum_required(VERSION 3.0)
set(CMAKE_CXX_STANDARD 17)


#Compile the files containing the tests to run
add_executable(test_utils test_utils.cpp)

target_link_libraries(test_utils PUBLIC src)

target_include_directories(test_utils PUBLIC  "${PROJECT_SOURCE_DIR}/include")

#Run the test
add_test(NAME test_utils_remove_index_test          COMMAND test_utils remove_index_test)
add_test(NAME test_utils_flatten_vector_test        COMMAND test_utils flatten_vector_test)









