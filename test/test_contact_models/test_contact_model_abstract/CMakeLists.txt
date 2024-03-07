cmake_minimum_required(VERSION 3.0)

#Use CPP 17
set(CMAKE_CXX_STANDARD 17)

#Indicate to main the path to project source directory
add_definitions(-DPROJECT_SOURCE_DIR="${PROJECT_SOURCE_DIR}")


#Compile the files containing the tests to run
add_executable(test_contact_model_abstract test_contact_model_abstract.cpp)

target_link_libraries(test_contact_model_abstract PUBLIC   
    src
    mesh
    contact_models
)

#Add all the header files required to run the root file of the project
target_include_directories(test_contact_model_abstract PUBLIC 
    "${PROJECT_SOURCE_DIR}/include/"
    "${PROJECT_SOURCE_DIR}/include/mesh"
    "${PROJECT_SOURCE_DIR}/include/contact_models"
)

#Run the test
add_test(NAME contact_models_abstract_compute_node_triangle_distance_test  COMMAND test_contact_model_abstract compute_node_triangle_distance_test)
