cmake_minimum_required(VERSION 3.0)

#Use CPP 17
set(CMAKE_CXX_STANDARD 17)

#Indicate to main the path to project source directory
add_definitions(-DPROJECT_SOURCE_DIR="${PROJECT_SOURCE_DIR}")


#Compile the files containing the tests to run
add_executable(test_simulation_initializer test_simulation_initializer.cpp)

target_link_libraries(test_simulation_initializer PUBLIC   src
                                                mesh
                                                io
)


#Add all the header files required to run the root file of the project
target_include_directories(test_simulation_initializer PUBLIC 
    "${PROJECT_SOURCE_DIR}/include/"
    "${PROJECT_SOURCE_DIR}/include/mesh"
    "${PROJECT_SOURCE_DIR}/include/io"
)

#Run the test
add_test(NAME io_simulation_initializer_simulation_initializer_constructor_test                 COMMAND test_simulation_initializer simulation_initializer_constructor_test)










