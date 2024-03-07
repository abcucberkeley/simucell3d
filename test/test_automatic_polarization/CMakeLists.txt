cmake_minimum_required(VERSION 3.0)

#Use CPP 17
set(CMAKE_CXX_STANDARD 17)

#Indicate to main the path to project source directory
add_definitions(-DPROJECT_SOURCE_DIR="${PROJECT_SOURCE_DIR}")


#Compile the files containing the tests to run
add_executable(test_automatic_polarization test_automatic_polarization.cpp)

target_link_libraries(test_automatic_polarization PUBLIC   
    src
    mesh
)

#Add all the header files required to run the root file of the project
target_include_directories(test_automatic_polarization PUBLIC 
    "${PROJECT_SOURCE_DIR}/include/"
    "${PROJECT_SOURCE_DIR}/include/mesh"
    "${PROJECT_SOURCE_DIR}/include/automatic_polarization"
    "${PROJECT_SOURCE_DIR}/include/io"
)






