cmake_minimum_required(VERSION 3.0)

#Use CPP 17
set(CMAKE_CXX_STANDARD 17)

#Indicate to main the path to the root directory and the build folder
add_definitions(-DCMAKE_BINARY_DIR="${CMAKE_BINARY_DIR}" -DPROJECT_SOURCE_DIR="${PROJECT_SOURCE_DIR}")

#Compile the files containing the tests to run
add_executable(test_nucleus_in_cell_contact_via_springs test_nucleus_in_cell_contact_via_springs.cpp)

#Add all the header files required to run the root file of the project
target_include_directories(test_nucleus_in_cell_contact_via_springs PUBLIC 
    "${PROJECT_SOURCE_DIR}/include/"
    "${PROJECT_SOURCE_DIR}/include/**"
)

#Link root library of the src directory with all the other libraries created
target_link_libraries(test_nucleus_in_cell_contact_via_springs PUBLIC    
    src
    uspg
    math_modules
    mesh
    io     
    triangulation_modules  
    time_integration  
    contact_models     
    automatic_polarization          
)
