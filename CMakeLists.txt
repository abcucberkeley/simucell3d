cmake_minimum_required(VERSION 3.0)

#Use CPP 17
set(CMAKE_CXX_STANDARD 17)


#The name of the project
project(DCM)

#Check if the python bindings should be compiled
if(NOT DEFINED ENABLE_PYTHON_BINDINGS)
    set(ENABLE_PYTHON_BINDINGS FALSE)
endif()

if(ENABLE_PYTHON_BINDINGS)
    message(STATUS "Python bindings enabled")
else()
    message(STATUS "Python bindings disabled")
endif()


#Activate unit testing
enable_testing()

#Contains the external library such as tinyxml (the xml reader)
add_subdirectory(lib)

#Contains the source files (.cpp)
add_subdirectory(src bin) 

#Library containing all the tests to run
add_subdirectory(test)

#Locate the openMP library to parallelize the code
find_package(OpenMP REQUIRED)

#Create the main executable
add_executable(simucell3d main.cpp)

#Indicate to main the path to project source directory
add_definitions(-DPROJECT_SOURCE_DIR="${PROJECT_SOURCE_DIR}")

#Link the lib, include and source libraries with the main library
target_link_libraries(simucell3d PUBLIC 
    src
    OpenMP::OpenMP_CXX
)

#In debog mode, we do not use openMP
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    message("-- Building in debug mode")
else()
    message("-- Building in release mode")
endif()





