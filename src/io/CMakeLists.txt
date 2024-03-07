#The version of cmake required
cmake_minimum_required(VERSION 3.0)

FILE(GLOB SRC "*.cpp")

#Create the source library
add_library(io ${SRC})

#Locate the openMP library to parallelize the code
find_package(OpenMP REQUIRED)

#Indicate to the modules the path to project source directory
add_definitions(-DPROJECT_SOURCE_DIR="${PROJECT_SOURCE_DIR}")

target_include_directories(io PUBLIC    "${PROJECT_SOURCE_DIR}/include/"
                                        "${PROJECT_SOURCE_DIR}/include/mesh" 
                                        "${PROJECT_SOURCE_DIR}/include/io"  
                                        "${PROJECT_SOURCE_DIR}/include/triangulation_modules"  
                                        "${PROJECT_SOURCE_DIR}/include/math_modules"                                                             
)

                                   
target_link_libraries(io PUBLIC 
    mesh
    tinyxml2    
    OpenMP::OpenMP_CXX  
)
