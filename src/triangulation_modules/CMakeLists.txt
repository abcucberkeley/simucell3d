#The version of cmake required
cmake_minimum_required(VERSION 3.0)

FILE(GLOB SRC "*.cpp")

#Indicate to main the path to project source directory
add_definitions(-DPROJECT_SOURCE_DIR="${PROJECT_SOURCE_DIR}")


#Create the source library
add_library(triangulation_modules ${SRC})


#Locate the openMP library to parallelize the code
find_package(OpenMP REQUIRED)


target_include_directories(triangulation_modules PUBLIC  
                                                "${PROJECT_SOURCE_DIR}/include/"
                                                "${PROJECT_SOURCE_DIR}/include/uspg"
                                                "${PROJECT_SOURCE_DIR}/include/io"
                                                "${PROJECT_SOURCE_DIR}/include/triangulation_modules"
                                                "${PROJECT_SOURCE_DIR}/include/mesh"
                                                "${PROJECT_SOURCE_DIR}/include/mesh/cell_types"
                                                "${PROJECT_SOURCE_DIR}/include/math_modules"
                                                "${PROJECT_SOURCE_DIR}/lib/delaunator/include"

)

target_link_libraries(triangulation_modules PUBLIC 
    uspg
    OpenMP::OpenMP_CXX    
)

