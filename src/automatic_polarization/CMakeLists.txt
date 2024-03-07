#The verscontact_modelsn of cmake required
cmake_minimum_required(VERSION 3.0)

FILE(GLOB SRC "*.cpp")

#Create the source library
add_library(automatic_polarization ${SRC})

#Locate the openMP library to parallelize the code
find_package(OpenMP REQUIRED)

target_include_directories(automatic_polarization PUBLIC    
    "${PROJECT_SOURCE_DIR}/include/"
    "${PROJECT_SOURCE_DIR}/include/math_modules"  
    "${PROJECT_SOURCE_DIR}/include/mesh" 
    "${PROJECT_SOURCE_DIR}/include/uspg"  
    "${PROJECT_SOURCE_DIR}/include/io"
    "${PROJECT_SOURCE_DIR}/include/automatic_polarization"                                                                                            
)
                                   
target_link_libraries(automatic_polarization PUBLIC 
    OpenMP::OpenMP_CXX     
)
