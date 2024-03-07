#The version of cmake required
cmake_minimum_required(VERSION 3.0)

FILE(GLOB SRC "*.cpp")

#Create the source library
add_library(uspg ${SRC})


target_include_directories(uspg PUBLIC  "${PROJECT_SOURCE_DIR}/include/"
                                        "${PROJECT_SOURCE_DIR}/include/uspg"
                                        "${PROJECT_SOURCE_DIR}/include/math_modules"

)



