#The version of cmake required
cmake_minimum_required(VERSION 3.0)

FILE(GLOB SRC "*.cpp")

#Create the source library
add_library(math_modules ${SRC})


target_include_directories(math_modules PUBLIC  "${PROJECT_SOURCE_DIR}/include/"
                                                "${PROJECT_SOURCE_DIR}/include/math_modules"                                             
)
