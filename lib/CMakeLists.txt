#The version of cmake required
cmake_minimum_required(VERSION 3.0)

#Add the xml parser library
add_subdirectory(tinyxml2)
add_subdirectory(delaunator)

if(ENABLE_PYTHON_BINDINGS)
    add_subdirectory(pybind11)
endif()







