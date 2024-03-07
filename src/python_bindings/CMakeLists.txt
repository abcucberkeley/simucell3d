cmake_minimum_required(VERSION 3.0)


FILE(GLOB SRC "*.cpp")

#Create the source library
#add_library(python_bindings ${SRC})


#Find the python library
#find_package(PythonLibs REQUIRED)

#execute_process(COMMAND python3 -c "import sysconfig; print(sysconfig.get_path('include'))" RESULT_VARIABLE PYTHON_INCLUDE_DIRS)
#execute_process(COMMAND python3 -c "import sysconfig; print(sysconfig.get_config_var('LIBDIR'))" RESULT_VARIABLE PYTHON_LIBRARY)

#find_package(Python REQUIRED COMPONENTS Interpreter Development)
#include_directories(${PYTHON_INCLUDE_DIRS})
message(STATUS "Python3_INCLUDE_DIRS: ${PYTHON_INCLUDE_DIRS}")


#Locate the openMP library to parallelize the code
find_package(OpenMP REQUIRED)


#Add the pybind11 library
pybind11_add_module(simucell3d_python_wrapper simucell3d_wrapper.cpp)

target_include_directories(simucell3d_python_wrapper PRIVATE 
    "${PROJECT_SOURCE_DIR}/include/**"
    "${PROJECT_SOURCE_DIR}/include/python_bindings"
    "${PROJECT_SOURCE_DIR}/include/io"
    "${PROJECT_SOURCE_DIR}/lib/pybind11/**"
    "${PROJECT_SOURCE_DIR}/lib/pybind11/include/"
    "${PYTHON_INCLUDE_DIRS}"
)


target_link_libraries(simucell3d_python_wrapper PRIVATE
                                    src    
                                    uspg
                                    math_modules
                                    mesh
                                    io     
                                    triangulation_modules  
                                    time_integration  
                                    contact_models       
                                    OpenMP::OpenMP_CXX          
)





