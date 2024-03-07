cmake_minimum_required(VERSION 3.0)

#Add the created libraries in the project
add_subdirectory(uspg)
add_subdirectory(math_modules)
add_subdirectory(mesh)
add_subdirectory(io)
add_subdirectory(triangulation_modules)
add_subdirectory(time_integration)
add_subdirectory(contact_models)
add_subdirectory(automatic_polarization)

if(ENABLE_PYTHON_BINDINGS)
    add_subdirectory(python_bindings)
endif()


FILE(GLOB SRC "*.cpp")

#Create the source library
add_library(src ${SRC})

#Add all the header files required to run the root file of the project
target_include_directories(src PUBLIC 
    "${PROJECT_SOURCE_DIR}/include/**"
    "${PROJECT_SOURCE_DIR}/include/io"
)

#Locate the openMP library to parallelize the code
find_package(OpenMP REQUIRED)

#Link root library of the src directory with all the other libraries created
target_link_libraries(src PUBLIC    
    uspg
    math_modules
    mesh
    io     
    triangulation_modules  
    time_integration  
    contact_models     
    automatic_polarization  
    OpenMP::OpenMP_CXX          
)




