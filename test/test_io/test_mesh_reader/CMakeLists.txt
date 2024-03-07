cmake_minimum_required(VERSION 3.0)

#Use CPP 17
set(CMAKE_CXX_STANDARD 17)

#Indicate to main the path to project source directory
add_definitions(-DPROJECT_SOURCE_DIR="${PROJECT_SOURCE_DIR}")


#Compile the files containing the tests to run
add_executable(test_mesh_reader test_mesh_reader.cpp)

target_link_libraries(test_mesh_reader PUBLIC   src
                                                mesh
                                                io
)

#Add all the header files required to run the root file of the project
target_include_directories(test_mesh_reader PUBLIC "${PROJECT_SOURCE_DIR}/include/")
target_include_directories(test_mesh_reader PUBLIC "${PROJECT_SOURCE_DIR}/include/mesh")
target_include_directories(test_mesh_reader PUBLIC "${PROJECT_SOURCE_DIR}/include/io")

#Run the test
add_test(NAME io_mesh_reader_instantiation          COMMAND test_mesh_reader instantiation_test)
add_test(NAME io_mesh_reader_get_node_pos           COMMAND test_mesh_reader get_node_pos_test)
add_test(NAME io_mesh_reader_read_cell_faces         COMMAND test_mesh_reader read_cell_faces_test)
add_test(NAME io_mesh_reader_get_cell_mesh         COMMAND test_mesh_reader get_cell_mesh_test)
add_test(NAME io_mesh_reader_get_cell_types_test        COMMAND test_mesh_reader get_cell_types_test)







