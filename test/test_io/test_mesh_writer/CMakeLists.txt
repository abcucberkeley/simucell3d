cmake_minimum_required(VERSION 3.0)

#Use CPP 17
set(CMAKE_CXX_STANDARD 17)

#Indicate to main the path to project source directory
add_definitions(-DPROJECT_SOURCE_DIR="${PROJECT_SOURCE_DIR}")


#Compile the files containing the tests to run
add_executable(test_mesh_writer test_mesh_writer.cpp)

target_link_libraries(test_mesh_writer PUBLIC   src
                                                mesh
                                                io
)

#Add all the header files required to run the root file of the project
target_include_directories(test_mesh_writer PUBLIC "${PROJECT_SOURCE_DIR}/include/")
target_include_directories(test_mesh_writer PUBLIC "${PROJECT_SOURCE_DIR}/include/mesh")
target_include_directories(test_mesh_writer PUBLIC "${PROJECT_SOURCE_DIR}/include/io")

#Run the test
add_test(NAME io_mesh_writer_get_node_pos              COMMAND test_mesh_writer write_point_data_test)
add_test(NAME io_mesh_writer_write_cell_data_from_mesh COMMAND test_mesh_writer write_cell_data_from_mesh_test)
add_test(NAME io_mesh_writer_write_cell_data_from_cell COMMAND test_mesh_writer write_cell_data_from_cell_test)










