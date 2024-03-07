#include <cassert>
#include <string>
#include <iostream>

#include "mesh_writer.hpp"
#include "mesh_reader.hpp"



//---------------------------------------------------------------------------------------------------------
//Create dummy mesh objects
std::vector<mesh> generate_arbitrary_mesh_lst(){

    //Create a mesh object
    mesh m1, m2; 

    //Create some arbitrary node positions
    m1.node_pos_lst = std::vector<double>{
        0,0,0,  1,0,0,  1,0,1,
        0,0,1,  0,1,0,  1,1,0,
        0,1,1,  1,1,1
    };
 

    //Create some arbitrary node positions
    m2.node_pos_lst = std::vector<double>{
        2,0,0,  3,0,0,  3,0,1,
        2,0,1,  2,1,0,  3,1,0,
        2,1,1,  3,1,1
    };


    //Create some arbitratry faces                         
    m1.face_point_ids = std::vector<std::vector<unsigned>>{
            {0, 1, 3},
            {2, 3, 1},
            {0, 4, 1},
            {5, 1, 4},
            {0, 3, 4},
            {6, 4, 3},
            {1, 5, 2},
            {7, 2, 5},
            {5, 4, 7},
            {6, 7, 4},
            {3, 2, 6},
            {7, 6, 2} 
    };
       
    m2.face_point_ids = std::vector<std::vector<unsigned>>(m1.face_point_ids);
    return {m1, m2};
}
//---------------------------------------------------------------------------------------------------------



//---------------------------------------------------------------------------------------------------------
//Create dummy cells
std::vector<cell_ptr> generate_arbitrary_cell_lst(){

    std::vector<double> cell_node_pos_1{
        0,0,0,  1,0,0,  1,0,1,
        0,0,1,  0,1,0,  1,1,0,
        0,1,1,  1,1,1
    };


    std::vector<double> cell_node_pos_2{
        2,0,0,  3,0,0,  3,0,1,
        2,0,1,  2,1,0,  3,1,0,
        2,1,1,  3,1,1
    };

    std::vector<unsigned> cell_face_connectivity{
        0, 1, 3,
        2, 3, 1,
        0, 4, 1,
        5, 1, 4,
        0, 3, 4,
        6, 4, 3,
        1, 5, 2,
        7, 2, 5,
        5, 4, 7,
        6, 7, 4,
        3, 2, 6,
        7, 6, 2 
    };

    cell_ptr c0 = std::make_shared<cell>(cell_node_pos_1, cell_face_connectivity, 0);
    cell_ptr c1 = std::make_shared<cell>(cell_node_pos_2, cell_face_connectivity, 1);

    return {c0, c1};
}
//---------------------------------------------------------------------------------------------------------




//---------------------------------------------------------------------------------------------------------
//Test that the point_data are correctly written in the file 
int write_point_data_test(){

    auto mesh_lst = generate_arbitrary_mesh_lst();
    auto m1 = mesh_lst.front();

    //Write the cell_data
    mesh_writer::write_cell_data_file(std::string(PROJECT_SOURCE_DIR) + "/test/test_io/test_mesh_writer/test_1.vtk", mesh_lst);

    //Use the mesh reader to check the coordinates
    mesh_reader reader(std::string(PROJECT_SOURCE_DIR) + "/test/test_io/test_mesh_writer/test_1.vtk");
    
    //Get the position of the nodes
    auto node_pos = reader.get_node_pos();

    //Check that the written positions are the correct ones
    return !std::equal(m1.node_pos_lst.begin(),  m1.node_pos_lst.end(), node_pos.begin());
}
//---------------------------------------------------------------------------------------------------------


//---------------------------------------------------------------------------------------------------------
//Write a mesh file where the input is a list of mesh objects
int write_cell_data_from_mesh_test(){

    auto mesh_lst = generate_arbitrary_mesh_lst();

    //Write the cell_data
    mesh_writer::write_cell_data_file(std::string(PROJECT_SOURCE_DIR) + "/test/test_io/test_mesh_writer/test_2.vtk", mesh_lst);

    //Use the mesh reader to check the coordinates
    mesh_reader reader(std::string(PROJECT_SOURCE_DIR) + "/test/test_io/test_mesh_writer/test_2.vtk");

    //Get the mesh object stored in the file
    std::vector<mesh> mesh_lst_written = reader.read();

    //Compare the generated meshes from the generate_arbitrary_mesh_lst() function with the meshes written in ./test_2.vtk
    mesh m1 = mesh_lst[0], m2 = mesh_lst[1];
    mesh m1_written = mesh_lst_written[0], m2_written = mesh_lst_written[1];

    //Compare the points of the true and written meshes
    bool t1 = std::equal(m1.node_pos_lst.begin(), m1.node_pos_lst.end(), m1_written.node_pos_lst.begin());
    bool t2 = std::equal(m2.node_pos_lst.begin(), m2.node_pos_lst.end(), m2_written.node_pos_lst.begin());

    //Compare the faces stored
    bool t3 = std::equal(m1.face_point_ids.begin(), m1.face_point_ids.end(), m1_written.face_point_ids.begin());
    bool t4 = std::equal(m2.face_point_ids.begin(), m2.face_point_ids.end(), m2_written.face_point_ids.begin());

    return !(t1 && t2 && t3 && t4);
}
//---------------------------------------------------------------------------------------------------------


//---------------------------------------------------------------------------------------------------------
int write_cell_data_from_cell_test(){

    //Generate a list of dummy cells
    auto cell_lst = generate_arbitrary_cell_lst();

    //Write the cell_data
    mesh_writer::write_cell_data_file(std::string(PROJECT_SOURCE_DIR) + "/test/test_io/test_mesh_writer/test_3.vtk", cell_lst);

    //Use the mesh reader to check the coordinates
    mesh_reader reader(std::string(PROJECT_SOURCE_DIR) + "/test/test_io/test_mesh_writer/test_3.vtk");

    //Read the data written in the VTK file
    std::vector<mesh> mesh_lst = reader.read();
    mesh m1 = mesh_lst[0], m2 = mesh_lst[1];

    //Check the coordinates of the cells
    std::vector<double> correct_cell_node_pos_1{
        0,0,0,  1,0,0,  1,0,1,
        0,0,1,  0,1,0,  1,1,0,
        0,1,1,  1,1,1
    };

    std::vector<double> correct_cell_node_pos_2{
        2,0,0,  3,0,0,  3,0,1,
        2,0,1,  2,1,0,  3,1,0,
        2,1,1,  3,1,1
    };

    bool t1 = std::equal(correct_cell_node_pos_1.begin(), correct_cell_node_pos_1.end(), m1.node_pos_lst.begin());
    bool t2 = std::equal(correct_cell_node_pos_2.begin(), correct_cell_node_pos_2.end(), m2.node_pos_lst.begin());

    //Check the connectivity of the faces
    std::vector<std::vector<unsigned>> correct_face_connectivity{
        {0, 1, 3},
        {2, 3, 1},
        {0, 4, 1},
        {5, 1, 4},
        {0, 3, 4},
        {6, 4, 3},
        {1, 5, 2},
        {7, 2, 5},
        {5, 4, 7},
        {6, 7, 4},
        {3, 2, 6},
        {7, 6, 2} 
    };

    bool t3 = std::equal(correct_face_connectivity.begin(), correct_face_connectivity.end(), m1.face_point_ids.begin());
    bool t4 = std::equal(correct_face_connectivity.begin(), correct_face_connectivity.end(), m2.face_point_ids.begin());

    return !(t1 && t2 && t3 && t4);
}
//---------------------------------------------------------------------------------------------------------





//---------------------------------------------------------------------------------------------------------
int write_face_data_test(){

    //Generate a list of dummy cells
    auto cell_lst = generate_arbitrary_cell_lst();

    //Write the cell_data
    mesh_writer::write_face_data_file(std::string(PROJECT_SOURCE_DIR) + "/test/test_io/test_mesh_writer/test_4.vtk", cell_lst);

    return 1;
}
//---------------------------------------------------------------------------------------------------------





//---------------------------------------------------------------------------------------------------------
// The main function
int main (int argc, char** argv){

    //Check that the command line input is correctly formatted
    assert(argc == 2); 

    //Get the name of the test to run
    std::string test_name = argv[1];
    
    //Run the selected test
    if (test_name == "write_point_data_test")           return write_point_data_test();
    if (test_name == "write_cell_data_from_mesh_test")  return write_cell_data_from_mesh_test();
    if (test_name == "write_cell_data_from_cell_test")  return write_cell_data_from_cell_test();
    if (test_name == "write_face_data_test")            return write_face_data_test();



    





    std::cout << "TEST NAME :" << test_name << " DOES NOT EXIST" << std::endl;
    return 1;
}
//---------------------------------------------------------------------------------------------------------



