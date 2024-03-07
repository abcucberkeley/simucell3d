#include <cassert>
#include <string>
#include <iostream>

#include "mesh_reader.hpp"


/*
Contains all the tests run on the mesh_reader class
*/


//---------------------------------------------------------------------------------------------------------
//Test that mesh reader is correctly instantiated as a singleton
int instantiation_test(){

    //Get the cwd
    std::string source_dir_path = PROJECT_SOURCE_DIR;

    //This path does not exist
    std::string inexistant_mesh_path = source_dir_path + "/test/test_io/test_mesh_reader/test_meshes/inexistant_mesh_file.vtk";

    //This mesh file has the correct version 
    std::string valid_mesh_test_path   = source_dir_path + "/test/test_io/test_mesh_reader/test_meshes/valid_mesh_file.vtk";

    bool t1 = false, t2 = true;

    //Try to open the mesh reader on an unexistant file
    try{
        mesh_reader reader_1(inexistant_mesh_path);
    }
    catch(const mesh_reader_exception& e){
        t1 = true;
    }

    //Try to open the mesh reader with a valid file
    try{
        mesh_reader reader_3(valid_mesh_test_path);
    }
    catch(const mesh_reader_exception& e){
        t2 = false;
    }

    std::cout << "t1: " << t1 << std::endl;
    std::cout << "t2: " << t2 << std::endl;

    return !(t1 && t2);

}
//---------------------------------------------------------------------------------------------------------



//---------------------------------------------------------------------------------------------------------
//Test the method off the mesh reader that reads the node positions
int get_node_pos_test(){

    //Get the cwd
    std::string source_dir_path = PROJECT_SOURCE_DIR;

    //This path does not exist
    std::string invalid_nb_nodes_mesh = source_dir_path + "/test/test_io/test_mesh_reader/test_meshes/invalid_nb_nodes.vtk";

    //This mesh file does not have the correct version
    std::string invalid_node_format_mesh = source_dir_path + "/test/test_io/test_mesh_reader/test_meshes/invalid_node_format.vtk";

    //This mesh file has the correct version 
    std::string valid_mesh_test_path   = source_dir_path + "/test/test_io/test_mesh_reader/test_meshes/valid_mesh_file.vtk";

    bool t1 = false, t2 = false, t3 = true;

    mesh_reader reader_1(invalid_nb_nodes_mesh);
    try{

        //Try to load nodes with missing number of node information
        auto v1 = reader_1.get_node_pos();

    }
    catch(const mesh_reader_exception& e){
        t1 = true;
    }
    
    mesh_reader reader_2(invalid_node_format_mesh);
    try{
        //Try to load nodes with wrong format
        auto v2 = reader_2.get_node_pos();
    }
    catch(const mesh_reader_exception& e){
        t2 = true;
    }
    
    //Try to load valid node data
    mesh_reader reader_3(valid_mesh_test_path);
    std::vector<double> v3;
    try{
        v3 = reader_3.get_node_pos();

    }
    catch(const mesh_reader_exception& e){
        t3 = false;
    }


    //Check that the node loaded have the good coordinates
    std::vector<double> correct_coordinates{
        4.5e-07,    0,      0,      7.45e-06,   0,      0,      7.45e-06,   7e-06,  0,
        4.5e-07,    7e-06,  0,      4.5e-07,    0,      7e-06,  7.45e-06,   0,      7e-06,
        7.45e-06,   7e-06,  7e-06,  4.5e-07,    7e-06,  7e-06,  7.5e-06,    0,      0,
        1.45e-05,   0,      0,      1.45e-05,   7e-06,  0,      7.5e-06,    7e-06,  0,
        7.5e-06,    0,      7e-06,  1.45e-05,   0,      7e-06,  1.45e-05,   7e-06,  7e-06,  
        7.5e-06,    7e-06,  7e-06};
    
    bool t4 = false;
    if(t3){
        t4  = std::equal(correct_coordinates.begin(), correct_coordinates.end(), v3.begin());
    }


    return !(t1 && t2 && t3 && t4);
}


//---------------------------------------------------------------------------------------------------------
int read_cell_faces_test(){

    //Get the cwd
    std::string source_dir_path = PROJECT_SOURCE_DIR;

    //Test the mesh reader on files wih correct and incorrect formats
    std::string valid_mesh_test_path   = source_dir_path + "/test/test_io/test_mesh_reader/test_meshes/hexagonal_cell_layer.vtk";
    std::string invalid_data_1  = source_dir_path + "/test/test_io/test_mesh_reader/test_meshes/invalid_cell_type_1.vtk";
    std::string invalid_data_2   = source_dir_path + "/test/test_io/test_mesh_reader/test_meshes/invalid_cell_data_1.vtk";
    std::string invalid_data_3   = source_dir_path + "/test/test_io/test_mesh_reader/test_meshes/invalid_cell_data_2.vtk";


    bool t1 = true, t2 = false, t3 = false, t4 = false;

    //Try to load valid cell data
    mesh_reader reader_1(valid_mesh_test_path);

    std::vector<std::vector<unsigned>> cell_face_lst_lst_1;
    try{
        cell_face_lst_lst_1 = reader_1.read_cell_faces();
    }catch(const mesh_reader_exception& e){
        t1 = false;
    }

    mesh_reader reader_2(invalid_data_1);

    try{
        auto cell_face_lst_lst_2 = reader_2.read_cell_faces();
    }catch(const mesh_reader_exception& e){
        t2 = true;
    }


    mesh_reader reader_3(invalid_data_2);

    try{
        auto cell_face_lst_lst_3 = reader_3.read_cell_faces();
    }catch(const mesh_reader_exception& e){
        t3 = true;
    }


    mesh_reader reader_4(invalid_data_3);

    try{
        auto cell_face_lst_lst_4 = reader_4.read_cell_faces();
    }catch(const mesh_reader_exception& e){
        t4 = true;
    }


    //Check that the content which is read is correct
    bool t5 = false;
    if(t1){

        std::vector<unsigned> correct_data{
            8,6,0,1,2,3,4,5,6,6,7,8,9,10,11,4,0,1,7,6,4,1,2,8,7,4,2,3,9,8,4,3,4,10,9,4,4,5,11,10,4,5,0,6,11
        };

        auto loaded_data = cell_face_lst_lst_1[0];


        t5 = std::equal(correct_data.begin(), correct_data.end(), loaded_data.begin());
    }




    return !(t1 && t2 && t3 && t4);

}





//---------------------------------------------------------------------------------------------------------
//Test if the mesh reader correclt connverts the data stored in the input mesh file
//into meshes
int get_cell_mesh_test(){

    //Create some artificial node position
    std::vector<double> point_coord_pos{
        6.6e-06,0,2e-05,3.3e-06,5.7158e-06,2e-05,-3.3e-06,5.7158e-06,2e-05,
        -6.6e-06,0,2e-05,-3.3e-06,-5.7158e-06,2e-05,3.3e-06,-5.7158e-06,2e-05,
        6.6e-06,0,-2e-05,3.3e-06,5.7158e-06,-2e-05,-3.3e-06,5.7158e-06,-2e-05,
        -6.6e-06,0,-2e-05,-3.3e-06,-5.7158e-06,-2e-05,3.3e-06,-5.7158e-06,-2e-05,
        1.65e-05,5.7158e-06,2e-05,1.32e-05,1.14316e-05,2e-05,6.6e-06,1.14316e-05,2e-05,
        3.3e-06,5.7158e-06,2e-05,6.6e-06,0,2e-05,1.32e-05,0,2e-05
    };


    //Create some artificial face data
    std::vector<std::vector<unsigned>> face_connectivity_data{
        {8,6,0,1,2,3,4,5,6,6,7,8,9,10,11,4,0,1,7,6,4,1,2,8,7,4,2,3,9,8,4,3,4,10,9,4,4,5,11,10,4,5,0,6,11}, 
        {8,6,0,1,2,3,4,5,6,6,7,8,9,10,11,4,0,1,7,6,4,1,2,8,7,4,2,3,9,8,4,3,4,10,9,4,4,5,11,10,4,5,0,6,11}
    };

    //Create some artificial face data
    std::vector<mesh>  mesh_lst_1;
    try{
        mesh_lst_1  = mesh_reader::get_cell_mesh(point_coord_pos, face_connectivity_data);
    }
    catch(const mesh_reader_exception& e){
        return 1;
    }

    //Check that there is the correct number of meshes
    bool t1 = mesh_lst_1.size() == 2;

    //Get the 2 meshes and run tests on them
    mesh m1 = mesh_lst_1[0];
    mesh m2 = mesh_lst_1[1];

    //Check that the two meshes have the correct faces
    std::vector<std::vector<unsigned>> correct_face_lst{
        {0,1,2,3,4,5}, {6,7,8,9,10,11}, {0,1,7,6}, {1,2,8,7}, {2,3,9,8}, {3,4,10,9}, {4,5,11,10}, {5,0,6,11}
    };

    bool t3 = std::equal(correct_face_lst.begin(), correct_face_lst.end(), m1.face_point_ids.begin());
    bool t4 = std::equal(correct_face_lst.begin(), correct_face_lst.end(), m2.face_point_ids.begin());


    //Make sure the nodes of the cells have the correct coordinates
    std::vector<double> correct_node_pos(point_coord_pos.begin(), point_coord_pos.begin() +12*3);
    bool t5 = std::equal(correct_node_pos.begin(), correct_node_pos.end(), m1.node_pos_lst.begin());
    bool t6 = std::equal(correct_node_pos.begin(), correct_node_pos.end(), m2.node_pos_lst.begin());



    //One of the faces have to much nodes
    std::vector<std::vector<unsigned>> incorrect_face_connectivity_data_1{
        {8,6,0,1,2,3,4,5,6,6,7,8,9,10,11,4,0,1,7,6,4,1,2,8,7,4,2,3,9,8,4,3,4,10,9,4,4,5,11,10,4,5,0,6,11, 6}, 
        {8,6,0,1,2,3,4,5,6,6,7,8,9,10,11,4,0,1,7,6,4,1,2,8,7,4,2,3,9,8,4,3,4,10,9,4,4,5,11,10,4,5,0,6,11}
    };


    //Test what happens when one of thhe cell has more nodes than anticipated
    bool t7 = false;
    try{
        std::vector<mesh>  mesh_lst_2 = mesh_reader::get_cell_mesh(point_coord_pos, incorrect_face_connectivity_data_1);
    }catch(const mesh_reader_exception& e){t7 = true;}


    //One of the faces have to much nodes
    std::vector<std::vector<unsigned>> incorrect_face_connectivity_data_2{
        {8,6,0,1,2,3,4,5,6,6,7,8,9,10,11,4,0,1,7,6,4,1,2,8,7,4,2,3,9,8,4,3,4,10,9,4,4,5,11,10}, 
        {8,6,0,1,2,3,4,5,6,6,7,8,9,10,11,4,0,1,7,6,4,1,2,8,7,4,2,3,9,8,4,3,4,10,9,4,4,5,11,10,4,5,0,6,11}
    };


    //Test what happens when one of thhe cell has less faces than anticipated
    bool t8 = false;
    try{
        auto  mesh_lst_3 = mesh_reader::get_cell_mesh(point_coord_pos, incorrect_face_connectivity_data_2);
    }catch(const mesh_reader_exception& e){t8 = true;}

    return !(t1 && t3 && t4 && t5 && t6 && t7 && t8);
}
//---------------------------------------------------------------------------------------------------------




//---------------------------------------------------------------------------------------------------------
int get_cell_types_test(){


    //This mesh file has the correct version 
    std::string valid_mesh_test_path   = std::string(PROJECT_SOURCE_DIR) + "/test/test_io/test_mesh_reader/test_meshes/valid_mesh_file.vtk";


    //Try to open the mesh reader on an unexistant file
    mesh_reader reader(valid_mesh_test_path);

    //Read the cell type array
    std::vector<short> cell_type_lst = reader.get_cell_types();
    std::vector<short> correct_cell_types_lst{0, 1, 2, 3, 4};

    return ! std::equal(cell_type_lst.begin(), cell_type_lst.end(), correct_cell_types_lst.begin());
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
    if (test_name == "instantiation_test")      return instantiation_test();
    if (test_name == "get_node_pos_test")       return get_node_pos_test();
    if (test_name == "read_cell_faces_test")    return read_cell_faces_test();
    if (test_name == "get_cell_mesh_test")      return get_cell_mesh_test();
    if (test_name == "get_cell_types_test")     return get_cell_types_test();


    std::cout << "TEST NAME :" << test_name << " DOES NOT EXIST" << std::endl;
    return 1;
}
//---------------------------------------------------------------------------------------------------------



