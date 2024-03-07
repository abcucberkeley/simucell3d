#include <cassert>
#include <string>
#include <iostream>
#include <vector>

#include "utils.hpp"
#include "mesh_writer.hpp"
#include "custom_structures.hpp"
#include "initial_triangulation.hpp"


int coarse_triangulation_test(){

    //Create the mesh of a cube
    mesh m;

    m.node_pos_lst = {
        0,0,0,
        1,0,0,
        1,1,0,
        0,1,0,
        0,0,1,
        1,0,1,
        1,1,1,
        0,1,1
    };

    m.face_point_ids = {
        {0,3,2,1},
        {4,5,6,7},
        {0,1,5,4},
        {1,2,6,5},
        {2,3,7,6},
        {3,0,4,7}
    };

    //Triangulate the mesh of the cube
    initial_triangulation::coarse_triangulation(m);

    //Make sure that the mesh has the expected number of nodes and faces
    bool t1 = m.get_nb_nodes() == 14;
    bool t2 = m.get_nb_faces() == 24;

    //Make sure that all the faces are triangles
    bool t3 = std::all_of(m.face_point_ids.begin(), m.face_point_ids.end(), [](const auto& f) -> bool {return f.size() == 3;});

    //It's hard to make sure that the coarse triangulation has been done properly, but what we can do is
    //check that the resulting triangulated mesh has the correct area and volume
    bool t4 = false, t5 = false;

    if(t3){
        //Convert the mesh to cell
        cell_ptr c = std::make_shared<cell>(m, 0);
        c->initialize_cell_properties();
        t4 = almost_equal(c->get_area(), 6.);
        t5 = almost_equal(c->get_volume(), 1.);
    }
    return !(t1 && t2 && t3 && t4 && t5);
}

int run_initial_triangulation(){

    //Create the mesh of a cube
    mesh m;

    m.node_pos_lst = {
        0,0,0,
        1,0,0,
        1,1,0,
        0,1,0,
        0,0,1,
        1,0,1,
        1,1,1,
        0,1,1
    };

    m.face_point_ids = {
        {0,3,2,1},
        {4,5,6,7},
        {0,1,5,4},
        {1,2,6,5},
        {2,3,7,6},
        {3,0,4,7}
    };

    mesh final_mesh  = initial_triangulation::triangulate_surface(0.05, 0.15, m, 0);

    //Write the mesh
    //Write the cell_data
    mesh_writer::write_cell_data_file("./test_bpa_2.vtk", {final_mesh});


    return 0;
}




//---------------------------------------------------------------------------------------------------------
// The main function
int main (int argc, char** argv){

    //Check that the command line input is correctly formatted
    assert(argc == 2); 

    //Get the name of the test to run
    std::string test_name = argv[1];
    
    if (test_name == "coarse_triangulation_test")     return coarse_triangulation_test();
    if (test_name == "run_initial_triangulation")     return run_initial_triangulation();



    

    std::cout << "TEST NAME :" << test_name << " DOES NOT EXIST" << std::endl;
    return 1;

}
//---------------------------------------------------------------------------------------------------------
