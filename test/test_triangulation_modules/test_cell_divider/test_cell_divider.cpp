#include <cassert>
#include <string>
#include <iostream>
#include <vector>

#include "utils.hpp"
#include "mesh_writer.hpp"
#include "mesh_reader.hpp"

#include "custom_structures.hpp"
#include "cell_divider.hpp"


int add_intersection_point_test(){


    std::vector<double> cell_node_pos_lst{
        0,0,0,  1,0,0,  1,0,1,
        0,0,1,  0,1,0,  1,1,0,
        0,1,1,  1,1,1
    };



    std::vector<unsigned> cell_face_connectivity{
        0, 1, 3, //0
        2, 3, 1, //1
        0, 4, 1, //2
        5, 1, 4, //3
        0, 3, 4, //4
        6, 4, 3, //5
        1, 5, 2, //6
        7, 2, 5, //7
        5, 4, 7, //8
        6, 7, 4, //9
        3, 2, 6, //10
        7, 6, 2  //11
    };

    cell_ptr c = std::make_shared<cell>(cell_node_pos_lst, cell_face_connectivity, 0);
    
    //Initialize the cell 
    c->initialize_cell_properties();


    mesh_writer::write_face_data_file("./cell_before_division.vtk", {c});

    //Get the centroid of the cell
    vec3 centroid = c->compute_centroid();
    vec3 division_plane_normal(0,0,1);

    //Add the points where the plane intersects with the edges of the cell
    mesh m = cell_divider::add_intersection_points(c, centroid, division_plane_normal);

    //Make sure that all the faces have either 3 or 5 nodes 
    bool t1 = std::all_of(m.face_point_ids.begin(), m.face_point_ids.end(), [](const std::vector<unsigned>& face){return face.size() == 3 || face.size() == 5;});

    mesh_writer::write_cell_data_file("./cell_divider_test_mesh_1.vtk", {m});

    return !t1; 
}



int divide_faces_test(){
    std::vector<double> cell_node_pos_lst{
        0,0,0,  1,0,0,  1,0,1,
        0,0,1,  0,1,0,  1,1,0,
        0,1,1,  1,1,1
    };



    std::vector<unsigned> cell_face_connectivity{
        0, 1, 3, //0
        2, 3, 1, //1
        0, 4, 1, //2
        5, 1, 4, //3
        0, 3, 4, //4
        6, 4, 3, //5
        1, 5, 2, //6
        7, 2, 5, //7
        5, 4, 7, //8
        6, 7, 4, //9
        3, 2, 6, //10
        7, 6, 2  //11
    };

    cell_ptr c = std::make_shared<cell>(cell_node_pos_lst, cell_face_connectivity, 0);
    
    //Initialize the cell 
    c->initialize_cell_properties();

    //Get the centroid of the cell
    vec3 centroid = c->compute_centroid();
    vec3 division_plane_normal(0,0,1);

    //All the points added during the division process will have an id greater than this value
    const unsigned intersection_point_ids_threshold = c->get_node_lst().size();

    //Add the points where the plane intersects with the edges of the cell
    mesh m = cell_divider::add_intersection_points(c, centroid, division_plane_normal);

    //Retriangulate all the faces
    cell_divider::divide_faces(m, intersection_point_ids_threshold);

    //Make sure that all the faces have either 3 or 5 nodes 
    //bool t1 = std::all_of(m.face_point_ids.begin(), m.face_point_ids.end(), [](const std::vector<unsigned>& face){return face.size() == 3 || face.size() == 5;});

    mesh_writer::write_cell_data_file("./cell_divider_test_mesh_2.vtk", {m});

    return !std::all_of(m.face_point_ids.begin(), m.face_point_ids.end(), [](const std::vector<unsigned>& face){return face.size() == 3;});
}
//---------------------------------------------------------------------------------------------------------



//---------------------------------------------------------------------------------------------------------
int add_division_interface_test(){

    std::vector<double> cell_node_pos_lst{
        0,0,0,  1,0,0,  1,0,1,
        0,0,1,  0,1,0,  1,1,0,
        0,1,1,  1,1,1
    };
    
    std::vector<unsigned> cell_face_connectivity{
        0, 1, 3, //0
        2, 3, 1, //1
        0, 4, 1, //2
        5, 1, 4, //3
        0, 3, 4, //4
        6, 4, 3, //5
        1, 5, 2, //6
        7, 2, 5, //7
        5, 4, 7, //8
        6, 7, 4, //9
        3, 2, 6, //10
        7, 6, 2  //11
    };

    cell_ptr c = std::make_shared<cell>(cell_node_pos_lst, cell_face_connectivity, 0);
    
    //Initialize the cell 
    c->initialize_cell_properties();

    //Get the centroid of the cell
    vec3 centroid = c->compute_centroid();
    vec3 division_plane_normal(0,0,1);

    //All the points added during the division process will have an id greater than this value
    const unsigned intersection_point_ids_threshold = c->get_node_lst().size();

    //Add the points where the plane intersects with the edges of the cell
    mesh m = cell_divider::add_intersection_points(c, centroid, division_plane_normal);

    //Retriangulate all the faces
    cell_divider::divide_faces(m, intersection_point_ids_threshold);
    
    //Form a polygon by connecting all the points where the plane intersects with the edges of the cell
    initial_triangulation::coarse_triangulation(m);

    mesh_writer::write_cell_data_file("./cell_divider_test_mesh_3.vtk", {m});

    //Make sure that all the faces are triangles
    return !std::all_of(m.face_point_ids.begin(), m.face_point_ids.end(), [](const std::vector<unsigned>& face){return face.size() == 3;});

}
//---------------------------------------------------------------------------------------------------------


//---------------------------------------------------------------------------------------------------------
int map_points_to_xy_plane_test(){

     std::vector<double> cell_node_pos_lst{
        0,0,0,  1,0,0,  1,0,1,
        0,0,1,  0,1,0,  1,1,0,
        0,1,1,  1,1,1
    };
    
    std::vector<unsigned> cell_face_connectivity{
        0, 1, 3, //0
        2, 3, 1, //1
        0, 4, 1, //2
        5, 1, 4, //3
        0, 3, 4, //4
        6, 4, 3, //5
        1, 5, 2, //6
        7, 2, 5, //7
        5, 4, 7, //8
        6, 7, 4, //9
        3, 2, 6, //10
        7, 6, 2  //11
    };

    cell_ptr c = std::make_shared<cell>(cell_node_pos_lst, cell_face_connectivity, 0);
    
    //Initialize the cell 
    c->initialize_cell_properties();

    //Get the centroid of the cell
    vec3 centroid = c->compute_centroid();
    vec3 division_plane_normal = vec3(0,0.1,0.9).normalize();

    //All the points added during the division process will have an id greater than this value
    const unsigned intersection_point_ids_threshold = c->get_node_lst().size();

    //Add the points where the plane intersects with the edges of the cell
    mesh m = cell_divider::add_intersection_points(c, centroid, division_plane_normal);

    //Retriangulate all the faces
    cell_divider::divide_faces(m, intersection_point_ids_threshold);

    //Form a polygon by connecting all the points where the plane intersects with the edges of the cell
    initial_triangulation::coarse_triangulation(m);

    //Map the points of the division face to the xy plane
    auto [translation, rotation] = cell_divider::map_points_to_xy_plane(m, intersection_point_ids_threshold, division_plane_normal);
    mesh_writer::write_cell_data_file("./cell_divider_test_mesh_4.vtk", {m});


    //Check that all the nodes that should have been mapped to the xy planes have been mapped correctly
    for(unsigned p_id = intersection_point_ids_threshold; p_id < m.node_pos_lst.size() / 3; p_id++){
        if(!almost_equal(m.node_pos_lst[p_id * 3 + 2], 0.)){
            std::cout << "ERROR: A node that should have been mapped to the xy plane has not been mapped correctly" << std::endl;
            return 1;
        }
    }
    return 0;
}

//---------------------------------------------------------------------------------------------------------




//---------------------------------------------------------------------------------------------------------
int triangulate_division_interface_test(){

    //Get the mesh of a cell to divide
    mesh_reader reader(std::string(PROJECT_SOURCE_DIR) + std::string("/test/test_triangulation_modules/test_cell_divider/test_cell.vtk"));
    std::vector<mesh> cell_mesh_lst = reader.read();

    if (cell_mesh_lst.size() != 1){
        std::cout << "ERROR: The test cell should contain only one cell" << std::endl;
        return 1;
    }


    //Transfom the mesh into a cell
    cell_ptr c = std::make_shared<cell>(cell_mesh_lst[0], 0);

    //Initialize the cell
    c->initialize_cell_properties();

    //Get the origin and the normal of the division plane
    vec3 centroid = c->compute_centroid();
    vec3 division_plane_normal = vec3(0,0.1,0.9).normalize();

    //All the points that will be part of the diivsion interface will have an id greater than this value
    const unsigned division_point_ids_threshold = c->get_node_lst().size();

    //Add the points where the plane intersects with the edges of the cell
    mesh m = cell_divider::add_intersection_points(c, centroid, division_plane_normal);
    mesh_writer::write_cell_data_file("./cell_divider_test_mesh_3.vtk", {m});

    //At this stage some faces are not triangles anymore and must be triangulated before the division can be performed
    cell_divider::divide_faces(m, division_point_ids_threshold);
    mesh_writer::write_cell_data_file("./cell_divider_test_mesh_4.vtk", {m});

    //All the faces that will be part of the division interface will have an id greater than this value
    const unsigned division_face_ids_threshold  = m.face_point_ids.size();

    //Form a polygon by connecting all the points where the plane intersects with the edges of the cell
    // and then divide this polygon in triangles
    const unsigned nb_division_face_points = m.node_pos_lst.size() / 3 - division_point_ids_threshold;
    std::vector<unsigned> division_face_point_ids(nb_division_face_points);
    std::iota(division_face_point_ids.begin(), division_face_point_ids.end(), division_point_ids_threshold);
    m.face_point_ids.push_back(division_face_point_ids);
    initial_triangulation::coarse_triangulation(m);
    mesh_writer::write_cell_data_file("./cell_divider_test_mesh_5.vtk", {m});

    //Map the points of the division face to the xy plane. This is done such that we can use the 2D
    //Delaunay triangulation and it also facilitates subsequent tests
    auto [translation, rotation] = cell_divider::map_points_to_xy_plane(m, division_point_ids_threshold, division_plane_normal);
    mesh_writer::write_cell_data_file("./cell_divider_test_mesh_6.vtk", {m});

    //Discard the triangles of the division interface and replace them with new triangles
    cell_divider::triangulate_division_interface(
       4e-7, 
       m, 
       division_point_ids_threshold, 
       division_face_ids_threshold, 
       division_plane_normal
    );

    mesh_writer::write_cell_data_file("./cell_divider_test_mesh_7.vtk", {m});

    //Map the points of the division face back to the original position
    cell_divider::map_points_to_division_plane(m, division_point_ids_threshold, translation, rotation);

    mesh_writer::write_cell_data_file("./cell_divider_test_mesh_8.vtk", {m});

    return !std::all_of(m.face_point_ids.begin(), m.face_point_ids.end(), [](const std::vector<unsigned>& face_point_ids){return face_point_ids.size() == 3;});
}

//---------------------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------------------
//Create a dummy face type parameter set
face_type_parameters create_face_type_parameters(){
    face_type_parameters face_parameters;
    face_parameters.name_ = "dummy_face_type";
    face_parameters.face_type_global_id_ = 0;
    face_parameters.surface_tension_ = 0.1;
    face_parameters.adherence_strength_ = 1.;
    face_parameters.repulsion_strength_ = 1.;
    return face_parameters;
}

//Create a dummy cell type parameter set
std::shared_ptr<cell_type_parameters> create_cell_type_parameters(){
    std::shared_ptr<cell_type_parameters>  cell_parameters = std::make_shared<cell_type_parameters>();
    cell_parameters->name_ = "dummy_cell_type";
    cell_parameters->global_type_id_ = 0;

    //Add a dummy face type to the cell type
    cell_parameters->face_types_.push_back(create_face_type_parameters());
    return cell_parameters;
}
//---------------------------------------------------------------------------------------------------------


//---------------------------------------------------------------------------------------------------------
int create_daughter_cell_meshes_test(){

    //Get the mesh of a cell to divide
    mesh_reader reader(std::string(PROJECT_SOURCE_DIR) + std::string("/test/test_triangulation_modules/test_cell_divider/test_cell.vtk"));
    std::vector<mesh> cell_mesh_lst = reader.read();

    if (cell_mesh_lst.size() != 1){
        std::cout << "ERROR: The test cell should contain only one cell" << std::endl;
        return 1;
    }

    auto cell_parameters = create_cell_type_parameters();

    cell_ptr c = std::make_shared<epithelial_cell>(cell_mesh_lst[0], 0, cell_parameters);
    c->initialize_cell_properties();


    //Get the origin and the normal of the division plane
    vec3 centroid = c->compute_centroid();
    vec3 division_plane_normal = vec3(0,0.1,0.9).normalize();


    //All the points that will be part of the division interface will have ids >= than this value
    const unsigned division_point_ids_threshold = c->get_node_lst().size();


    //Add the points where the plane intersects with the edges of the cell
    mesh m = cell_divider::add_intersection_points(c, centroid, division_plane_normal);


    //At this stage some faces are not triangles anymore and must be triangulated before the division can be performed
    cell_divider::divide_faces(m, division_point_ids_threshold);


    //All the faces that will be part of the division interface will have an id greater than this value
    const unsigned division_face_ids_threshold  = m.face_point_ids.size();



    //Form a polygon by connecting all the points where the plane intersects with the edges of the cell
    // and then divide this polygon in triangles
    const unsigned nb_division_face_points = m.node_pos_lst.size() / 3 - division_point_ids_threshold;
    std::vector<unsigned> division_face_point_ids(nb_division_face_points);
    std::iota(division_face_point_ids.begin(), division_face_point_ids.end(), division_point_ids_threshold);
    m.face_point_ids.push_back(division_face_point_ids);
    initial_triangulation::coarse_triangulation(m);
 

    //Map the points of the division face to the xy plane. This is done such that we can use the 2D
    //Delaunay triangulation and it also facilitates subsequent tests
    auto [translation, rotation] = cell_divider::map_points_to_xy_plane(m, division_point_ids_threshold, division_plane_normal);


    //Triangulate the the division interface in a way that respects the constraints on the edge lengths
    cell_divider::triangulate_division_interface(
       4e-7, 
       m, 
       division_point_ids_threshold, 
       division_face_ids_threshold, 
       division_plane_normal
    );


    //Map the points of the division interface back to their original positions
    cell_divider::map_points_to_division_plane(m, division_point_ids_threshold, translation, rotation);

    //Split the mesh of the mother cell in 2 and create the daughter cells
    auto [daughter_cell_1, daughter_cell_2] = cell_divider::create_daughter_cells(
        c, 
        m, 
        division_point_ids_threshold, 
        division_face_ids_threshold, 
        division_plane_normal,
        centroid
    );


    mesh_writer::write_cell_data_file("./cell_divider_test_mesh_9.vtk",  {daughter_cell_1});
    mesh_writer::write_face_data_file("./cell_divider_test_mesh_10.vtk", {daughter_cell_2});


 
   return 0;
}



//---------------------------------------------------------------------------------------------------------
// The main function
int main (int argc, char** argv){

    //Check that the command line input is correctly formatted
    assert(argc == 2); 

    //Get the name of the test to run
    std::string test_name = argv[1];
    
    if (test_name == "add_intersection_point_test")             return add_intersection_point_test();
    if (test_name == "divide_faces_test")                       return divide_faces_test();
    if (test_name == "add_division_interface_test")             return add_division_interface_test();
    if (test_name == "map_points_to_xy_plane_test")             return map_points_to_xy_plane_test();
    if (test_name == "triangulate_division_interface_test")     return triangulate_division_interface_test();
    if (test_name == "create_daughter_cell_meshes_test")        return create_daughter_cell_meshes_test();


    




    

    std::cout << "TEST NAME :" << test_name << " DOES NOT EXIST" << std::endl;
    return 1;

}
//---------------------------------------------------------------------------------------------------------
