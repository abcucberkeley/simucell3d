
#include "test_bpa.hpp"


//---------------------------------------------------------------------------------------------------------
int bpa_tester::compute_oriented_normal_test() const{

    //Create 3 oriented points with normals pointing in the same direction
    oriented_point p1(vec3(0., 0., 0.), vec3(0., 0., 1.));
    oriented_point p2(vec3(1., 0., 0.), vec3(0., 0., 1.));
    oriented_point p3(vec3(0., 1., 0.), vec3(0., 0., 1.));
    
    //Get the normal of the plane formed by the 3 points. This 
    //normal should point in the same direction as the average normal of the points
    vec3 plane_normal_1 = ball_pivoting_algorithm::compute_oriented_normal(p1, p2 ,p3);
    bool t1 = plane_normal_1 == vec3(0., 0., 1.);

    //Make the normal of the plane point in the opposite direction
    vec3 plane_normal_2 = ball_pivoting_algorithm::compute_oriented_normal(p1, p3 ,p2);
    bool t2 = plane_normal_2 == vec3(0., 0., 1.);


    return !(t1 && t2); 
}
//---------------------------------------------------------------------------------------------------------



//---------------------------------------------------------------------------------------------------------
int bpa_tester::points_coherently_oriented_test() const{
    
    //Create 3 oriented points with normals pointing in the same direction
    oriented_point p1(vec3(0., 0., 0.), vec3(0., 0., 1.));
    oriented_point p2(vec3(1., 0., 0.), vec3(0., 0., 1.));
    oriented_point p3(vec3(0., 1., 0.), vec3(0., 0., 1.));
    bool t1 = ball_pivoting_algorithm::points_coherently_oriented(p1, p2, p3);

    //A case where the points are not coherently oriented
    oriented_point p4(vec3(0., 1., 0.), vec3(0., 0., -1.));
    bool t2 = ball_pivoting_algorithm::points_coherently_oriented(p1, p2, p4) == false;

    return !(t1 && t2); 
}
//---------------------------------------------------------------------------------------------------------



//---------------------------------------------------------------------------------------------------------
//Make sure that the normal of the face {f1_a, f1_b, p} points in the same direction as the normal
//of the face {f1_a, f1_b, f1_c}
int bpa_tester::node_is_compatible_with_edge_test() const{
    
    //Create 4 nodes
    oriented_point f1_a(vec3(1., 0., 0.), vec3(0., 0., 1.));
    oriented_point f1_b(vec3(0., 1., 0.), vec3(0., 0., 1.));
    oriented_point f1_c(vec3(0., 0., 0.), vec3(0., 0., 1.));
    oriented_point p_1( vec3(1., 1., 0.), vec3(0., 0., 1.));

    //In this first case the normals of the face {f1_a, f1_b, p_1} and {f1_a, f1_b, f1_c} should be pointing in the same direction
    //Therefore the point p_1 should be compatible with the edge {f1_a, f1_b}
    bool t1 = ball_pivoting_algorithm::node_is_compatible_with_edge(p_1, f1_a, f1_b, f1_c) == true;

    //In this second case the normals of the face {f1_a, f1_b, p_2} and {f1_a, f1_b, f1_c} should be not
    //be pointing in the same direction
    oriented_point p_2(vec3(0.5, 0.5, 0.), vec3(0., 0., 1.));
    bool t2 = ball_pivoting_algorithm::node_is_compatible_with_edge(p_2, f1_a, f1_b, f1_c) == false;

    return !(t1 && t2);
}
//---------------------------------------------------------------------------------------------------------



//---------------------------------------------------------------------------------------------------------
int bpa_tester::compute_angle_between_circumspheres_test() const{

    //Create an edge along the y axis
    const oriented_point p1(vec3(0., -1., 0.), vec3());
    const oriented_point p2(vec3(0.,  1., 0.), vec3());
    const vec3 edge_center_vec = (p2.position_ - p1.position_) / 2. + p1.position_;

    //The 2 circumsphere centers in this example lie in the x-z plane
    const vec3 face_a_circum_center(-1., 0., 1.);
    const vec3 face_b_circum_center( 1., 0., 1.);

    //The angle between the 2 circumspheres should be 90 degrees
    const double angle_1 = ball_pivoting_algorithm::compute_angle_between_circumspheres(face_a_circum_center, face_b_circum_center, edge_center_vec);
    bool t1 = almost_equal(angle_1, M_PI / 2.);


    //The angle between the 2 circumspheres should be 180 degrees
    const vec3 face_c_circum_center(-1., 0., 0.);
    const vec3 face_d_circum_center( 1., 0., 0.);
    const double angle_2 = ball_pivoting_algorithm::compute_angle_between_circumspheres(face_c_circum_center, face_d_circum_center, edge_center_vec);
    bool t2 = almost_equal(std::abs(angle_2), M_PI);

    return !(t1 && t2);
}
//---------------------------------------------------------------------------------------------------------


//---------------------------------------------------------------------------------------------------------                  
int bpa_tester::compute_circum_sphere_center_test() const{

    //Create an instance of the ball pivoting algorithm
    ball_pivoting_algorithm bpa(1.);

    //Compute the circumsphere values returned by the function against known values. 
    oriented_point p1(vec3(0., 0., 0.), vec3(0., 0., 1.));
    oriented_point p2(vec3(1., 0., 0.), vec3(0., 0., 1.));
    oriented_point p3(vec3(0., 1., 0.), vec3(0., 0., 1.));
    std::optional<vec3> center_1 = bpa.compute_circum_sphere_center(p1, p2, p3);
    if (!center_1) return 1;

    oriented_point p4(vec3(0., 0., 0.),  vec3(0., 0., 1.));
    oriented_point p5(vec3(-1., 0., 0.), vec3(0., 0., 1.));
    oriented_point p6(vec3(0., -1., 0.), vec3(0., 0., 1.));
    std::optional<vec3> center_2 = bpa.compute_circum_sphere_center(p4, p5, p6);
    if (!center_2) return 1;


    oriented_point p7(vec3(0., 0., 0.), vec3(0., 0., 1.));
    oriented_point p8(vec3(2., 0., 0.), vec3(0., 0., 1.));
    oriented_point p9(vec3(0., 2., 0.), vec3(0., 0., 1.));
    std::optional<vec3> center_3 = bpa.compute_circum_sphere_center(p7, p8, p9);
    if (center_3) return 1;

    //Check that the computed values are correct
    vec3 vec_center_1 = center_1.value();
    vec3 vec_center_2 = center_2.value();

    bool t1 = (vec_center_1 - vec3( 0.5 ,  0.5 , 0.969536)).norm() < 1e-3;
    bool t2 = (vec_center_2 - vec3(-0.5 , -0.5 , 0.969536)).norm() < 1e-3;

    return !(t1 && t2);
}
//---------------------------------------------------------------------------------------------------------



//---------------------------------------------------------------------------------------------------------
//Unit test of the function ball_pivoting_algorithm::circum_sphere_is_empty()
int bpa_tester::circum_sphere_is_empty_test() const{

    //Create a few random points
    oriented_point p1(0, vec3(0., 0., 0.), vec3(0., 0., 1.));
    oriented_point p2(1, vec3(1., 0., 0.), vec3(0., 0., 1.));
    oriented_point p3(2, vec3(0., 1., 0.), vec3(0., 0., 1.));
    oriented_point p4(3, vec3(2., 2., 0.), vec3(0., 0., 1.));
    std::vector<oriented_point> points_1 = {p1, p2, p3, p4};

    //Create an instance of the ball pivoting algorithm
    ball_pivoting_algorithm bpa_1(
        points_1,        //The point cloud
        1.,             //The min edge length l_min
        0., 0., 0.,     //The min x, y, z coordinates
        2., 2., 1.      //The max x, y, z coordinates
    );

    //Compute the circumcenter of a random face
    const auto circum_center_opt_1 = bpa_1.compute_circum_sphere_center(p1, p2, p3);
    if(!circum_center_opt_1) return 1;
    const vec3 circum_center_1 = circum_center_opt_1.value();

    //Check that their is no point in the circum sphere except the 3 making the 
    //up the triangular face
    bool t1 = bpa_1.circum_sphere_is_empty(circum_center_1, p1, p2, p3);

    //Reproduce the same test with a different set of points
    oriented_point p5(5, vec3(0.5, 0.5, 0.), vec3(0., 0., 1.));
    std::vector<oriented_point> points_2 = {p1, p2, p3, p5};

    //Create an instance of the ball pivoting algorithm
    ball_pivoting_algorithm bpa_2(
        points_2,       //The point cloud
        1.,             //The min edge length l_min
        0., 0., 0.,     //The min x, y, z coordinates
        2., 2., 1.      //The max x, y, z coordinates
    );


    const auto circum_center_opt_2 = bpa_2.compute_circum_sphere_center(p1, p2, p3);
    if(!circum_center_opt_2) return 1;
    const vec3 circum_center_2 = circum_center_opt_2.value();
    bool t2 = bpa_2.circum_sphere_is_empty(circum_center_2, p1, p2, p3) == false;

    return !(t1 && t2);
}
//---------------------------------------------------------------------------------------------------------



//---------------------------------------------------------------------------------------------------------
//Unit test of the function ball_pivoting_algorithm::get_edge()
int bpa_tester::get_edge_test() const{

    //Create an instance of the bpa
    ball_pivoting_algorithm bpa(1.);

    //Create some random edges
    bpa.edge_lst_.emplace_back(0, 1);
    bpa.edge_lst_.emplace_back(1, 2);
    bpa.edge_lst_.emplace_back(2, 0);
    bpa.edge_lst_.emplace_back(0, 3);
    bpa.edge_lst_.emplace_back(3, 1);

    //Try to get an edge that exists
    const unsigned edge_id_0  = bpa.get_edge(0, 1);
    const unsigned edge_id_2  = bpa.get_edge(0, 2);

    //Try to get an edge that does not exist
    const unsigned edge_id_3 = bpa.get_edge(3, 2);

    bool t1 = edge_id_0 == 0;
    bool t2 = edge_id_2 == 2;

    bool t3 = edge_id_3 == bpa.edge_lst_.size() -1;
    bool t4 = bpa.edge_lst_.size() == 6;


    return !(t1 && t2 && t3 && t4);
}
//---------------------------------------------------------------------------------------------------------




//---------------------------------------------------------------------------------------------------------
int bpa_tester::node_is_inner_vertex_test() const{

    //Create an instance of the bpa
    ball_pivoting_algorithm bpa_1(1.);

    //Create a few random points
    oriented_point p0(0, vec3(0., 0., 0.), vec3(0., 0., 1.));
    oriented_point p1(1, vec3(1., 0., 0.), vec3(0., 0., 1.));
    oriented_point p2(2, vec3(0., 1., 0.), vec3(0., 0., 1.));
    oriented_point p3(3, vec3(0.5,0.5, 0.),vec3(0., 0., 1.));

    //Set the point cloud used by the bpa_1 
    bpa_1.node_lst_ = {p0, p1, p2, p3};

    //Create a few random faces
    bpa_1.face_lst_.emplace_back(0, 1, 3, 0);
    bpa_1.face_lst_.emplace_back(1, 2, 3, 1);
    bpa_1.face_lst_.emplace_back(0, 3, 2, 2);

    //Create the edges connected to the node 3
    bpa_1.edge_lst_.emplace_back(0, 3, 0, 2);
    bpa_1.edge_lst_.emplace_back(1, 3, 0, 1);
    bpa_1.edge_lst_.emplace_back(2, 3, 1, 2);

    //Indicate to node 3 that it is part of these edges
    bpa_1.node_edge_lst_[3].push_back(0);
    bpa_1.node_edge_lst_[3].push_back(1);
    bpa_1.node_edge_lst_[3].push_back(2);

    bpa_1.node_face_lst_[3].push_back(0);
    bpa_1.node_face_lst_[3].push_back(1);
    bpa_1.node_face_lst_[3].push_back(2);


    //Check that node 3 is an inner vertex
    bool t1 = bpa_1.node_is_inner_vertex(3);



    //Repeat the same test with a different set of points where the point should not be an inner vertex
    //Create an instance of the bpa
    ball_pivoting_algorithm bpa_2(1.);

    //Set the point cloud used by the bpa_2 
    bpa_2.node_lst_ = {p0, p1, p2, p3};

    //Create a few random faces
    bpa_2.face_lst_.emplace_back(0, 1, 3, 0);
    bpa_2.face_lst_.emplace_back(1, 2, 3, 1);
    bpa_2.face_lst_.emplace_back(0, 3, 2, 2);

    //Create the edges connected to the node 3
    bpa_2.edge_lst_.emplace_back(0, 3, 0, 2);
    bpa_2.edge_lst_.emplace_back(1, 3);
    bpa_2.edge_lst_.emplace_back(2, 3);

    //Indicate to node 3 that it is part of these edges
    bpa_2.node_edge_lst_[3].push_back(0);
    bpa_2.node_edge_lst_[3].push_back(1);
    bpa_2.node_edge_lst_[3].push_back(2);

    bpa_2.node_face_lst_[3].push_back(0);
    bpa_2.node_face_lst_[3].push_back(1);
    bpa_2.node_face_lst_[3].push_back(2);

    //Check that node 3 is an inner vertex
    bool t2 = bpa_2.node_is_inner_vertex(3) == false;
   
    return !(t1 && t2);
}
//---------------------------------------------------------------------------------------------------------



//---------------------------------------------------------------------------------------------------------
int bpa_tester::fill_surface_holes_test() const{
    /*
    Start by recreating the geometry shown below

       E____C___F
       |   /\   |  
       |  /  \  |   
       | /    \ |  
       A/______\B
       |\      /|   
       | \    / |  
       |  \  /  |   
       |___\/___|   
       H   D    G
    */


    //Create the nodes
    std::vector<double> node_pos_lst{
        0.,  0., 0.,    //Node A 0
        1.,  0., 0.,    //Node B 1
        0.5, 1., 0.,    //Node C 2
        0.5, -1., 0.,   //Node D 3
        0.,  1.,  0.,   //Node E 4      
        1.,  1.,  0.,   //Node F 5
        1., -1., 0.,    //Node G 6
        0., -1., 0.,    //Node H 7


        0.,  0., -1.,    //Node I 8
        1.,  0., -1.,    //Node J 9
        0.5, 1., -1.,    //Node K 10
        0.5, -1., -1.,   //Node L 11
        0.,  1.,  -1.,   //Node M 12      
        1.,  1.,  -1.,   //Node N 13
        1., -1., -1.,    //Node O 14
        0., -1., -1.,    //Node P 15
    };

    //Create the faces
    std::vector<std::vector<unsigned>> face_conn_lst{

        //The top 
//        {0, 1, 2},  //Face ABC 0
//        {0, 1, 3},  //Face ABD 1
        {0, 2, 4},  //Face ACE 2
        {1, 2, 5},  //Face BCF 3
        {0, 3, 7},  //Face ADH 4
        {1, 6, 3},  //Face BGD 5

        //The bottom
        {8, 9, 10},  //Face IJK 6
        {8, 9, 11},  //Face IJL 7
        {8, 10, 12},  //Face IKM 8
        {9, 10, 13},  //Face JKN 9
        {8, 11, 15},  //Face IPL 10
        {9, 14, 11},  //Face JLO 11

        //The sides
        {0, 8, 4 }, // Face AIE 
        {4, 12, 8}, // Face EMI
        {0, 8, 15}, // Face AIP
        {15, 7, 0}, // Face PHA

        {1, 9, 5 }, // Face BJF 
        {5, 13, 9}, // Face FNJ
        {1, 9, 14}, // Face BJO
        {14, 6, 1}, // Face OGB

        {4, 2, 10 }, // Face ECK
        {10, 12, 4}, // Face KME
        {5, 2, 10 }, // Face FCK
        {10, 13, 5}, // Face KNF

        {7, 3, 11 }, // Face HDL
        {11, 15, 7}, // Face LPH
        {6, 3, 11 }, // Face GDL
        {11, 14, 6}, // Face LOG

    };

    //Write the mesh in a file to check that its correct
    mesh m;
    m.node_pos_lst = node_pos_lst;
    m.face_point_ids = face_conn_lst;
    mesh_writer::write_cell_data_file(std::string(PROJECT_SOURCE_DIR) + "/test/test_triangulation_modules/test_bpa/test_fil_holes.vtk", {m});


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

    //Create the test object
    auto tester = bpa_tester();
    
    if (test_name == "points_coherently_oriented_test")             return tester.points_coherently_oriented_test();
    if (test_name == "compute_oriented_normal_test")                return tester.compute_oriented_normal_test();
    if (test_name == "node_is_compatible_with_edge_test")           return tester.node_is_compatible_with_edge_test();
    if (test_name == "compute_angle_between_circumspheres_test")    return tester.compute_angle_between_circumspheres_test();
    if (test_name == "compute_circum_sphere_center_test")           return tester.compute_circum_sphere_center_test();
    if (test_name == "circum_sphere_is_empty_test")                 return tester.circum_sphere_is_empty_test();
    if (test_name == "get_edge_test")                               return tester.get_edge_test();
    if (test_name == "node_is_inner_vertex_test")                   return tester.node_is_inner_vertex_test();
    if (test_name == "fill_surface_holes_test")                     return tester.fill_surface_holes_test();

    std::cout << "TEST NAME :" << test_name << " DOES NOT EXIST" << std::endl;
    return 1;

}
//---------------------------------------------------------------------------------------------------------
