#include "test_contact_model_abstract.hpp"



int tester_contact_model_abstract::compute_node_triangle_distance_test(){

    /*
            
       \  c |  
        \   |
         \  |
          \ |
           \|C
            |\
            | \
            |  \
       ac   |   \     bc
            |    \
            |     \
            | abc  \ 
            |       \
    ________|________\__________________
            |A       B\
      a     |   ab     \  b
            |           \
    

    The CPA of a point can be located in one of the following 7 regions:
    1) Region a
    2) Region b
    3) Region c
    4) Region ab
    5) Region bc
    6) Region ac
    7) Region abc

    //The following code will place nodes in the 7 different regions and check that the 
    CPA is correctly computed as well as the distance to the triangle
    */


   //Those are the 3 points of the triangle
    vec3 A(0., 0., 0.); 
    vec3 B(1., 0., 0.);
    vec3 C(0., 1., 0.);

    //First test the region a
    vec3 p_a1(-1., -1.,  0.);
    vec3 p_a2( 0.,  0.,  0.);

    const auto [p_a1_dist_sq, p_a1_cpa_bary] = contact_model_abstract::compute_node_triangle_distance(p_a1, A, B, C);
    const auto [p_a2_dist_sq, p_a2_cpa_bary] = contact_model_abstract::compute_node_triangle_distance(p_a2, A, B, C);

    bool t1 = p_a1_dist_sq == 2.;
    bool t2 = p_a1_cpa_bary == vec3(1., 0. , 0.);

    bool t3 = p_a2_dist_sq == 0.;
    bool t4 = p_a2_cpa_bary == vec3(1., 0. , 0.);

    std::cout << "t1 " << t1 << std::endl;
    std::cout << "t2 " << t2 << std::endl;
    std::cout << "t3 " << t3 << std::endl;
    std::cout << "t4 " << t4 << std::endl;

    //Test the region b
    vec3 p_b1( 2., -1.,  0.);
    vec3 p_b2( 1.,  0.,  0.);

    const auto [p_b1_dist_sq, p_b1_cpa_bary] = contact_model_abstract::compute_node_triangle_distance(p_b1, A, B, C);
    const auto [p_b2_dist_sq, p_b2_cpa_bary] = contact_model_abstract::compute_node_triangle_distance(p_b2, A, B, C);


    bool t5 = p_b1_dist_sq == 2.;
    bool t6 = p_b1_cpa_bary == vec3(0., 1. , 0.);

    bool t7 = p_b2_dist_sq == 0.;
    bool t8 = p_b2_cpa_bary == vec3(0., 1. , 0.);

    std::cout << "t5 " << t5 << std::endl;
    std::cout << "t6 " << t6 << std::endl;
    std::cout << "t7 " << t7 << std::endl;
    std::cout << "t8 " << t8 << std::endl;

    //Test the region C
    vec3 p_c1(-1.,  2.,  0.);
    vec3 p_c2( 0.,  1.,  0.);

    const auto [p_c1_dist_sq, p_c1_cpa_bary] = contact_model_abstract::compute_node_triangle_distance(p_c1, A, B, C);
    const auto [p_c2_dist_sq, p_c2_cpa_bary] = contact_model_abstract::compute_node_triangle_distance(p_c2, A, B, C);

    bool t9  = p_c1_dist_sq == 2.;
    bool t10 = p_c1_cpa_bary == vec3(0., 0. , 1.);

    bool t11 = p_c2_dist_sq == 0.;
    bool t12 = p_c2_cpa_bary == vec3(0., 0. , 1.);

    std::cout << "t9  " << t9  << std::endl;
    std::cout << "t10 " << t10 << std::endl;
    std::cout << "t11 " << t11 << std::endl;
    std::cout << "t12 " << t12 << std::endl;

    //Test the region AB
    vec3 p_ab1( 0.5, -1.,  0.);
    vec3 p_ab2( 0.5,  0.,  0.);

    const auto [p_ab1_dist_sq, p_ab1_cpa_bary] = contact_model_abstract::compute_node_triangle_distance(p_ab1, A, B, C);
    const auto [p_ab2_dist_sq, p_ab2_cpa_bary] = contact_model_abstract::compute_node_triangle_distance(p_ab2, A, B, C);


    bool t13 = p_ab1_dist_sq == 1.;
    bool t14 = p_ab1_cpa_bary == vec3(0.5, 0.5 , 0.);

    bool t15 = p_ab2_dist_sq == 0.;
    bool t16 = p_ab2_cpa_bary == vec3(0.5, 0.5 , 0.);

    std::cout << "t13 " << t13 << std::endl;
    std::cout << "t14 " << t14 << std::endl;
    std::cout << "t15 " << t15 << std::endl;
    std::cout << "t16 " << t16 << std::endl;

    //Test the region AC
    vec3 p_ac1(-1.,  0.5,  0.);
    vec3 p_ac2( 0.,  0.5,  0.);

    const auto [p_ac1_dist_sq, p_ac1_cpa_bary] = contact_model_abstract::compute_node_triangle_distance(p_ac1, A, B, C);
    const auto [p_ac2_dist_sq, p_ac2_cpa_bary] = contact_model_abstract::compute_node_triangle_distance(p_ac2, A, B, C);

    bool t17 = p_ac1_dist_sq == 1.;
    bool t18 = p_ac1_cpa_bary == vec3(0.5, 0. , 0.5);

    bool t19 = p_ac2_dist_sq == 0.;
    bool t20 = p_ac2_cpa_bary == vec3(0.5, 0. , 0.5);

    std::cout << "t17 " << t17 << std::endl;
    std::cout << "t18 " << t18 << std::endl;
    std::cout << "t19 " << t19 << std::endl;
    std::cout << "t20 " << t20 << std::endl;

    //Test the region BC
    vec3 p_bc1( 1.,    1.,  0.);
    vec3 p_bc2( 0.5,  0.5,  0.);

    const auto [p_bc1_dist_sq, p_bc1_cpa_bary] = contact_model_abstract::compute_node_triangle_distance(p_bc1, A, B, C);
    const auto [p_bc2_dist_sq, p_bc2_cpa_bary] = contact_model_abstract::compute_node_triangle_distance(p_bc2, A, B, C);

    bool t21 = p_bc1_dist_sq == 0.5;
    bool t22 = p_bc1_cpa_bary == vec3(0., 0.5 , 0.5);

    bool t23 = p_bc2_dist_sq == 0.;
    bool t24 = p_bc2_cpa_bary == vec3(0., 0.5 , 0.5);

    std::cout << "t21 " << t21 << std::endl;
    std::cout << "t22 " << t22 << std::endl;
    std::cout << "t23 " << t23 << std::endl;
    std::cout << "t24 " << t24 << std::endl;

    //Finally test the region ABC
    vec3 p_abc1( 0.25,  0.25,   1.);
    vec3 p_abc2( 0.25,  0.25,   0.);

    const auto [p_abc1_dist_sq, p_abc1_cpa_bary] = contact_model_abstract::compute_node_triangle_distance(p_abc1, A, B, C);
    const auto [p_abc2_dist_sq, p_abc2_cpa_bary] = contact_model_abstract::compute_node_triangle_distance(p_abc2, A, B, C);

    bool t25 = p_abc1_dist_sq == 1.;
    bool t26 = p_abc1_cpa_bary == vec3(0.5, 0.25 , 0.25);

    bool t27 = p_abc2_dist_sq == 0.;
    bool t28 = p_abc2_cpa_bary == vec3(0.5, 0.25 , 0.25);

    std::cout << "t25 " << t25 << std::endl;
    std::cout << "t26 " << t26 << std::endl;

    std::cout << "t27 " << t27 << std::endl;
    std::cout << "t28 " << t28 << std::endl;


    return !(t1 && t2 && t3 && t4 && t5 && t6 && t7 && t8 && t9 && t10 && t11 && t12 && t13 && t14 && t15 && t16 && t17 && t18 && t19 && t20 && t21 && t22 && t23 && t24 && t25 && t26 && t27 && t28);
}




//---------------------------------------------------------------------------------------------------------
// The main function
int main (int argc, char** argv){

    //Check that the command line input is correctly formatted
    assert(argc == 2); 

    //Get the name of the test to run
    std::string test_name = argv[1];
    
    //Run the selected test
    tester_contact_model_abstract tester;

    if (test_name == "compute_node_triangle_distance_test")               return tester.compute_node_triangle_distance_test();

    std::cout << "TEST NAME :" << test_name << " DOES NOT EXIST" << std::endl;
    return 1;

}