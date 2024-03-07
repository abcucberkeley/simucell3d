#include <cassert>
#include <string>
#include "utils.hpp"

#include "vec3.hpp"




/*
Contains all the tests run on the vec3 class which is contained in the vec3.hpp file
*/


//---------------------------------------------------------------------------------------------------------
//Test all the constructors of vec3
int test_constructors(){
    //Test the default constructor
    vec3 v1; 
    bool t1 = v1.dx() == 0. && v1.dy() == 0. && v1.dz() == 0.;

    //Test the trivail default constructor
    vec3 v2(1., 2., 3.); 
    bool t2 = v2.dx() == 1. && v2.dy() == 2. && v2.dz() == 3.;

    //Test the copy constructor
    vec3 v3(v2);
    bool t3 = v3.dx() == 1. && v3.dy() == 2. && v3.dz() == 3.;

    //Test the move constructor
    vec3 v4( vec3(1., 2., 3.) );
    bool t4 = v4.dx() == 1. && v4.dy() == 2. && v4.dz() == 3.;

    //Test the copy assignment operator
    vec3 v5 = v2;
    bool t5 = v5.dx() == 1. && v5.dy() == 2. && v5.dz() == 3.;

    //Test the move assignment operator
    vec3 v6 = vec3(1,2,3);
    bool t6 = v6.dx() == 1. && v6.dy() == 2. && v6.dz() == 3.;

    //Check that they are all equal to 0
    return !(t1 && t2 && t3 && t4 && t5 && t6);
}
//---------------------------------------------------------------------------------------------------------




//---------------------------------------------------------------------------------------------------------
//Test the coordinates methods
int coordinate_getter_test(){
    vec3 v(1., 2., 3.); 
    return !(v.dx() == 1. && v.dy() == 2. && v.dz() == 3.);
}

//---------------------------------------------------------------------------------------------------------



//---------------------------------------------------------------------------------------------------------
//Test addition with rvalue
int sum_rvalue_test(){

    //Create a default vec3 object
    vec3 v1(1., 2., 3.); 

    vec3 v2 = v1 + vec3(1., 2., 3.); 

    //Extract the coordinates of the vector
    const auto [dx, dy, dz] = v2.to_array();

    //Check that they are all equal to 0
    return !(dx == 2. && dy == 4. && dz == 6.);
}
//---------------------------------------------------------------------------------------------------------





//---------------------------------------------------------------------------------------------------------
//Test addition with lvalue
int sum_lvalue_test(){

    //Create a default vec3 object
    vec3 v0(1., 2., 3.); 
    vec3 v1(1., 2., 3.); 

    vec3 v2 = v0 + v1; 

    //Extract the coordinates of the vector
    const auto [dx, dy, dz] = v2.to_array();

    //Check that they are all equal to 0
    return !(dx == 2. && dy == 4. && dz == 6.);
}
//---------------------------------------------------------------------------------------------------------



//---------------------------------------------------------------------------------------------------------
int cap_test(){
    vec3 v1(2,0,0);
    v1.cap(1.);
    return !(v1 == vec3(1., 0., 0.));
}
//---------------------------------------------------------------------------------------------------------





//---------------------------------------------------------------------------------------------------------
//Test substraction with rvalue
int substraction_rvalue_test(){

    //Create a default vec3 object
    vec3 v1(1., 2., 3.); 

    vec3 v2 = v1 - vec3(1., 2., 3.); 

    //Extract the coordinates of the vector
    const auto [dx, dy, dz] = v2.to_array();

    //Check that they are all equal to 0
    return !(dx == 0. && dy == 0. && dz == 0.);
}
//---------------------------------------------------------------------------------------------------------




//---------------------------------------------------------------------------------------------------------
//Test substraction with lvalue
int substraction_lvalue_test(){

    //Create a default vec3 object
    vec3 v0(1., 2., 3.); 
    vec3 v1(1., 2., 3.); 

    vec3 v2 = v0 - v1; 

    //Extract the coordinates of the vector
    const auto [dx, dy, dz] = v2.to_array();

    //Check that they are all equal to 0
    return !(dx == 0. && dy == 0. && dz == 0.);
}
//---------------------------------------------------------------------------------------------------------





//---------------------------------------------------------------------------------------------------------
//Test multiplication of a vector with scalar
int multiplication_test(){

    //Create a default vec3 object
    vec3 v0(1., 2., 3.); 

    vec3 v1 = v0 * 2.; 

    //Extract the coordinates of the vector
    const auto [dx, dy, dz] = v1.to_array();

    //Check that they are all equal to 0
    return !(dx == 2. && dy == 4. && dz == 6.);
}
//---------------------------------------------------------------------------------------------------------


//---------------------------------------------------------------------------------------------------------
//Check the equality and non equale operators
int equality_test(){

    //Create a default vec3 object
    vec3 v0(1., 2., 3.); 
    vec3 v1(1., 2., 3.); 
    vec3 v2(2., 2., 3.); 

    bool t1 = v0 == v1;
    bool t2 = v0 != v2;
    bool t3 = (v0 == v2) == false;

    return !(t1 && t2 && t3);
}
//---------------------------------------------------------------------------------------------------------





//---------------------------------------------------------------------------------------------------------
//Test multiplication of a vector with scalar
int division_test(){

    //Create a default vec3 object
    vec3 v0(1., 2., 3.); 

    vec3 v1 = v0 / 2.; 

    //Extract the coordinates of the vector
    const auto [dx, dy, dz] = v1.to_array();

    //Check that they are all equal to 0
    return !(dx == 0.5 && dy == 1. && dz == 1.5);
}
//---------------------------------------------------------------------------------------------------------


//---------------------------------------------------------------------------------------------------------
//Test the norm method
int norm_test(){

    //Create a default vec3 object
    vec3 v0(1., 1., 1.); 

    //Allow a certain numeric tolerance
    return!(almost_equal(v0.norm(), std::sqrt(3)));
}
//---------------------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------------------
//Test the normalization function
int normalize_test(){

    //Create a default vec3 object
    vec3 v0(1., 2., 3.);
    vec3 v1 = v0.normalize(); 

    //Allow a certain numeric tolerance
    return!(almost_equal(v1.norm(), 1.));
}
//---------------------------------------------------------------------------------------------------------



//---------------------------------------------------------------------------------------------------------
//Test the dot product function
int dot_product_test(){

    //Create a default vec3 object
    vec3 v0(1., 2., 3.);
    vec3 v1(1., 2., 3.);

    //Allow a certain numeric tolerance
    return!(almost_equal(v0.dot(v1), 14.));
}
//---------------------------------------------------------------------------------------------------------




//---------------------------------------------------------------------------------------------------------
//Test the cross product function
int cross_product_test(){
    vec3 v1(1., 0., 0.);
    vec3 v2(0., 1., 0.);
    vec3 v3 = v1.cross(v2);

    return !(v3 == vec3(0. , 0., 1.));
}
//---------------------------------------------------------------------------------------------------------



//---------------------------------------------------------------------------------------------------------
//test that the values of the vector are all reset to 0
int reset_test(){
    vec3 v1(1., 2., 3.);
    vec3 v2(1., 2., 3.);

    //Reset v1 to (0., 0., 0.)
    v1.reset();

    //Reset v2 to a defined value
    v2.reset(2., 4., 6.);

    //Check the 2 reset operations
    bool t1 = (v1.dx() == 0. && v1.dy() == 0. && v1.dz() == 0.);
    bool t2 = (v2.dx() == 2. && v2.dy() == 4. && v2.dz() == 6.);

    return !(t1 && t2);
} 
//---------------------------------------------------------------------------------------------------------


//---------------------------------------------------------------------------------------------------------
//test that the values of the vector are all reset to 0
int translate_test(){
    vec3 v1(1., 2., 3.);
    v1.translate(1., 2., 3.);
    bool t1 = (v1.dx() == 2. && v1.dy() == 4. && v1.dz() == 6.);

    vec3 v2(1., 2., 3.);
    v2.translate(1., 2., 3.);
    bool t2 = (v2.dx() == 2. && v2.dy() == 4. && v2.dz() == 6.);

    return !(t1 && t2);
} 
//---------------------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------------------
int get_angle_with_test(){

    //Create 2 orthogonals vectors
    vec3 v1(1., 0., 0.);
    vec3 v2(0., 1., 0.);

    //Check that the angle between them is 90 degrees
    bool t1 = almost_equal(v1.get_angle_with(v2), M_PI/2.);

    //Create 2 colinear vectors
    vec3 v3(1., 0., 0.);

    //Check that the angle between them is 0 degrees
    bool t2 = almost_equal(v1.get_angle_with(v3), 0.);

    return !(t1 && t2);
}
//---------------------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------------------
int rotate_around_axis_test(){
    const double angle = M_PI/2.;
    vec3 v1(1., 0., 0.);
    vec3 v2(0., 1., 0.);
    vec3 v3 = v1.rotate_around_axis(v2, angle);

    v3.print();
    bool t1 = (v3 - vec3(0., 0., -1.)).norm() < 1e-10;

    std::cout << "t1 = " << t1 << std::endl;
    return !(t1);
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
    if (test_name == "test_constructors")               return test_constructors();
    if (test_name == "coordinate_getter_test")          return coordinate_getter_test(); 
    if (test_name == "sum_lvalue_test")                 return sum_lvalue_test();
    if (test_name == "sum_rvalue_test")                 return sum_rvalue_test();
    if (test_name == "substraction_lvalue_test")        return substraction_lvalue_test();
    if (test_name == "substraction_rvalue_test")        return substraction_rvalue_test();
    if (test_name == "multiplication_test")             return multiplication_test();
    if (test_name == "division_test")                   return division_test();
    if (test_name == "equality_test")                   return equality_test();
    if (test_name == "norm_test")                       return norm_test();
    if (test_name == "normalize_test")                  return normalize_test();
    if (test_name == "dot_product_test")                return dot_product_test();
    if (test_name == "cross_product_test")              return cross_product_test();
    if (test_name == "reset_test")                      return reset_test();
    if (test_name == "translate_test")                  return translate_test();
    if (test_name == "get_angle_with_test")             return get_angle_with_test();
    if (test_name == "rotate_around_axis_test")         return rotate_around_axis_test();
    if (test_name == "cap_test")                        return cap_test();



    



    std::cout << "TEST NAME :" << test_name << " DOES NOT EXIST" << std::endl;
    return 1;
}
//---------------------------------------------------------------------------------------------------------



