#include <cassert>
#include <string>
#include "utils.hpp"

#include "mat33.hpp"




/*
Contains all the tests run on the vec3 class which is contained in the vec3.hpp file
*/


//---------------------------------------------------------------------------------------------------------
//Test the trivial constructor
int test_constructors(){
    mat33 m1(
        {1., 2., 3.},
        {4., 5., 6.},
        {7., 8., 9.}
    );


    bool t1 = m1[0][0] == 1 && m1[0][1] == 2 && m1[0][2] == 3;
    bool t2 = m1[1][0] == 4 && m1[1][1] == 5 && m1[1][2] == 6;
    bool t3 = m1[2][0] == 7 && m1[2][1] == 8 && m1[2][2] == 9;

    return !(t1 && t2 && t3);
}
//---------------------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------------------
int get_row_col_test(){

 
    //Test the default constructor
    mat33 m1(
        {1., 2., 3.},
        {4., 5., 6.},
        {7., 8., 9.}
    );


    bool t1 = m1.get_row(0) == vec3(1., 2., 3.);
    bool t2 = m1.get_row(1) == vec3(4., 5., 6.);
    bool t3 = m1.get_row(2) == vec3(7., 8., 9.);

    bool t4 = m1.get_col(0) == vec3(1., 4., 7.);
    bool t5 = m1.get_col(1) == vec3(2., 5., 8.);
    bool t6 = m1.get_col(2) == vec3(3., 6., 9.);

    return !(t1 && t2 && t3 && t4 && t5 && t6);
}
//---------------------------------------------------------------------------------------------------------


//---------------------------------------------------------------------------------------------------------
int matrix_dot_test(){

    //Test the default constructor
    mat33 m1(
        {1., 2., 3.},
        {4., 5., 6.},
        {7., 8., 9.}
    );

    mat33 m2(
        {1., 2., 3.},
        {4., 5., 6.},
        {7., 8., 9.}
    );

    mat33 m3 = m1.dot(m2);

    bool t1 = m3.get_col(0) == m1.dot(m2.get_col(0));
    bool t2 = m3.get_col(1) == m1.dot(m2.get_col(1));
    bool t3 = m3.get_col(2) == m1.dot(m2.get_col(2));

    std::cout << t1 << std::endl;
    std::cout << t2 << std::endl;
    std::cout << t3 << std::endl;


    return !(t1 && t2 && t3);
}
//---------------------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------------------
int vector_dot_test(){
    
    mat33 m1(        
        {1., 2., 3.},
        {4., 5., 6.},
        {7., 8., 9.}
    );

    vec3 v1(1., 2., 3.);
    vec3 v2 = m1.dot(v1);


    bool t1 = v2 == m1.get_col(0) * v1.dx() + m1.get_col(1) * v1.dy() + m1.get_col(2) * v1.dz();
    return !t1;
}
//---------------------------------------------------------------------------------------------------------


//---------------------------------------------------------------------------------------------------------
int transpose_test(){
    mat33 m1(        
        {1., 2., 3.},
        {4., 5., 6.},
        {7., 8., 9.}
    );

    mat33 m2 = m1.transpose();


    bool t1 = m2.get_col(0) == m1.get_row(0);
    bool t2 = m2.get_col(1) == m1.get_row(1);
    bool t3 = m2.get_col(2) == m1.get_row(2);

    return !(t1 && t2 && t3);
}
//---------------------------------------------------------------------------------------------------------


//---------------------------------------------------------------------------------------------------------
int determinant_test(){

    mat33 m1(        
        {6., 1. , 1.},
        {4., -2., 5.},
        {2., 8. , 7.}
    );

    double det = m1.determinant();
    return !(det == -306.);
}
//---------------------------------------------------------------------------------------------------------




//---------------------------------------------------------------------------------------------------------
int inverse_test(){

    mat33 m1(        
        {6., 1. , 1.},
        {4., -2., 5.},
        {2., 8. , 7.}
    );

    mat33 m2 = m1.inverse();

    m2.print();



    bool t1 = m2.get_row(0) == vec3( 3./17.,	 -1./306.,	-7./306.);
    bool t2 = m2.get_row(1) == vec3( 1./17.,	-20./153.,    13./153.);
    bool t3 = m2.get_row(2) == vec3(-2./17.,	 23./153.,	 8./153.);

    std::cout << t1 << std::endl;
    std::cout << t2 << std::endl;
    std::cout << t3 << std::endl;

    return !(t1 && t2 && t3);
}
//---------------------------------------------------------------------------------------------------------


//---------------------------------------------------------------------------------------------------------

int operator_minus_test(){
    mat33 m1(
        {1., 2., 3.},
        {4., 5., 6.},
        {7., 8., 9.}
    );
    mat33 m2 = m1;
    mat33 m3 = m1 - m2;


    bool t1 = m3.get_row(0) == vec3(0., 0., 0.);
    bool t2 = m3.get_row(1) == vec3(0., 0., 0.);
    bool t3 = m3.get_row(2) == vec3(0., 0., 0.);

    return !(t1 && t2 && t3);
}

//---------------------------------------------------------------------------------------------------------




//---------------------------------------------------------------------------------------------------------

int operator_plus_test(){
    mat33 m1(
        {1., 2., 3.},
        {4., 5., 6.},
        {7., 8., 9.}
    );
    mat33 m2 = m1;
    mat33 m3 = m1 + m2;


    bool t1 = m3.get_row(0) == vec3(2., 4., 6.);
    bool t2 = m3.get_row(1) == vec3(8., 10., 12.);
    bool t3 = m3.get_row(2) == vec3(14., 16., 18.);

    return !(t1 && t2 && t3);
}

//---------------------------------------------------------------------------------------------------------


//---------------------------------------------------------------------------------------------------------
int operator_multiply_test(){
    mat33 m1(
        {1., 2., 3.},
        {4., 5., 6.},
        {7., 8., 9.}
    );
    mat33 m2 = m1 * 2.;


    bool t1 = m2.get_row(0) == vec3(2., 4., 6.);
    bool t2 = m2.get_row(1) == vec3(8., 10., 12.);
    bool t3 = m2.get_row(2) == vec3(14., 16., 18.);

    return !(t1 && t2 && t3);
}
//---------------------------------------------------------------------------------------------------------


//---------------------------------------------------------------------------------------------------------
int identity_test(){
    mat33 m1 = mat33::identity();

    bool t1 = m1.get_row(0) == vec3(1., 0., 0.);
    bool t2 = m1.get_row(1) == vec3(0., 1., 0.);
    bool t3 = m1.get_row(2) == vec3(0., 0., 1.);
    return !(t1 && t2 && t3);
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
    if (test_name == "get_row_col_test")                return get_row_col_test();
    if (test_name == "matrix_dot_test")                 return matrix_dot_test();
    if (test_name == "vector_dot_test")                 return vector_dot_test();
    if (test_name == "transpose_test")                  return transpose_test();
    if (test_name == "determinant_test")                return determinant_test();
    if (test_name == "inverse_test")                    return inverse_test();
    if (test_name == "operator_minus_test")             return operator_minus_test();
    if (test_name == "operator_plus_test")              return operator_plus_test();
    if (test_name == "operator_multiply_test")          return operator_multiply_test();
    if (test_name == "identity_test")                   return identity_test();







    

    


    




    
    std::cout << "TEST NAME :" << test_name << " DOES NOT EXIST" << std::endl;
    return 1;
}
//---------------------------------------------------------------------------------------------------------



