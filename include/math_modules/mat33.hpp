#ifndef DEF_MAT33
#define DEF_MAT33


#include <iostream>
#include <limits>
#include <iomanip>
#include <cassert>
#include <cmath>
#include <array>

#include "utils.hpp"
#include "vec3.hpp"
#include "eigen_solver.hpp"



/*
    Class for a 3x3 matrix object
*/

class mat33 
{

    private: 

        //Initialize the 3 rows of the matrix
        std::array<double, 3> row_1_ = {0., 0., 0.};
        std::array<double, 3> row_2_ = {0., 0., 0.};
        std::array<double, 3> row_3_ = {0., 0., 0.};

    public:

        //Matrix is a zero matrix with the default consyucor 
        mat33() = default;                        //Default constructor
        mat33(const mat33& m) = default;           //copy constructor
        mat33(mat33&& m) = default;                //move constructor
        mat33& operator=(const mat33& m) = default;//copy assignment operator
        mat33& operator=(mat33&& m) = default;     //move assignment operator 


        //Trivial constructor
        mat33(const std::array<double, 3>& row_1,
              const std::array<double, 3>& row_2,
              const std::array<double, 3>& row_3
            )noexcept;


        mat33(
            std::array<double, 3>&& row_1,
            std::array<double, 3>&& row_2,
            std::array<double, 3>&& row_3
        )noexcept;


        static mat33 identity() noexcept{
            return mat33(
                {1., 0., 0.},
                {0., 1., 0.},
                {0., 0., 1.}
            );
        }

        //Print the content of the vector
        void print() const noexcept; 


        vec3 get_row(const int i) const noexcept;
        vec3 get_col(const int i) const noexcept;

        //Return the dot product
        mat33 dot(const mat33& m) const noexcept;
        mat33 dot(      mat33&& m) const noexcept;

        vec3 dot(const vec3& v) const noexcept;
        vec3 dot(      vec3&& v) const noexcept;
       

        mat33 transpose() const noexcept;
        double determinant() const noexcept;

        std::pair<vec3, mat33> eigen_decomposition() const noexcept;


        mat33 inverse() const noexcept;


        //Overload substraction operator
        mat33 operator-(mat33&& m) const noexcept;         //rvalue arg
        mat33 operator-(const mat33& m) const noexcept;    //lvalue arg

        //Overload addition operator
        mat33 operator+(mat33&& m) const noexcept;         //rvalue arg
        mat33 operator+(const mat33& m) const noexcept;    //lvalue arg

        mat33 operator*(const double scalar) const noexcept;         //rvalue arg
        std::array<double, 3> operator[](const int i) const noexcept;    //lvalue arg


};

#endif