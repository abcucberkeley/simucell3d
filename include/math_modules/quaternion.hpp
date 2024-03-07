#ifndef DEF_QUATERNION
#define DEF_QUATERNION

#include <iostream>
#include <cmath>

#include "mat33.hpp"


//-----------------------------------------------------------------------------------------------------
class quaternion{
    //Quaternions are used during the cell divisions to normalize the rotation matrix

    private:
        double w=0., i=0., j=0., k=0.;

    public:
        quaternion() = default;                        //Default constructor
        quaternion(const quaternion& q) = default;           //copy constructor
        quaternion(quaternion&& q) = default;                //move constructor
        quaternion& operator=(const quaternion& q) = default;//copy assignment operator
        quaternion& operator=(quaternion&& q) = default;     //move assignment operator 

        //Construct the quaternion from 4 doubles
        quaternion(double nw, double ni, double nj, double nk): w(nw), i(ni), j(nj) , k(nk){}

        //Construct the quaternion from a matrix
        quaternion(const mat33& M1){
            // Get the real part of the quaternion first
            w = std::sqrt( 1.0 + M1[0][0]+M1[1][1]+ M1[2][2]) *0.5;
            i = (M1[2][1]-M1[1][2])/(4.0*w);
            j = (M1[0][2]-M1[2][0])/(4.0*w);
            k = (M1[1][0]-M1[0][1])/(4.0*w);

            assert(std::isfinite(w));
        }

        quaternion operator*(const double scalar) const{
            double nw = w * scalar;
            double ni = i * scalar;
            double nj = j * scalar;
            double nk = k * scalar;
            return quaternion(nw, ni, nj, nk);
        }

        quaternion operator/(double scalar) const{
            double nw = w / scalar;
            double ni = i / scalar;
            double nj = j / scalar;
            double nk = k / scalar;
            return quaternion(nw, ni, nj, nk);
        }

      

        quaternion normalize() const{
            double norm = std::sqrt(w*w + i*i + j*j + k*k);
            return quaternion(w / norm, i / norm, j / norm, k / norm);
        }

        quaternion inverse() const{
            return quaternion(w, -i, -j, -k);
        }

        void print() const{
            printf("[%.2e, %.2ei, %.2ej, %.2ek]\n", w, i, j, k);
        }

        static quaternion from_matrix(const mat33& M1){
            // Get the real part of the quaternion first
            double r = std::sqrt( 1.0 + M1[0][0]+M1[1][1]+ M1[2][2]) *0.5;
            assert(!std::isnan(r));
            double i = (M1[2][1]-M1[1][2])/(4.0*r);
            double j = (M1[0][2]-M1[2][0])/(4.0*r);
            double k = (M1[1][0]-M1[0][1])/(4.0*r);
            return quaternion(r,i,j,k);
        }


        mat33 to_matrix() const{
            double qw = w;
            double qx = i; 
            double qy = j; 
            double qz = k; 

            double I11 = 1.0- 2.0*qz*qz- 2.0*qy*qy;
            double I12 = -2.0* qz*qw+2.0*qy*qx;
            double I13 = 2.0*qy*qw +2.0* qz*qx;
            double I21 = 2.0*qx*qy+ 2.0*qw*qz;
            double I22 = 1.0 - 2.0*qz*qz - 2.0*qx*qx;
            double I23 = 2.0*qz*qy- 2.0*qx*qw;
            double I31 = 2.0*qx*qz- 2.0*qw*qy;
            double I32 = 2.0*qy*qz + 2.0*qw*qx;
            double I33 = 1.0- 2.0*qy*qy- 2.0*qx*qx;

            //Calculate the skew rotational displacement mat33
            mat33 matrix(
                {I11, I12, I13}, 
                {I21, I22, I23}, 
                {I31, I32, I33}
            );
            return matrix;
        }
};

#endif