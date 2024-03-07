#include "vec3.hpp"



vec3::vec3(const double dx,  const double dy, const double dz): 
    dx_(dx), dy_(dy), dz_(dz)
    
{
    //Check the components of the vector
    assert(std::isfinite(dx) && std::isfinite(dy) && std::isfinite(dz));

}


//--------------------------------------------------------------------------------------------------------------------------
//Print the vec3 coordinates in scientific notation
void vec3::print() const{
    std::cout << std::setprecision(2) << std::scientific  << "(" << dx_ << ", " << dy_ << ", "<< dz_ << ")" <<std::endl; 
}
//--------------------------------------------------------------------------------------------------------------------------


//--------------------------------------------------------------------------------------------------------------------------
//Return the dot product between this vector and the one given in argument
double vec3::dot(const vec3& v) const{
    return dx_* v.dx() + dy_* v.dy() + dz_* v.dz(); 
}
//--------------------------------------------------------------------------------------------------------------------------

//--------------------------------------------------------------------------------------------------------------------------

//The translate function acts on the current vector, it avoids allocating new memory slots
void vec3::translate(const vec3& v){
    assert(std::isfinite(v.dx()) && std::isfinite(v.dy()) && std::isfinite(v.dz()));

    #pragma omp atomic update
    dx_ += v.dx();

    #pragma omp atomic update
    dy_ += v.dy();

    #pragma omp atomic update
    dz_ += v.dz();
}



//Translate the values by a certain amount
void vec3::translate(const double dx, const double dy, const double dz){
    assert(std::isfinite(dx) && std::isfinite(dy) && std::isfinite(dz));

    #pragma omp atomic update
    dx_ += dx;

    #pragma omp atomic update
    dy_ += dy;

    #pragma omp atomic update
    dz_ += dz;
}
//--------------------------------------------------------------------------------------------------------------------------


//--------------------------------------------------------------------------------------------------------------------------
//Compute the cross product of 2 vectors
vec3 vec3::cross(const vec3& v) const{
    return vec3(
        dy_*v.dz() - dz_*v.dy(),
        dz_*v.dx() - dx_*v.dz(),
        dx_*v.dy() - dy_*v.dx()
    );
}
//--------------------------------------------------------------------------------------------------------------------------


//--------------------------------------------------------------------------------------------------------------------------
//Sum 2 vectors
vec3 vec3::operator+(vec3&& v) const{
    return vec3(dx_ + v.dx(), dy_ + v.dy(), dz_ + v.dz());
}

//Sum two vec3s
vec3 vec3::operator+(const vec3&  v) const{
    return vec3(dx_ + v.dx(), dy_ + v.dy(), dz_ + v.dz());
}
//--------------------------------------------------------------------------------------------------------------------------

//--------------------------------------------------------------------------------------------------------------------------
//Return the unit vector of the given vector
vec3 vec3::normalize() const{
    const double norm = this->norm();

    //If the norm of the vector is 0, return a null vector
    return (norm != 0.) ? *(this) / norm : vec3(0., 0., 0.);
}
//--------------------------------------------------------------------------------------------------------------------------


//--------------------------------------------------------------------------------------------------------------------------
//Cap the vector to a maximum norm
bool vec3::cap(const double max_norm){
    const double norm = this->norm();
    if(norm > max_norm){
        dx_ *= max_norm / norm;
        dy_ *= max_norm / norm;
        dz_ *= max_norm / norm;
    }
    return norm > max_norm;
}
//--------------------------------------------------------------------------------------------------------------------------



//--------------------------------------------------------------------------------------------------------------------------
//Get the angle between 2 vectors, the angle is in the range [0, pi]
double vec3::get_angle_with(const vec3& v) const{

    //Get the norm of the 2 vectors
    const double norm1 = this->norm();
    const double norm2 = v.norm();
    double angle = std::acos(dot(v) / (norm1 * norm2));

    //Prevent erros due to numerical precision
    if(std::isnan(angle)) angle = 1.;

    return angle;
}

double vec3::get_angle_with(vec3&& v) const{

    //Get the norm of the 2 vectors
    const double norm1 = this->norm();
    const double norm2 = v.norm();
    double angle = std::acos(dot(v) / (norm1 * norm2));

    //Prevent erros due to numerical precision
    if(!std::isfinite(angle)) angle = 1.;

    return angle;
}
//--------------------------------------------------------------------------------------------------------------------------


//--------------------------------------------------------------------------------------------------------------------------
//Rotate this vector around the given axis by the given angle
vec3 vec3::rotate_around_axis(const vec3& axis, const double angle) const{

    //Use Rodrigus formula to perform the rotation. 
    //This formula creates a non negligible numerical approximation error
    return (*this)  * std::cos(angle) +  axis.cross(*this) * std::sin(angle) + axis * (1. - std::cos(angle)) * (axis.dot(*this));

}






//--------------------------------------------------------------------------------------------------------------------------




//--------------------------------------------------------------------------------------------------------------------------
//Substract two vec3s
vec3 vec3::operator-(vec3 const &v) const{
    return vec3(dx_ - v.dx(), dy_ - v.dy(), dz_ - v.dz());
}

vec3 vec3::operator-(vec3&& v) const{
    return vec3(dx_ - v.dx(), dy_ - v.dy(), dz_ - v.dz());
}
//--------------------------------------------------------------------------------------------------------------------------


//--------------------------------------------------------------------------------------------------------------------------
//Multiply a vector with a scalar
vec3 vec3::operator*(const double scalar) const{
    return vec3(dx_ * scalar, dy_ * scalar, dz_ * scalar);
}
//--------------------------------------------------------------------------------------------------------------------------


//--------------------------------------------------------------------------------------------------------------------------
//Divide a vector with a scalar
vec3 vec3::operator/(const double scalar) const{
    return vec3(dx_ / scalar, dy_ / scalar, dz_ / scalar);
}
//--------------------------------------------------------------------------------------------------------------------------


//--------------------------------------------------------------------------------------------------------------------------
//Check if 2 vectors are equal
bool vec3::operator==(const vec3& v) const{
    if (!(almost_equal(dx_, v.dx()))) return false;
    if (!(almost_equal(dy_, v.dy()))) return false;
    if (!(almost_equal(dz_, v.dz()))) return false;
    return true;
}

//Check if 2 vectors are different
bool vec3::operator!=(const vec3& v) const{return !(*(this) == v);}
//--------------------------------------------------------------------------------------------------------------------------
