#include "mat33.hpp"

//----------------------------------------------------------------------
mat33::mat33(const std::array<double, 3>& row_1,
             const std::array<double, 3>& row_2,
             const std::array<double, 3>& row_3
            ) noexcept
{
    row_1_ = row_1;
    row_2_ = row_2;
    row_3_ = row_3;
}

mat33::mat33(std::array<double, 3>&& row_1,
             std::array<double, 3>&& row_2,
             std::array<double, 3>&& row_3
            ) noexcept
{
    row_1_ = std::move(row_1);
    row_2_ = std::move(row_2);
    row_3_ = std::move(row_3);
}
//----------------------------------------------------------------------


//----------------------------------------------------------------------


//----------------------------------------------------------------------
void mat33::print() const noexcept{
    printf("|%.1e, %.1e, %.1e|\n", row_1_[0], row_1_[1], row_1_[2]);
    printf("|%.1e, %.1e, %.1e|\n", row_2_[0], row_2_[1], row_2_[2]);
    printf("|%.1e, %.1e, %.1e|\n", row_3_[0], row_3_[1], row_3_[2]);
}
//----------------------------------------------------------------------

//----------------------------------------------------------------------
vec3 mat33::get_row(const int i) const noexcept{
    assert(i >= 0 && i < 3);
    if(i == 0) return       vec3(row_1_[0], row_1_[1], row_1_[2]);
    else if(i == 1) return  vec3(row_2_[0], row_2_[1], row_2_[2]);
    else return             vec3(row_3_[0], row_3_[1], row_3_[2]);
}
vec3 mat33::get_col(const int i) const noexcept{
    assert(i >= 0 && i < 3);
    if(i == 0) return       vec3(row_1_[0], row_2_[0], row_3_[0]);
    else if(i == 1) return  vec3(row_1_[1], row_2_[1], row_3_[1]);
    else return             vec3(row_1_[2], row_2_[2], row_3_[2]);

}
//----------------------------------------------------------------------



//----------------------------------------------------------------------
std::array<double, 3> mat33::operator[](const int i) const noexcept{
    assert(i >= 0 && i < 3);
    if(i == 0) return row_1_;
    else if(i == 1) return row_2_;
    else return row_3_;
}
//----------------------------------------------------------------------



//----------------------------------------------------------------------
mat33 mat33::dot(const mat33& m) const noexcept
{
    mat33 result;
    result.row_1_[0] = row_1_[0]*m.row_1_[0] + row_1_[1]*m.row_2_[0] + row_1_[2]*m.row_3_[0];
    result.row_1_[1] = row_1_[0]*m.row_1_[1] + row_1_[1]*m.row_2_[1] + row_1_[2]*m.row_3_[1];
    result.row_1_[2] = row_1_[0]*m.row_1_[2] + row_1_[1]*m.row_2_[2] + row_1_[2]*m.row_3_[2];

    result.row_2_[0] = row_2_[0]*m.row_1_[0] + row_2_[1]*m.row_2_[0] + row_2_[2]*m.row_3_[0];
    result.row_2_[1] = row_2_[0]*m.row_1_[1] + row_2_[1]*m.row_2_[1] + row_2_[2]*m.row_3_[1];
    result.row_2_[2] = row_2_[0]*m.row_1_[2] + row_2_[1]*m.row_2_[2] + row_2_[2]*m.row_3_[2];

    result.row_3_[0] = row_3_[0]*m.row_1_[0] + row_3_[1]*m.row_2_[0] + row_3_[2]*m.row_3_[0];
    result.row_3_[1] = row_3_[0]*m.row_1_[1] + row_3_[1]*m.row_2_[1] + row_3_[2]*m.row_3_[1];
    result.row_3_[2] = row_3_[0]*m.row_1_[2] + row_3_[1]*m.row_2_[2] + row_3_[2]*m.row_3_[2];

    return result;
}


mat33 mat33::dot(mat33&& m) const noexcept
{
    mat33 result;
    result.row_1_[0] = row_1_[0]*m.row_1_[0] + row_1_[1]*m.row_2_[0] + row_1_[2]*m.row_3_[0];
    result.row_1_[1] = row_1_[0]*m.row_1_[1] + row_1_[1]*m.row_2_[1] + row_1_[2]*m.row_3_[1];
    result.row_1_[2] = row_1_[0]*m.row_1_[2] + row_1_[1]*m.row_2_[2] + row_1_[2]*m.row_3_[2];

    result.row_2_[0] = row_2_[0]*m.row_1_[0] + row_2_[1]*m.row_2_[0] + row_2_[2]*m.row_3_[0];
    result.row_2_[1] = row_2_[0]*m.row_1_[1] + row_2_[1]*m.row_2_[1] + row_2_[2]*m.row_3_[1];
    result.row_2_[2] = row_2_[0]*m.row_1_[2] + row_2_[1]*m.row_2_[2] + row_2_[2]*m.row_3_[2];

    result.row_3_[0] = row_3_[0]*m.row_1_[0] + row_3_[1]*m.row_2_[0] + row_3_[2]*m.row_3_[0];
    result.row_3_[1] = row_3_[0]*m.row_1_[1] + row_3_[1]*m.row_2_[1] + row_3_[2]*m.row_3_[1];
    result.row_3_[2] = row_3_[0]*m.row_1_[2] + row_3_[1]*m.row_2_[2] + row_3_[2]*m.row_3_[2];

    return result;
}
//----------------------------------------------------------------------


//----------------------------------------------------------------------
vec3 mat33::dot(const vec3& v) const noexcept{
    double dx = row_1_[0]*v.dx() + row_1_[1]*v.dy() + row_1_[2]*v.dz();
    double dy = row_2_[0]*v.dx() + row_2_[1]*v.dy() + row_2_[2]*v.dz();
    double dz = row_3_[0]*v.dx() + row_3_[1]*v.dy() + row_3_[2]*v.dz();
    return vec3(dx, dy, dz);
}

vec3 mat33::dot(vec3&& v) const noexcept{
    double dx = row_1_[0]*v.dx() + row_1_[1]*v.dy() + row_1_[2]*v.dz();
    double dy = row_2_[0]*v.dx() + row_2_[1]*v.dy() + row_2_[2]*v.dz();
    double dz = row_3_[0]*v.dx() + row_3_[1]*v.dy() + row_3_[2]*v.dz();
    return vec3(dx, dy, dz);
}
//----------------------------------------------------------------------



//----------------------------------------------------------------------
mat33 mat33::operator-(mat33&& m) const noexcept{
    mat33 result;
    result.row_1_[0] = row_1_[0] - m.row_1_[0];
    result.row_1_[1] = row_1_[1] - m.row_1_[1];
    result.row_1_[2] = row_1_[2] - m.row_1_[2];

    result.row_2_[0] = row_2_[0] - m.row_2_[0];
    result.row_2_[1] = row_2_[1] - m.row_2_[1];
    result.row_2_[2] = row_2_[2] - m.row_2_[2];

    result.row_3_[0] = row_3_[0] - m.row_3_[0];
    result.row_3_[1] = row_3_[1] - m.row_3_[1];
    result.row_3_[2] = row_3_[2] - m.row_3_[2];

    return result;
}


mat33 mat33::operator-(const mat33& m) const noexcept{
    mat33 result;
    result.row_1_[0] = row_1_[0] - m.row_1_[0];
    result.row_1_[1] = row_1_[1] - m.row_1_[1];
    result.row_1_[2] = row_1_[2] - m.row_1_[2];

    result.row_2_[0] = row_2_[0] - m.row_2_[0];
    result.row_2_[1] = row_2_[1] - m.row_2_[1];
    result.row_2_[2] = row_2_[2] - m.row_2_[2];

    result.row_3_[0] = row_3_[0] - m.row_3_[0];
    result.row_3_[1] = row_3_[1] - m.row_3_[1];
    result.row_3_[2] = row_3_[2] - m.row_3_[2];

    return result;
}




//----------------------------------------------------------------------

mat33 mat33::operator+(mat33&& m) const noexcept{
    mat33 result;
    result.row_1_[0] = row_1_[0] + m.row_1_[0];
    result.row_1_[1] = row_1_[1] + m.row_1_[1];
    result.row_1_[2] = row_1_[2] + m.row_1_[2];

    result.row_2_[0] = row_2_[0] + m.row_2_[0];
    result.row_2_[1] = row_2_[1] + m.row_2_[1];
    result.row_2_[2] = row_2_[2] + m.row_2_[2];

    result.row_3_[0] = row_3_[0] + m.row_3_[0];
    result.row_3_[1] = row_3_[1] + m.row_3_[1];
    result.row_3_[2] = row_3_[2] + m.row_3_[2];

    return result;
}


mat33 mat33::operator+(const mat33& m) const noexcept{
    mat33 result;
    result.row_1_[0] = row_1_[0] + m.row_1_[0];
    result.row_1_[1] = row_1_[1] + m.row_1_[1];
    result.row_1_[2] = row_1_[2] + m.row_1_[2];

    result.row_2_[0] = row_2_[0] + m.row_2_[0];
    result.row_2_[1] = row_2_[1] + m.row_2_[1];
    result.row_2_[2] = row_2_[2] + m.row_2_[2];

    result.row_3_[0] = row_3_[0] + m.row_3_[0];
    result.row_3_[1] = row_3_[1] + m.row_3_[1];
    result.row_3_[2] = row_3_[2] + m.row_3_[2];

    return result;
}

//----------------------------------------------------------------------


//----------------------------------------------------------------------
mat33 mat33::transpose() const noexcept{
    mat33 result;
    result.row_1_[0] = row_1_[0];
    result.row_1_[1] = row_2_[0];
    result.row_1_[2] = row_3_[0];

    result.row_2_[0] = row_1_[1];
    result.row_2_[1] = row_2_[1];
    result.row_2_[2] = row_3_[1];

    result.row_3_[0] = row_1_[2];
    result.row_3_[1] = row_2_[2];
    result.row_3_[2] = row_3_[2];

    return result;
}
//----------------------------------------------------------------------


//----------------------------------------------------------------------
double mat33::determinant() const noexcept{
    return row_1_[0]*(row_2_[1]*row_3_[2] - row_2_[2]*row_3_[1]) -
           row_1_[1]*(row_2_[0]*row_3_[2] - row_2_[2]*row_3_[0]) +
           row_1_[2]*(row_2_[0]*row_3_[1] - row_2_[1]*row_3_[0]);
}
//----------------------------------------------------------------------


//----------------------------------------------------------------------
//Compute the eigenvalues and eigenvectors of a symmetric 3x3 matrix
//Return the sorted eigen values in a vector and the corresponding eigenvectors in a matrix
std::pair<vec3, mat33> mat33::eigen_decomposition() const noexcept{
    assert(row_1_[1] == row_2_[0]);
    assert(row_1_[2] == row_3_[0]);
    assert(row_2_[2] == row_3_[1]);

    //Use a dedicated solver for 3x3 symmetric matrices
    gte::SymmetricEigensolver3x3<double> eigen_solver;
    std::array<double, 3> eigen_values;
    std::array<std::array<double, 3>, 3> eigen_vectors;

    eigen_solver(row_1_[0], row_1_[1], row_1_[2],
                           row_2_[1],  row_2_[2],
                                       row_3_[2],
                true, 1, eigen_values, eigen_vectors);

    //Each column of the matrix is an eigenvector
    mat33 eigen_vector_matrix(  {eigen_vectors[0][0], eigen_vectors[0][1], eigen_vectors[0][2]},
                                {eigen_vectors[1][0], eigen_vectors[1][1], eigen_vectors[1][2]},
                                {eigen_vectors[2][0], eigen_vectors[2][1], eigen_vectors[2][2]}); 
                                
    eigen_vector_matrix = eigen_vector_matrix.transpose();              
    vec3 eigen_values_vec = vec3(eigen_values[0], eigen_values[1], eigen_values[2]);

    return std::make_pair(eigen_values_vec, eigen_vector_matrix);
}

//----------------------------------------------------------------------


//----------------------------------------------------------------------
mat33 mat33::inverse() const noexcept{
    mat33 result;
    double a = row_1_[0];
    double b = row_1_[1];
    double c = row_1_[2];
    double d = row_2_[0];
    double e = row_2_[1];
    double f = row_2_[2];
    double g = row_3_[0];
    double h = row_3_[1];
    double i = row_3_[2];

    double det = determinant();
    result.row_1_[0] = (e*i - f*h)/det;
    result.row_1_[1] = (c*h - b*i)/det;
    result.row_1_[2] = (b*f - c*e)/det;

    result.row_2_[0] = (f*g - d*i)/det;
    result.row_2_[1] = (a*i - c*g)/det;
    result.row_2_[2] = (c*d - a*f)/det;

    result.row_3_[0] = (d*h - e*g)/det;
    result.row_3_[1] = (b*g - a*h)/det;
    result.row_3_[2] = (a*e - b*d)/det;

    return result;


}
//----------------------------------------------------------------------


//----------------------------------------------------------------------
mat33 mat33::operator*(const double scalar) const noexcept{
    mat33 result;
    result.row_1_[0] = row_1_[0] * scalar;
    result.row_1_[1] = row_1_[1] * scalar;
    result.row_1_[2] = row_1_[2] * scalar;

    result.row_2_[0] = row_2_[0] * scalar;
    result.row_2_[1] = row_2_[1] * scalar;
    result.row_2_[2] = row_2_[2] * scalar;

    result.row_3_[0] = row_3_[0] * scalar;
    result.row_3_[1] = row_3_[1] * scalar;
    result.row_3_[2] = row_3_[2] * scalar;

    return result;
}
//----------------------------------------------------------------------
