#include <cassert>
#include <string>
#include <iostream>
#include <vector>

#include "utils.hpp"
#include "custom_structures.hpp"
#include "initial_triangulation.hpp"


//-------------------------------------------------------------------------------------------
//Generate the poisson point cloud of a unit cube, the l_min is set to 0.1
std::vector<oriented_point> generate_poisson_point_cloud(const double l_min){

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


    //Triangulate the untriangulated mesh, acts in place on the mesh object
    initial_triangulation::coarse_triangulation(m);

    //Convert the mesh to a cell object
    cell_ptr c1 = initial_triangulation::convert_mesh_to_cell(m);

    //Generate the poisson point cloud of the surface, with a l_min of 0.1
    std::vector<oriented_point> poisson_point_cloud = initial_triangulation::generate_poisson_point_cloud(l_min, c1);

    return poisson_point_cloud;
}
//-------------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------------
// Test if point P lies inside the counterclockwise triangle f
int point_in_triangle(const mesh& m, const vec3& p, const std::vector<unsigned>& f){

    //Get the points making up the triangle
    vec3 a(m.node_pos_lst[f[0]*3], m.node_pos_lst[f[0]*3 + 1], m.node_pos_lst[f[0]*3 + 2]);
    vec3 b(m.node_pos_lst[f[1]*3], m.node_pos_lst[f[1]*3 + 1], m.node_pos_lst[f[1]*3 + 2]);
    vec3 c(m.node_pos_lst[f[2]*3], m.node_pos_lst[f[2]*3 + 1], m.node_pos_lst[f[2]*3 + 2]);

    // Translate point and triangle so that point lies at origin
    a = a - p;
    b = b - p;
    c = c - p;

    // Compute normal vectors for triangles pab and pbc
    vec3 u = b.cross(c);
    vec3 v = c.cross(a);

    // Make sure they are both pointing in the same direction
    if (u.dot(v) < 0.) return 0;

    // Compute normal vector for triangle pca
    vec3 w = a.cross(b);
    // Make sure it points in the same direction as the first two
    if (u.dot(w) < 0.) return 0;

    // Otherwise P must be in (or on) the triangle
    return 1;
}
//-------------------------------------------------------------------------------------------




//-------------------------------------------------------------------------------------------
//Test that all the points are separated by l_min
int point_min_distance_test(const double l_min = 0.1){

    //Generate a poisson point cloud
    std::vector<oriented_point> poisson_point_cloud = generate_poisson_point_cloud(l_min);

    //Place all the points in uspg_4D
    const double voxel_size = l_min / std::sqrt(3);
    uspg_4d<oriented_point> grid(0., 0., 0., 1., 1., 1., l_min, poisson_point_cloud.size());

    //Place the point in the space partitionning grid
    for(auto& p: poisson_point_cloud){grid.place_object(p, p.position_);}

    //Loop over the points and 
    for(auto& p: poisson_point_cloud){

        //Get the points nearby
        const auto nearby_point_lst = grid.get_neighborhood(p.position_);

        //Calculate the distance between the points
        for(const auto& nearby_point: nearby_point_lst){

            if(nearby_point.id_ != p.id_){
                const double distance = (nearby_point.position_ - p.position_).norm();
                if(distance < l_min) return 1;
            }
        }
    }
    return 0;
}
//-------------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------------
//Test that all the points are located on the cell surface
int point_on_surface_test(const double l_min = 0.1){

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


    //Triangulate the untriangulated mesh, acts in place on the mesh object
    initial_triangulation::coarse_triangulation(m);

    //Convert the mesh to a cell object
    cell_ptr c1 = initial_triangulation::convert_mesh_to_cell(m);

    //Generate the poisson point cloud of the surface, with a l_min of 0.1
    std::vector<oriented_point> poisson_point_cloud = initial_triangulation::generate_poisson_point_cloud(l_min, c1);

    //The following section has a n^2 complexity, so we cannot run the test with too many triangles
    for(auto& p: poisson_point_cloud){
        bool point_is_on_surface = false;
    
        //Loop over the faces of the mesh
        for(auto& f: m.face_point_ids){

            if(point_in_triangle(m, p.position_, f)){
                point_is_on_surface = true;
                break;
            }
        }
        if(!point_is_on_surface) return 1;
    }

    return 0;
}
//-------------------------------------------------------------------------------------------




//---------------------------------------------------------------------------------------------------------
// The main function
int main (int argc, char** argv){

    //Check that the command line input is correctly formatted
    assert(argc == 2); 

    //Get the name of the test to run
    std::string test_name = argv[1];
    
    if (test_name == "point_min_distance_test")     return point_min_distance_test();
    if (test_name == "point_on_surface_test")       return point_on_surface_test();

    


    std::cout << "TEST NAME :" << test_name << " DOES NOT EXIST" << std::endl;
    return 1;

}
//---------------------------------------------------------------------------------------------------------
