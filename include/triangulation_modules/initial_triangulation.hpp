#ifndef DEF_INITIAL_TRIANGULATION
#define DEF_INITIAL_TRIANGULATION


#include <cassert>
#include <iostream>
#include <vector>

#include "utils.hpp"
#include "custom_structures.hpp"
#include "custom_exception.hpp"
#include "vec3.hpp"
#include "cell.hpp"

#include "poisson_sampling.hpp"
#include "ball_pivoting_algorithm.hpp"

/*
    This class has no constructor and all its methods are static. They can therefore
    be called without instantiating the class. 

    Given an arbitrary surface, this class returns a corresponding triangulated surface
    where the lengths of the edges are comprised in the range [l_min, l_max]. This class 
    does it by following the given procedure:

    1 - Coarsely triangulate the surface, by connecting the edges of each face to the face center
    2 - Sample the newly triangulated surface with the Poisson Disc Sampling algorithm
    3 - Reconstruct a finely triangulated surface with the Ball Pivoting Algorithm
    4 - Call the edge_manager to make sure that the lengths of the edges are in the prescribed range

    This class will throw an error if it receives a surface that is not closed

*/



class initial_triangulation
{

    public:
        //Make sure the initial_triangulation cannot be default instantiated
        initial_triangulation() = delete;                                           //default constructor
        initial_triangulation(const initial_triangulation& c) = delete;           //copy constructor
        initial_triangulation(initial_triangulation&& c) = delete;                //move constructor
        initial_triangulation& operator=(const initial_triangulation& c) = delete;//copy assignment operator
        initial_triangulation& operator=(initial_triangulation&& c) = delete;     //move assignment operator 

        //The main function from which the surface is triangulated
        static mesh triangulate_surface(
            double const l_min, 
            double const l_max, 
            const mesh& surface, 
            const unsigned cell_id
        ) noexcept(false); 

        //Coarsely triangulate the mesh by connecting the edges of each face to the center of the face
        static void coarse_triangulation(mesh& surface) noexcept; 

        static cell_ptr convert_mesh_to_cell(const mesh& surface) noexcept(false); //throw mesh_integrity_exception

        static std::vector<oriented_point> generate_poisson_point_cloud(const double l_min, const cell_ptr c) noexcept;

};

#endif