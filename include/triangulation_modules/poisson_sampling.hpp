#ifndef DEF_POISSON_SAMPLING
#define DEF_POISSON_SAMPLING

#define _USE_MATH_DEFINES


#include <cassert>
#include <chrono> 
#include <cmath>  
#include <random> 
#include <vector>
#include <forward_list>

#include "custom_exception.hpp"
#include "custom_structures.hpp"

#include "utils.hpp"
#include "vec3.hpp"

//Why the linker fails if you don't give the uspg absolute path? This is a mistery
//beyond the reach of human knowledge
#include "../uspg/uspg_3d.hpp"
#include "../uspg/uspg_4d.hpp"

#include "cell.hpp"
#include "face.hpp"





/*
    Given a triangulated surface, this class returns a point cloud where the points are located on the surface
    and are separated from each other by at least a given distance of l_min. The theory on the technique can be 
    found in the paper: "Parallel Poisson disk sampling with spectrum analysis on surfaces" Bowers et al, SIGGRAPH ASIA 2010.
    DOI: 978-1-4503-0439-9.

    The class is a static class and cannot be istantiated
*/

class poisson_sampling
{



    public:
        //Make sure the poisson_sampling cannot be default instantiated
        poisson_sampling() = delete;                                           //default constructor
        poisson_sampling(const poisson_sampling& c) = delete;           //copy constructor
        poisson_sampling(poisson_sampling&& c) = delete;                //move constructor
        poisson_sampling& operator=(const poisson_sampling& c) = delete;//copy assignment operator
        poisson_sampling& operator=(poisson_sampling&& c) = delete;     //move assignment operator 

        //Return a poisson point cloud of the cell surface
        static std::vector<oriented_point> compute_poisson_point_cloud(double const l_min, const cell_ptr input_cell) noexcept(false); 

        //The input mesh should already be triangulated, this version of the poisson point cloud 
        //is used during the cell division. It makes sure that the nodes of the mesh that are located at the 
        //division interface will be included in the poisson point cloud 
        static std::vector<oriented_point> compute_poisson_point_cloud(
            double const l_min, 
            mesh& m,
            const unsigned division_point_ids_threshold,  // All the points with an id equal or greater than this value should be added to the sampling
            const unsigned division_face_ids_threshold,   //All the faces with an id equal or greater than this value will be sampled
            const vec3 division_plane_normal //The normal of the division plane
        ) noexcept(false);

        //Uniformly sample the cell surface
        static std::vector<oriented_point> uniform_sampling(cell_ptr input_cell, const double l_min) noexcept(false); 
        static std::vector<oriented_point> uniform_sampling(const mesh& m, const unsigned division_face_ids_threshold, const double l_min) noexcept(false);

        //Carry out the Poisson Disk Sampling of the cell surface based on the preceeding uniform sampling
        static std::vector<oriented_point> poisson_disk_sampling(
            const uspg_4d<oriented_point>& grid_1, 
            uspg_4d<oriented_point>& grid_2,
            const double l_min
        ) noexcept;

        //Utility function to get the AABB of the set of faces of a mesh
        static std::array<double, 6> get_surface_aabb(
            const mesh& m, 
            const std::vector<unsigned> face_to_sample_ids
        ) noexcept;

        //Utility function to get the AABB of the set of faces of a mesh
        static std::array<double, 6> get_surface_aabb(
            const mesh& m, 
            const unsigned division_face_ids_threshold // All the faces with id >= division_face_ids_threshold will be included in the AABB
        ) noexcept;
};

#endif