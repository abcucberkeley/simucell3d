#ifndef DEF_CELL_DIVIDER
#define DEF_CELL_DIVIDER


#include <cassert>
#include <iostream>
#include <vector>
#include <optional>

#include "utils.hpp"
#include "custom_structures.hpp"
#include "custom_exception.hpp"

#include "vec3.hpp"
#include "mat33.hpp"
#include "quaternion.hpp"

#include "cell.hpp"
#include "face.hpp"

#include "epithelial_cell.hpp"


#include "initial_triangulation.hpp"
#include "poisson_sampling.hpp"
#include "ball_pivoting_algorithm.hpp"
#include "local_mesh_refiner.hpp"


#include "delaunator.hpp" //Defined in lib/delaunator/include/delaunator.hpp


/*
    This class contains all the functions that are used to divide a cell into two daughter cells. Most 
    of the time is spent on retriangulating the daughter cells in a way that respects the constraints on
    the cell edge lengths.

    The methods of this class are static, so there is no need to instantiate the class, it only acts as a facade (see code pattern). 

*/



class cell_divider
{

    public:
        //Make sure the cell_divider cannot be default instantiated
        cell_divider() = delete;                                           //default constructor
        cell_divider(const cell_divider& c) = delete;           //copy constructor
        cell_divider(cell_divider&& c) = delete;                //move constructor
        cell_divider& operator=(const cell_divider& c) = delete;//copy assignment operator
        cell_divider& operator=(cell_divider&& c) = delete;     //move assignment operator 

        //If a cell is ready to be divided, this function will divide it into two daughter cells
        static void run(
            std::vector<cell_ptr>& cell_lst,
            const double l_min, //The min edge length
            const local_mesh_refiner& lmr, // Will be used to refine the meshes of the daughter cells
            unsigned& max_cell_id_, //Use this number to assign unique ids to the daughter cells
            bool verbose = true
        ) noexcept;

        //If the cell division fails, this function will return an empty optional
        static std::optional<std::pair<cell_ptr, cell_ptr>> divide_cell(
            cell_ptr c, 
            const double l_min, //The min edge length
            const local_mesh_refiner& lmr // Will be used to refine the meshes of the daughter cells
        ) noexcept; 

        //Find where the plane intersect with the edges of the cell and returns the mesh of the 
        //cell with the points where the plane has intersected with the edges
        static mesh add_intersection_points(
            cell_ptr c, 
            const vec3& p,  //Position of a point on the plane
            const vec3& n   //Normal vector of the plane
        )noexcept(false);


        //Find the intersection point between an edge (line segment) and a plane.
        //If there is no intersection returns an empty optional
        static std::optional<vec3> find_edge_plane_intersection(
            const vec3& e1, //Position of the 1st point of the edge
            const vec3& e2, //Position of the 2nd point of the edge

            const vec3& p,  //Position of a point on the plane
            const vec3& n   //Normal vector of the plane
        ) noexcept(false);


        //Based on the side of a face wrt the plane, return 1, -1 or 0
        static bool face_side_wrt_plane(
            const face& f,
            const cell_ptr c,
            const vec3& p,  //Position of a point on the plane
            const vec3& n   //Normal vector of the plane
        ) noexcept;

        //Add the id of the point p to the given face and make sure that 
        //the point p is inserted between the nodes n_a and n_b
        static void add_point_to_face(
            mesh& m, 
            const unsigned face_id, 
            const unsigned n_a_id, 
            const unsigned n_b_id,
            const unsigned n_p_id
        ) noexcept(false);


        //At this stage of the division process, some faces of the mesh have 5 nodes. This is because the plane
        //has intersected with 2 of their edges. This function splits those faces in 3 triangular subfaces. 
        //The argument "intersection_point_ids_threshold" corresponds to the id of the first intersection point
        //that was added to the mesh. The other intersection points that were added to the mesh have ids greater
        static void divide_faces(mesh& m, const unsigned intersection_point_ids_threshold) noexcept(false);

        //As it names indicates, this function maps the points in intersection_points to the xy plane
        static std ::pair<vec3, mat33> map_points_to_xy_plane(
            mesh& m, 
            const unsigned division_point_ids_threshold,
            const vec3& n
        ) noexcept;

        //Triangulate the interface between the 2 daughter cells in a way that respects the constraints on the edge lengths.
        static void triangulate_division_interface(
            const double l_min,                             //The minimum edge length
            mesh& m,                                        //The mesh containing all the faces and nodes
            const unsigned division_point_ids_threshold,    //All the nodes that belong to the division interface have an id greater than this threshold
            const unsigned division_face_ids_threshold,     //All the faces that belong to the division interface have an id greater than this threshold
            const vec3& division_face_normal                //The normal vector of the division interface
        ) noexcept(false);


        //Maps back the points of the division face to their original plane
        static void map_points_to_division_plane(
            mesh& m, 
            const unsigned division_point_ids_threshold,
            const vec3& translation,  
            const mat33& rotation_matrix
        ) noexcept;

        //Check if a point is inside a polygon
        static bool point_is_in_polygon(
            const std::vector<vec3>& polygon, 
            const vec3& point
        ) noexcept;

        //Split the mesh of the mother cell into 2 meshes, one for each daughter cell and returnt the daughter cell meshes
        static std::pair<cell_ptr, cell_ptr> create_daughter_cells(
            cell_ptr mother_cell, 
            mesh& mother_cell_mesh,                           //The mesh containing all the faces and nodes
            const unsigned division_point_ids_threshold,        //All the nodes that belong to the division interface have an id greater than this threshold
            const unsigned division_face_ids_threshold,         //All the faces that belong to the division interface have an id greater than this threshold
            const vec3& division_plane_normal,                  //The normal vector of the division interface
            const vec3& division_plane_origin                   // A point located on the division plane
        ) noexcept(false);




};

#endif