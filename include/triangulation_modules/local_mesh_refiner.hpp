#ifndef DEF_local_mesh_refiner
#define DEF_local_mesh_refiner

#define _USE_MATH_DEFINES


#include <cassert>
#include <cmath>  
#include <set>  


#include "utils.hpp"


#include "custom_exception.hpp"
#include "custom_structures.hpp"

#include "mesh_writer.hpp"
#include "cell.hpp"
#include "face.hpp"
#include "edge.hpp"



/*
    This class makes sure that all the edges of the cell meshes have lengths comprised in the range [l_min, l_max]. 
    When an edge is too long, it is subdivided into two edges. When an edge is too short, it is merged with into one node.
    Triangles with very high isoperimetric ratio are also subdivided into two triangles.
*/

class local_mesh_refiner
{

    private:

        //The minimum and maximum authorized edge lengths
        const double l_min_, l_max_, l_min_squared_, l_max_squared_;

        //If set to true, the local mesh refiner will remove the triangles with high aspect ratio
        //by swapping their longest edge
        const bool enable_edge_swap_operation_;

        //The isoperimetric ratio of an equliateral triangle, which
        //is the smmallest possible isoperimetric ratio
        static constexpr double q_min_ = 36. / std::sqrt(3.);

        //The minimum allowed triangle score, below that, the edge swap operation is triggered
        static constexpr double triangle_score_min_ = 0.2;

        //Wrap the call to the local_mesh_refiner::refine_mesh() method into a lambda function
        const std::function<void(cell_ptr)> refine_mesh_func_ = [=](cell_ptr c) -> void {refine_mesh(c);};

        friend class local_mesh_refiner_tester;


    public:
        //Make sure the local_mesh_refiner cannot be default instantiated
        local_mesh_refiner() = default;                                     //default constructor
        local_mesh_refiner(const local_mesh_refiner& c) = delete;           //copy constructor
        local_mesh_refiner(local_mesh_refiner&& c) = delete;                //move constructor
        local_mesh_refiner& operator=(const local_mesh_refiner& c) = default;//copy assignment operator
        local_mesh_refiner& operator=(local_mesh_refiner&& c) = default;     //move assignment operator 

        double get_l_min() const noexcept {return l_min_;};
        double get_l_max() const noexcept {return l_max_;};

        double get_l_min_squared() const noexcept {return l_min_squared_;};
        double get_l_max_squared() const noexcept {return l_max_squared_;};

        //Trivial constructor
        local_mesh_refiner(const double l_min, const double l_max, const bool enable_edge_swap_operation = true) noexcept;

        //Refine the meshes of all the cells in the vector
        void refine_meshes(const std::vector<cell_ptr> cell_lst) const noexcept(false);

        //Refine the mesh of a given cell
        void refine_mesh(cell_ptr c) const noexcept(false);

        //Remove all the elongated triangles
        void remove_elongated_triangles(cell_ptr c) const noexcept(false);

        //Return the triangle score and also the longest edge of the triangle
        std::pair<double, edge> get_triangle_score(cell_ptr c, const face& f) const noexcept;

        //Check if an edge can be merged into a node
        bool can_be_merged(edge& e_ab, cell_ptr c) const noexcept(false);
 
        void swap_edge(edge& e_ab, cell_ptr c) const noexcept(false);

        void split_edge(edge& e_ab, cell_ptr c, edge_set& edge_to_check_set) const noexcept(false);

        void merge_edge(edge& e_ab, cell_ptr c, edge_set& edge_to_check_set) const noexcept(false);


        




};

#endif