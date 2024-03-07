#ifndef DEF_TEST_BPA
#define DEF_TEST_BPA

#include <cassert>
#include <string>
#include <iostream>
#include <vector>

#include "utils.hpp"
#include "custom_structures.hpp"
#include "ball_pivoting_algorithm.hpp"
#include "mesh_writer.hpp"


//Class to test the method of the Ball pivoting algorithm
class bpa_tester
{   

    public:
        bpa_tester() = default;

        int compute_oriented_normal_test() const;
        int points_coherently_oriented_test() const;
        int node_is_compatible_with_edge_test() const;
        int compute_angle_between_circumspheres_test() const;
        int compute_circum_sphere_center_test() const;
        int circum_sphere_is_empty_test() const;
        int get_edge_test() const;
        int node_is_inner_vertex_test() const;
        int fill_surface_holes_test() const;

};

#endif