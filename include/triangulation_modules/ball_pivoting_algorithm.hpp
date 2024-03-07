#ifndef DEF_BALL_PIVOTING_ALGORITHM
#define DEF_BALL_PIVOTING_ALGORITHM

#define _USE_MATH_DEFINES

#include <map>
#include <random>
#include <optional>
#include <unordered_set>

#include "custom_exception.hpp"
#include "custom_structures.hpp"

#include "utils.hpp"
#include "vec3.hpp"
#include "../uspg/uspg_4d.hpp"

#include "cell.hpp"
#include "face.hpp"
#include "edge.hpp"
#include "node.hpp"

#include "mesh_writer.hpp"





/*
    Given a poisson point cloud, this class reconstructs a triangulated surface by using the 
    ball pivoting algorithm. The theory behind this algorithm can be found in the original paper:
    "The ball-pivoting algorithm for surface reconstruction", by Bernardini et al, DOI: 10.1109/2945.817351
*/



class ball_pivoting_algorithm
{
    private:
        double ball_radius_;
        double ball_radius_squared_;

        std::vector<oriented_point> node_lst_;

        //Store all the points in a uspg_4D
        uspg_4d<oriented_point> grid_;

        //All the created edges and faces are stored in vectors
        std::vector<edge> edge_lst_;
        std::vector<face> face_lst_;

        //Some faces will have to be deleted during the hole filling process. The indexes of the faces
        //to be deleted are stored in this vector. This faces can only be deleted at the end of the triagulation
        //std::vector<unsigned> face_to_delete_lst_;
        //std::vector<unsigned> edge_to_delete_lst_;


        //The nodes keep track of the edges or faces in which they are involved. All this tracking is 
        //needed to determine if a node is already in the triangulated region of the cell surface or not.
        //The key of the map is the node_id, the value is the list of the indexes of the faces using the node
        std::map<unsigned, std::vector<unsigned>> node_face_lst_; 
        std::map<unsigned, std::vector<unsigned>> node_edge_lst_;

        //Keep track of the edges at the front of the trangulation
        std::vector<unsigned> edge_id_front_lst_;

        //The tester class used to test the methods of this bpa class
        friend class bpa_tester;

        //This constructor is only for testing purposes
        ball_pivoting_algorithm(const double l_min): ball_radius_(1.2 * l_min), ball_radius_squared_(1.44 * l_min * l_min){};
        

    public:
        //Make sure the ball_pivoting_algorithm cannot be default instantiated                                        
        ball_pivoting_algorithm() = delete;                                           //default constructor
        ball_pivoting_algorithm(const ball_pivoting_algorithm& c) = delete;           //copy constructor
        ball_pivoting_algorithm(ball_pivoting_algorithm&& c) = delete;                //move constructor
        ball_pivoting_algorithm& operator=(const ball_pivoting_algorithm& c) = delete;//copy assignment operator
        ball_pivoting_algorithm& operator=(ball_pivoting_algorithm&& c) = delete;     //move assignment operator 

        //The trivial construtor
        ball_pivoting_algorithm(
            const std::vector<oriented_point>& poisson_point_cloud, 
            const double l_min,
            const double min_x, const double min_y, const double min_z, 
            const double max_x, const double max_y, const double max_z 
        ) noexcept(false);

        //Run the algorithm
        void run() noexcept(false){
            find_seed_triangle();
            expand_triangulation();
            fill_surface_holes();
        }

        //Find a triangle from which the surface expansion can be started
        void find_seed_triangle() noexcept(false);

        //Getters, return copies (it's expensive but only run once at the initialization)
        std::vector<face> get_face_lst() const noexcept{return face_lst_;}
        std::vector<edge> get_edge_lst() const noexcept{return edge_lst_;}
        std::vector<oriented_point> get_node_lst() const noexcept{return node_lst_;}

        void expand_triangulation() noexcept;

        //Check if the points have normals pointing in the same direction
        static bool points_coherently_oriented(
            const oriented_point& a, 
            const oriented_point& b,
            const oriented_point& c
        ) noexcept;
        
        //Compute the normal of the plane formed by the 3 points and make 
        //that this normal is oriented in the same direction as the points' normals
        //The retruned normal is a unit vector
        static vec3 compute_oriented_normal(
            const oriented_point& point_a, 
            const oriented_point& point_b,
            const oriented_point& point_c
        ) noexcept;

        //Get the center of the sphere of radius = ball_radius_ passing by the 3 points
        //If this method doesn't return anything then this circumsphere does not exist
        std::optional<vec3> compute_circum_sphere_center(
            const oriented_point& point_a, 
            const oriented_point& point_b,
            const oriented_point& point_c
        ) const noexcept;

        std::optional<vec3> compute_circum_sphere_center(
            const face& f
        ) const noexcept;

        //Check if the circum_sphere with center circum_sphere_center and radius ball_radius_
        //does not overlap with points that are point_a, point_b, or point_c..
        bool circum_sphere_is_empty(
            const vec3& circum_sphere_center, 
            const oriented_point& point_a, 
            const oriented_point& point_b, 
            const oriented_point& point_c
        ) const noexcept;

        //Returns the id of the edge formed by 2 points.
        unsigned get_edge(const unsigned point_1_id, const unsigned point_2_id) noexcept;
        unsigned get_edge(const oriented_point& point_a, const oriented_point& point_b) noexcept;

        //Returns the id of the best node with which to connect the edge to form a triangle
        std::optional<oriented_point> find_candidate_node(const edge& edge) const noexcept;

        //Returns if the point is already part of the triangulated region of the cell surface
        bool node_is_inner_vertex(const unsigned point_id)  const noexcept;
        bool node_is_inner_vertex(const oriented_point& p) const noexcept;

        //Make sure that the normal of the face {f1_a, f1_b, p} points in the same direction as the normal
        //of the face {f1_a, f1_b, f1_c}
        static bool node_is_compatible_with_edge(
            const oriented_point& p, 
            const oriented_point& f1_a,
            const oriented_point& f1_b,
            const oriented_point& f1_c
        ) noexcept;

        //Compute the angle between the center of the circumspheres and the center of the edge
        //the candidate point, with the smallest angle is the one choosen. Basically the 
        //algorithm tries to create a flat triangulated surface. The returned angle should be in the 
        //range [0, 2pi]
        static double compute_angle_between_circumspheres(
                const vec3& face_a_circum_center,
                const vec3& face_b_circum_center,
                const vec3& edge_center_vec
            ) noexcept;



        //Fill the remaining holes in the surface with triangles
        void fill_surface_holes() noexcept(false);

        //Get the next node which is part of the hole
        unsigned get_next_node_of_the_hole(const unsigned previous_node, const unsigned current_node) const noexcept(false);





};

#endif