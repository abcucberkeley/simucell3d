#ifndef DEF_CONTACT_MODEL_ABSTRACT
#define DEF_CONTACT_MODEL_ABSTRACT


#include <iostream>
#include <vector>
#include <array>
#include <omp.h>

#include "global_configuration.hpp"

#include "utils.hpp"
#include "custom_structures.hpp"
#include "custom_exception.hpp"

#include "cell.hpp"
#include "face.hpp"
#include "uspg_4d.hpp"





//--------------------------------------------------------------------------------------------------------------------------------
//This is an abstract class that incorporate all the attributes and methods common to the different contact models
class contact_model_abstract{

        private:
                friend class tester_contact_model_abstract;

        protected:

                //The length size of each voxel used to instantiate the grid
                double voxel_size_;
                
                //The cutoff distance above wich the interaction forces between 2 faces are not calculated anymore
                double interaction_cutoff_repulsion_, interaction_cutoff_square_repulsion_; 
                double interaction_cutoff_adhesion_, interaction_cutoff_square_adhesion_; 
                double max_interaction_cutoff_square_;

                
                //Keep track of the min and max coordinates of the mesh to update the partitionning grid
                double global_min_x_, global_min_y_, global_min_z_;
                double global_max_x_, global_max_y_, global_max_z_;

                //The uniform space partitionning grid that will store the faces of the mesh
                uspg_4d<face*> grid_;

                //Store the Axis Aligned Bounding Boxes of the face in this vector
                std::vector<double> face_aabb_lst_;

                //Store pointers to all the faces of the mesh in this vector
                std::vector<face*> face_lst_;

                //The extra padding added to the AABB of the faces, to check for potential contacts
                double aabb_padding_;

        public:
                contact_model_abstract() = default;                                         //default constructor
                contact_model_abstract(const contact_model_abstract& v) = delete;           //copy constructor
                contact_model_abstract(contact_model_abstract&& v) = delete;                //move constructor
                contact_model_abstract& operator=(const contact_model_abstract& v) = delete;//copy assignment operator
                contact_model_abstract& operator=(contact_model_abstract&& v) = default;    //move assignment operator 

                //The trivial constructor of the contact model
                contact_model_abstract(const global_simulation_parameters& sim_parameters) noexcept(false);

                //Find the contacts between the faces of the mesh and apply the contact forces
                virtual void run(const std::vector<cell_ptr>& cell_lst) noexcept = 0;

                //Create the axis aligned bounding box of the faces
                void update_face_aabbs(
                        const std::vector<cell_ptr>& cell_lst
                ) noexcept;

                //Store the faces in the uniform space partitionning grid
                void store_face_in_uspg() noexcept;

                //Check if the node is in the AABB of the face
                bool aabb_intersection_check(const size_t face_aabb_pos, const vec3& node_pos) const noexcept;

                //Find the closest point on a triangle to a given point, and return as well
                //the minimal squared distance between the point and the triangle
                static std::pair<double, vec3> compute_node_triangle_distance(
                        const vec3& p, 
                        const vec3& a, 
                        const vec3& b, 
                        const vec3& c
                ) noexcept;

};
//--------------------------------------------------------------------------------------------------------------------------------








#endif