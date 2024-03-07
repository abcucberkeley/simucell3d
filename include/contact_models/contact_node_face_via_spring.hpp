#ifndef DEF_CONTACT_NODE_FACE_VIA_SPRING
#define DEF_CONTACT_NODE_FACE_VIA_SPRING


#include "global_configuration.hpp"


#if CONTACT_MODEL_INDEX == 0

     #include "contact_model_abstract.hpp"


     /*
     This class contains the methods to compute the contact forces between pairs of nodes and faces belonging to adjacent cells. 
     The contact forces are calculated in the following way: 

     1 -> All the axias-aligned bounding boxes of the faces are computed and stored in a vector
     2 -> The faces are stored in an uniform space partitioning grid (USPG) based on their AABB
     3 -> For each point, we used the USPG to find the faces that are close to it
     4 -> We check if the point is within the AABB of the face
     5 -> If so, we compute the minimal distance between the point and the face
     6 -> The contact force is the computed based on the minimum distance and the contact parameters of the face

     */


     class contact_node_face_via_spring: public contact_model_abstract{

     private:
          double interaction_cutoff_, interaction_cutoff_square_; 

          //The distance between 2 surfaces where the adhesion force is maximal, by default it's equal to interaction_cutoff_ / 2
          double hardening_distance_;

     public:
          contact_node_face_via_spring() = default;                                          //default constructor
          contact_node_face_via_spring(const contact_node_face_via_spring& v) = delete;           //copy constructor
          contact_node_face_via_spring(contact_node_face_via_spring&& v) = delete;                //move constructor
          contact_node_face_via_spring& operator=(const contact_node_face_via_spring& v) = delete;//copy assignment operator
          contact_node_face_via_spring& operator=(contact_node_face_via_spring&& v) = default;     //move assignment operator 

          //The trivial constructor of the contact model
          contact_node_face_via_spring(const global_simulation_parameters& sim_parameters) noexcept(false);

          //Find the contacts between the faces of the mesh and apply the contact forces
          void run(const std::vector<cell_ptr>& cell_lst) noexcept override;

          //Find the faces that are within a distance below the contact cutoff and apply adhesive or repulsive forces 
          void resolve_contacts(const std::vector<cell_ptr>& cell_lst) noexcept;

          //Compute the distance between the 2 faces and apply the contact forces
          void apply_contact_forces(cell_ptr c1, node& n, face* f) const noexcept;

     };

     #endif

#endif