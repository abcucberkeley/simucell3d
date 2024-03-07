#ifndef DEF_NODE
#define DEF_NODE


#include <iostream>
#include <memory>
#include <optional>
#include <map>

#include <omp.h>

#include "global_configuration.hpp"


class face;
class cell;

#include "vec3.hpp"


/*
    Class for a node of the cell mesh
*/


typedef std::shared_ptr<cell> cell_ptr;

class node final 
{

    private: 



        //Each node has a unique identifier
        unsigned node_id_;

        //Wether or not it is used by a cell
        bool is_used_ = true;

        //cell is a friend class because its the only one that can reset the node id
        friend class cell;
        friend class epithelial_cell;
        friend class local_mesh_refiner;
        friend class time_integration_scheme;
        friend class dynamic_motion_semi_implicit_euler;       //The time integration scheme
        friend class overdamped_motion_forward_euler;          //The time integration scheme
        friend class overdamped_motion_improved_euler;         //The time integration scheme
        
        //The class that computes the contact forces between the cells
        friend class contact_model_abstract;
        friend class contact_node_node_via_coupling;
        friend class contact_face_face_via_coupling;
        friend class contact_node_face_via_spring;

        friend class cell_divider;
        friend class node_tester;
        friend class cell_tester;


        
    protected:

        ///If we use a contact model that mechanically link the nodes or the faces of adjacent cells
        #if CONTACT_MODEL_INDEX == 1 ||  CONTACT_MODEL_INDEX == 2
            //The normal of the node
            vec3 normal_;

            //Curvature of the cell surface at the node
            double curvature_;

            //Create a lock
            omp_lock_t lock_;


        #endif


        //If we link pairs of adjacent nodes
        #if CONTACT_MODEL_INDEX == 1

            //This the node that is coupled to this node. A node can only be coupled to one other node
            std::optional<std::pair<unsigned, unsigned>> coupled_node_;

            //The squared distance to the closest node on another cell
            double squared_distance_to_closest_node_ = 0.;
        
        //If we link pairs of adjacent faces
        #elif CONTACT_MODEL_INDEX == 2

            //Each node can be coupled to one node on each adjacent cell. I.E. a node can be coupled to multiple nodes
            //if they belong to different cells

            //The key of this map is the local id of the adjacent cell
            //The value is the local id of the coupled node as well as the distance to this coupled node
            std::map<unsigned, std::pair<unsigned, double>> coupled_nodes_map_;

        #endif


        //Based on the dynamic model chosen to run the simulation the nodes need to store different state variables
        vec3 pos_, force_;

        #if DYNAMIC_MODEL_INDEX == 0 //Full equations of motion solved with the semi implicit euler scheme
            vec3 momentum_;
        #endif 


        void set_local_id(const unsigned node_id) noexcept {node_id_ = node_id;}

        //Reset the vectors owned by the node
        void reset() noexcept;

        //Indicate wether or not the node is used by a cell
        void set_is_used(const bool is_used) noexcept {is_used_  = is_used;}

    public:

        //Default constructor, everything is set to 0
        node() = default;
        node(const node& n) = default;           //copy constructor
        node(node&& n) = default;                //move constructor
        node& operator=(const node& n) = default;//copy assignment operator
        node& operator=(node&& n) = default;     //move assignment operator 

        //Leave the position of the node to (0., 0., 0.)
        explicit node(const unsigned node_id) noexcept : node_id_(node_id){};

        //Specify node coordinates at instantiation
        explicit node(const double dx, const double dy, const double dz, const unsigned node_id) noexcept;

        explicit node(const vec3& pos, const unsigned node_id) noexcept;
        
        //Get the unique identifier of this node
        unsigned get_local_id() const noexcept {return node_id_;}

        //Get a constant reference of the node position, force and momentum
        inline const vec3& pos()       const noexcept {return pos_;}
        inline const vec3& force()     const noexcept {return force_;}
        inline const bool  is_used()   const noexcept {return is_used_;}

        //If we use a contact model that mechanically link the nodes or the faces of adjacent cells
        #if CONTACT_MODEL_INDEX == 1 ||  CONTACT_MODEL_INDEX == 2

            double get_curvature() const noexcept {return curvature_;}

            vec3 get_normal() const noexcept {return normal_;}

            // Destroy the lock in the destructor
            ~node(){omp_destroy_lock(&lock_);}

        #endif

        //If we link pairs of adjacent nodes
        #if CONTACT_MODEL_INDEX == 1

            bool is_coupled() const noexcept {return coupled_node_.has_value();}

            std::pair<unsigned, unsigned> get_coupled_node() const noexcept {return coupled_node_.value();}

            //Indicate to this node that it is coupled to another node
            void set_coupled_node_and_min_distance(const std::pair<unsigned, unsigned>& coupled_node, const double min_dist) noexcept {

                //Make sure that one thread at the time can access this function
                omp_set_lock(&lock_);

                coupled_node_ = coupled_node;
                squared_distance_to_closest_node_ = min_dist;

                omp_unset_lock(&lock_); 
            }

        //If we link pairs of adjacent faces
        #elif CONTACT_MODEL_INDEX == 2

            //Return wether or not this node is coupled to one or more nodes on other cells
            bool is_coupled() const noexcept {return !coupled_nodes_map_.empty();}

            //Get the number of nodes on other cells that are coupled to this node
            int get_nb_coupled_nodes() const noexcept {return coupled_nodes_map_.size();}

            void set_coupled_node_and_min_distance(const unsigned coupled_cell_local_id, const unsigned coupled_node_local_id,  const double squared_distance) noexcept {
                /*
                Store the information that this node is coupled to a node with the local id `coupled_node_local_id` on the cell with the local id `coupled_cell_local_id`.
                */
                
                //Make sure that one thread at the time can access this function
                omp_set_lock(&lock_);

                //Check if this node has already been coupled to a node on the second cell
                auto it = coupled_nodes_map_.find(coupled_cell_local_id);

                //If this node has already been coupled to a node on the second cell
                if(it != coupled_nodes_map_.end()){

                    //Change the node to which this node is coupled on the second cell
                    it->second = std::make_pair(coupled_node_local_id, squared_distance);
                }

                //If this node has not yet been coupled to nodes on the cell c2
                else{
                    coupled_nodes_map_[coupled_cell_local_id] = std::make_pair(coupled_node_local_id, squared_distance);
                }
                omp_unset_lock(&lock_); 
            }
        #endif




        //Add a force on the node
        void set_force(const vec3& force) noexcept {force_.reset(force);}
        void add_force(const vec3& force) noexcept {force_.translate(force);}

        //Overload the comparion operators
        vec3 operator-(const node& n) const noexcept; 
        bool operator==(const node& n) const noexcept; 
        bool operator!=(const node& n) const noexcept;

        //Full equations of motion solved with the semi implicit euler scheme
        #if DYNAMIC_MODEL_INDEX == 0 

            //Methods to get and set the momentum of the node
            inline const vec3& momentum()  const noexcept {return momentum_;}
            void set_momentum(const vec3& momentum) noexcept {momentum_ = momentum;}
        #endif 


};

#endif