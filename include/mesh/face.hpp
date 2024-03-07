#ifndef DEF_FACE
#define DEF_FACE

#include <cassert>
#include <map>
#include <iostream>
#include <optional>
#include <memory>
#include <omp.h>

class node;
class cell;

#include "global_configuration.hpp"

#include "vec3.hpp"
#include "node.hpp"



typedef std::shared_ptr<cell> cell_ptr;

/*
    Class for a triangular face of the cell mesh
*/

class face final
{

    private: 
        vec3 normal_;                    // Face outward normal
        double area_ = -1.;              // Face area

        //The id of the face in the cell lst
        unsigned local_face_id_;

        //The id of the face in the mesh. This global id is only used
        //by the contact detection modules
        unsigned global_face_id_;

        //A pointer toward the cell that owns this face
        cell_ptr owner_cell_ = nullptr;

        //The type of the face
        unsigned short type_id_ = 0;

        //Only compiles this section if the faces store their contact energies
        #if FACE_STORE_CONTACT_ENERGY
            //Store the adhesion and repulsion energies associated to this face, float are good enough
            //since these values are only used for plotting
            float adhesion_energy_ = 0.;
            float repulsion_energy_ = 0.;
        #endif

        //cell is a friend class because its the only one that can reset the face id
        friend class cell;
        friend class local_mesh_refiner;

        //All the contact model classes have direct access to the face attributes 
        //to compute more efficiently the contact forces
        friend class contact_model_abstract;
        friend class contact_node_node_via_coupling;
        friend class contact_face_face_via_coupling;
        friend class contact_node_face_via_spring;

        friend class cell_divider;
        friend class face_tester;
        friend class cell_tester;



    //All this functions are protected to make sure that only the cell owning this face
    //can call them.
    protected:

        //The id of the node composing the face
        //A face might be left in a state where it is not connected to 
        //any node.
        unsigned n1_id_ = 0;
        unsigned n2_id_ = 0;
        unsigned n3_id_ = 0; 


        //Indicate if the face is used by the cell to store nodes
        bool is_used_ = true;

        //The ids of the nodes might change, so we need to update them in each face
        void update_node_ids(std::map<unsigned, unsigned>& node_id_correspondence) noexcept;

        //Only the cell should modify the local ids of its faces
        void set_local_id(const unsigned id) noexcept {local_face_id_ = id;}; 
   
        //Empties the face content
        void reset() noexcept;

        //Swap the second and third node to correctly orient the face normal
        void swap_nodes() noexcept;

        //Only the cell can set the area of a face
        void set_area(double area) noexcept;
        void set_normal(const vec3& normal) noexcept{normal_ = normal;}


        //Only compiles this section if the faces store their contact energies
        #if FACE_STORE_CONTACT_ENERGY
            void add_adhesion_energy(float adhesion_energy) noexcept {
                #pragma omp atomic update
                adhesion_energy_ += adhesion_energy;
            }

            void add_repulsion_energy(float repulsion_energy) noexcept {
                #pragma omp atomic update
                repulsion_energy_ += repulsion_energy;
            }
        #endif




    public:

        //Make sure a face cannot be instantiated without its nodes
        face() = default;
        face(const face& f) = default;            //copy constructor
        face(face&& f) = default;                  //move constructor
        face& operator=(const face& f) = default; //copy assignment operator
        face& operator=(face&& f) = default;      //move assignment operator 

        //Instantiate a face with the unique id of its nodes
        explicit face(const unsigned node_id_1, const unsigned node_id_2, const unsigned node_id_3, const unsigned face_id) noexcept;

        //Instantiate a face with the reference of its nodes 
        explicit face(const node& n1, const node& n2, const node& n3, const unsigned face_id) noexcept;

        //Get the unique id of this face
        const unsigned get_local_id() const noexcept {return local_face_id_;};
 
        //Return the face area
        inline double get_area() const noexcept {return area_;}

        //Set the type of the face 
        void set_face_type_id(const unsigned short type_id) noexcept {

            //Protect against race conditions
            #pragma omp atomic write
            type_id_ = type_id;
        };

        //WARNING: CAN RETURN A NULL POINTER
        cell_ptr get_owner_cell() const noexcept {return owner_cell_;}

        inline const vec3& get_normal() const noexcept {return normal_;}

        //Indicate if nodes are using this face
        bool is_used() const noexcept {return is_used_;}
        void set_is_used(bool is_used) noexcept {is_used_ = is_used;}


        bool has_node(const unsigned node_id) const noexcept{
            return (node_id == n1_id_ || node_id == n2_id_ || node_id == n3_id_);
        }

        //Get the nodes of the face in a array
        std::array<unsigned, 3> get_node_ids() const noexcept; 

        //Return the ids of the nodes of the face
        unsigned n1_id() const noexcept {assert(is_used_); return n1_id_;}
        unsigned n2_id() const noexcept {assert(is_used_); return n2_id_;}
        unsigned n3_id() const noexcept {assert(is_used_); return n3_id_;}


        //Replace one node of the face by another one
        void replace_node(const unsigned old_node_id, const unsigned new_node_id) noexcept;

        //Get the id of the node in the triangle opposite to the given edge
        unsigned get_opposite_node(const unsigned n1, const unsigned n2) const noexcept;

        //This is the position of the face type in the cell_type.face_types_ vector
        unsigned short get_local_face_type_id() const noexcept {return type_id_;}
      
        //Only compiles this section if the faces store their contact energies
        #if FACE_STORE_CONTACT_ENERGY
            float get_adhesion_energy()  const noexcept {return adhesion_energy_;}
            float get_repulsion_energy() const noexcept {return repulsion_energy_;}
        #endif

        bool operator==(const face& f) const noexcept; 
        bool operator!=(const face& f) const noexcept;

};

#endif