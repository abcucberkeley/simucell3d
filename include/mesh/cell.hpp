#ifndef DEF_CELL
#define DEF_CELL

#define _USE_MATH_DEFINES

#include <cassert>
#include <cmath>
#include <vector>
#include <map>
#include <iostream>
#include <algorithm>
#include <type_traits>
#include <memory>
#include <list>
#include <set>
#include <numeric>
#include <chrono> 
#include <random>

#include "global_configuration.hpp"

class face;
class edge;
class node;


#include "custom_structures.hpp"
#include "custom_exception.hpp"
#include "node.hpp"
#include "edge.hpp"
#include "face.hpp"
#include "vec3.hpp"
#include "mat33.hpp"



/*
    The cell class is a mediator (see mediator pattern design). It manages the relationships between the nodes, the edges
    and the faces and make sure that the integrity of the cell surface is maintained. This limits
    the dependencies between the nodes and faces, and therefore reduces the memory usage.
*/

//typedef std::unordered_set<edge, edge_hasher,  edge_comparator> edge_set;
typedef std::set<edge> edge_set;



typedef std::shared_ptr<cell_type_parameters> cell_type_param_ptr;

class cell: public std::enable_shared_from_this<cell> {

    protected: 
        
        //Store all the nodes in an vector
        std::vector<node> node_lst_;

        //Store all the faces in an vector
        std::vector<face> face_lst_;

        //The centroid of the cell which is equal to the center of mass
        vec3 centroid_;

        double area_ = 0.;   //  Cell area
        double volume_ = 0.; //  Cell volume
        double pressure_ = 0.; //Pressure difference between the inside and outside of the cell
        double target_volume_ = 0.; //Target volume of the cell
        double target_area_ = 0.; //Target area of the cell

        //Each cell has a unique identifier, note that the cell_id might not be the position of the cell
        //in the cell_lst vector. 
        unsigned cell_id_;

        //The index of the cell in the cell_lst vector
        unsigned local_id_; 

        //Pointer to the structure that stores all the cell type parameters
        cell_type_param_ptr cell_type_; 

        //Stores in this vector the IDs of all the nodes or faces that are not used to represent the mesh of the cell
        std::vector<unsigned> free_node_queue_;
        std::vector<unsigned> free_face_queue_;

        //The different contributions to the potential energy of the cell
        //We use float values since they are only meant for plotting
        float surface_tension_energy_ = 0.;
        float membrane_elasticity_energy_ = 0.;
        float bending_energy_ = 0.;
        float pressure_energy_ = 0.;
        float kinetic_energy_ = 0.;

        //The growth rate and division volumes of the cells are drawn from a normal distribution
        double growth_rate_ = 0.;
        double division_volume_ = 0.;

        //If this flag is set to true, the cell will not move
        bool is_static_ = false; 
        
        //Store all the edges of the cell in this unordered set
        edge_set edge_set_; 

        //We use a lot of friend classes too directly acces the attributes of the cell class
        //and thus avoid the use of getters and setters which would slow down the code.
        friend class node;
        friend class time_integration_scheme;
        friend class dynamic_motion_semi_implicit_euler;       //The time integration scheme
        friend class overdamped_motion_forward_euler;             //The time integration scheme
        friend class overdamped_motion_improved_euler;             //The time integration scheme
        friend class contact_model_abstract;
        friend class contact_node_node_via_coupling;
        friend class contact_face_face_via_coupling;
        friend class contact_node_face_via_spring;
        friend class cell_divider;
        friend class automatic_polarizer;

        //The following classes are used for testing purposes
        friend class cell_tester;
        friend class node_tester;
        friend class local_mesh_refiner;
        friend class local_mesh_refiner_tester;

 
        void add_free_node(const unsigned node_id) noexcept;
        void add_free_face(const unsigned face_id) noexcept;

        void update_target_volume(const double time_step) noexcept;

        //Apply the pressure on the nodes of the cell mesh
        void apply_pressure_on_surface() noexcept;

        //Apply at the same time the surface tension forces and the membrane elasticity forces
        void apply_surface_tension_and_membrane_elasticity() noexcept;

        void apply_bending_forces() noexcept;


        node& get_node(const unsigned node_id) noexcept{
            assert(node_id < node_lst_.size()); 
            assert(node_lst_[node_id].is_used());
            return node_lst_[node_id]; 
        }
        face& get_face(const unsigned face_id) noexcept{
            assert(face_id < face_lst_.size()); 
            assert(face_lst_[face_id].is_used());
            return face_lst_[face_id]; 
        }
        
        //Translate the whole cell
        void translate(const vec3& translation) noexcept;
       
    public:

        //Make sure the cell cannot be default instantiated
        cell() = delete;                         //default constructor
        cell(const cell& c) = default;           //copy constructor
        cell(cell&& c) = delete;                 //move constructor
        cell& operator=(const cell& c) = delete; //copy assignment operator
        cell& operator=(cell&& c) = delete;      //move assignment operator 


        
        //Construct the cell with a node_lst and a face_lst
        cell(
            const std::vector<node>& node_lst,
            const std::vector<face>& face_lst,
            unsigned cell_id, 
            cell_type_param_ptr cell_type_parameters = nullptr
        ) noexcept;



        //Construct the cell based on the node position and the node  ids of the faces
        cell(
            const std::vector<double>&  node_position,
            const std::vector<unsigned>& face_node_ids,
            unsigned cell_id,
            cell_type_param_ptr cell_type_parameters = nullptr
        ) noexcept;


        //Convert a mesh object into a cell object
        cell(
            const mesh& m,
            unsigned cell_id,
            cell_type_param_ptr cell_type_parameters = nullptr
        ) noexcept;


        //Set the position of the cell in the cell_lst vector
        void set_local_id(const unsigned local_id) noexcept {local_id_ = local_id;}
        unsigned get_local_id() const noexcept {return local_id_;}


        //Returns the number of nodes used by the cell mesh
        size_t get_nb_of_nodes() const noexcept {return node_lst_.size() - free_node_queue_.size();}
        size_t get_nb_of_faces() const noexcept {return face_lst_.size() - free_face_queue_.size();}

        //Remove all the nodes and faces of the cell
        void clear_data() noexcept;

        //Compute the pressure difference between the cell and the environment
        void update_pressure() noexcept;

        //Get the face made up of the 3 nodes. If the face doesn't exist return a nullopt
        face* get_face(const unsigned n1_id, const unsigned n2_id, const unsigned n3_id) noexcept;

        //Load all the cell data e.g. cell edges, face areas, cell volume etc
        void initialize_cell_properties(bool check_cell_integrity = true) noexcept(false);

        //Initialize the cell growth rate and division volume
        void initialize_random_properties() noexcept;

        //Update the different cell properties and then apply the inner forces
        virtual void apply_internal_forces(const double time_step) noexcept;

        //Make sure that all the nodes are used by the faces and that the cell is not storing unsued nodes
        void remove_unused_nodes() noexcept;

        //Returns wether or not the cell is below the minimum allowed cell volume 
        bool is_below_min_vol() const noexcept{assert(cell_type_ != nullptr); return volume_ < cell_type_->min_vol_;}
        bool is_static() const noexcept {return is_static_;}

        //Apply a force on all the nodes of the cell
        void add_force(const vec3& force) noexcept{
            std::for_each(node_lst_.begin(), node_lst_.end(), [&force](node& n){
                if(n.is_used()) n.add_force(force);
            });
        }

        //A static method to check the winding order of the nodes of a given face is correct compare to an adjacent 
        //reference face. If the winding order is not correct the nodes of the face are swapped
        void check_face_winding_order(const face& ref_face, face& f) noexcept;

        void set_id(const unsigned cell_id) noexcept {cell_id_ = cell_id;}
        unsigned get_id() const noexcept {return cell_id_;};

        //Return the type of the cell
        cell_type_param_ptr get_cell_type() const noexcept {return cell_type_;}
        unsigned short get_cell_type_id() const noexcept {assert(cell_type_ != nullptr); return cell_type_->global_type_id_;}

        //Generate the set of edges of the cell from the faces
        void generate_edge_set() noexcept(false);

        //Returns a copy of the edge if it exists, otherwise returns an empty optional
        std::optional<edge> get_edge(const unsigned n1_id, const unsigned n2_id) const noexcept;

        const node& get_const_ref_node(const unsigned node_id) const noexcept;
        const face& get_const_ref_face(const unsigned face_id) const noexcept;

        //Get the mass of the cell
        double get_mass() const noexcept {
            assert(cell_type_ != nullptr);
            return cell_type_->mass_density_ * volume_;
        }

        //Return the mass of a node in this cell
        double get_node_mass() const noexcept{
            assert(get_nb_of_nodes() > 0);
            return get_mass() / static_cast<double>(get_nb_of_nodes());
        }

        //Get the type of a given face of the cell
        const face_type_parameters& get_face_type(const unsigned face_id) const noexcept{
            assert(face_id < face_lst_.size());
            assert(cell_type_ != nullptr);
            assert(face_lst_[face_id].type_id_ < cell_type_->face_types_.size());
            return cell_type_->face_types_[face_lst_[face_id].type_id_];
        }


        //Used after the cell division to maintain the cell inner pressure
        void set_target_volume(const double target_volume) noexcept {target_volume_ = target_volume;}

        //Remove a face from the cell mesh, the not is not deleted per see but marked as unused
        void delete_face(const unsigned face_id) noexcept(false);
        void delete_face(face& f) noexcept(false);

        //Remove a node from the cell mesh, the not is not deleted per see but marked as unused
        void delete_node(const unsigned node_id) noexcept;
        void delete_node(node& n) noexcept;

        //Add a new node to the cell mesh and return the node id
        unsigned create_node(const vec3& pos) noexcept;
        unsigned create_node(const double x, const double y, const double z) noexcept;
        unsigned add_node(const node& n) noexcept;

        //Replace the old node with the new node in all its faces and edges, return a vector of the edges that have been deleted
        std::pair<std::vector<edge>, std::vector<edge>> replace_node(const edge& start_edge, const unsigned old_node_id, const unsigned new_node_id) noexcept;

        //Add a new face to the cell mesh, and return the face id
        unsigned add_face(const face& f) noexcept(false);

        //Create a new face from the node ids, and return the face id
        unsigned create_face(
            const unsigned n1_id, 
            const unsigned n2_id, 
            const unsigned n3_id  
        ) noexcept(false);

        //Get a constant reference to the edge_set
        const edge_set& get_edge_set() const noexcept {return edge_set_;}

        //Return a copy of the node_lst or face_lst
        const std::vector<node>& get_node_lst() const noexcept {return node_lst_;};
        const std::vector<face>& get_face_lst() const noexcept {return face_lst_;};

        //Remove all the unused slots of the cell in its node_lst and face_lst
        void rebase() noexcept(false);


        //Returns a mesh object of the cell
        mesh get_mesh() const noexcept;

        //Get the ids of the nodes composing the faces of the cell
        std::vector<std::vector<unsigned>> get_face_point_ids() const noexcept{
            std::vector<std::vector<unsigned>> face_point_ids;
            std::for_each(face_lst_.begin(), face_lst_.end(), [&face_point_ids](const face& f){
                auto [n1, n2, n3] = f.get_node_ids();
                face_point_ids.push_back({n1, n2, n3});
            });
            return face_point_ids;
        }



        //Return in a flattened 1D vector the coordinates of the nodes
        std::vector<double> get_flat_node_coord_lst() const noexcept;

        //Return the coordinates of the nodes in a 2D vector
        std::vector<std::vector<double>> get_node_coord_lst() const noexcept{
            std::vector<std::vector<double>> node_coord_lst;
            std::for_each(node_lst_.begin(), node_lst_.end(), [&node_coord_lst](const node& n){
                const vec3& node_pos = n.pos();
                node_coord_lst.push_back({node_pos.dx(), node_pos.dy(), node_pos.dz()});
            });
            return node_coord_lst;
        }

        //Orient all the face normals such that they point in the outward direction
        void check_face_normal_orientation() noexcept;

        //Compute the normal of the given face. WARNING: the normal is not a unit vector!!
        vec3 get_face_normal(const unsigned face_id) const noexcept;

        //Given 3 nodes, that form a face with edge vectors a = j - i, b = k - i, 
        //return the gradient of the angle between a and b with respect to the coordinates of the 3 nodes
        static std::array<vec3, 3> get_angle_gradient(const vec3& i, const vec3& j, const vec3& k) noexcept;

        //Apply a force to make sure that the angle of the faces are all aproximately 60 degrees
        void regularize_all_face_angles() noexcept;
        void regularize_face_angles(const face& f) noexcept;

        //Update the area and normal of each face
        void update_all_face_normals_and_areas() noexcept;

        //Update at the ssame time the normal and the area of the given face
        void update_face_normal_and_area(const unsigned face_id) noexcept;
        void update_face_normal_and_area(face& f) noexcept;

        //Update the centroid of the cell
        void update_centroid() noexcept{centroid_ = compute_centroid();}

        //Return the volume of the cell
        double compute_volume() const noexcept;

        //Return the area of the cell const noexcept;
        double compute_area() const noexcept;

        //Compute the centroid of the cell 
        vec3 compute_centroid() const noexcept;

        //Check if the cell has a hole in its surface
        bool is_manifold() const noexcept;

        //Indicate on which side of a face a point lies
        bool point_is_on_positive_side_of_face(const unsigned face_id, const vec3& point) const noexcept;

        //Get the axis aligned bounding box of the cell
        std::array<double, 6> get_aabb() const noexcept;

        //Return a cell of the same type as the current cell (epithelial, ecm, static, lumen, nucleus)
        //If called on the base class you'll get an error message
        virtual cell_ptr get_cell_same_type(const mesh& m) noexcept(false){
            throw std::runtime_error("Cannot call get_cell_same_type on the base class cell");
        };


        //Indicate to all the faces of the cell that they belong to this cell
        void set_face_owner_cell() noexcept{
            std::for_each(face_lst_.begin(), face_lst_.end(), [this](auto& f){f.owner_cell_ = shared_from_this();});
        }

        //Set the local id of the faces and the nodes
        void set_local_ids() noexcept{
            for(unsigned i = 0; i < face_lst_.size(); ++i){face_lst_[i].local_face_id_ = i;}
            for(unsigned i = 0; i < node_lst_.size(); ++i){node_lst_[i].node_id_ = i;}
        }

        //Returns the ids of all the nodes that are connected to the given node by an edge
        std::vector<unsigned> get_connected_nodes(const unsigned node_id, const edge& e_start) const noexcept;

        virtual bool is_ready_to_divide() const noexcept {return false;}

        //By default return the cell longest axis to divide the cells
        virtual vec3 get_cell_division_axis() const noexcept {return get_cell_longest_axis();}

        //Return the cell longest axis
        vec3 get_cell_longest_axis() const noexcept;

        //Indicate to the cell that one of its face is in contact with another face
        virtual void face_is_in_contact(const unsigned face_local_id, const cell_ptr other_cell, const face& other_face) noexcept {}
        virtual void face_is_in_contact(const unsigned face_local_id, const cell_ptr other_cell) noexcept {}

    

        //Ask to the cell to update the type of its faces
        virtual void update_face_types() noexcept {}
  
        //Return the fraction of the cell surface that is in contact with other cells
        virtual double get_contact_area_fraction() {return 0.;}

        //Bunch of getters
        const vec3& get_centroid() const noexcept {return centroid_;}
        double get_division_volume() const noexcept {return division_volume_;}
        double get_growth_rate() const noexcept {return growth_rate_;}
        void set_growth_rate(const double growth_rate) noexcept {growth_rate_ = growth_rate;}
        double get_area() const noexcept {return area_;}
        double get_volume() const noexcept {return volume_;}
        double get_pressure() const noexcept {return pressure_;}
        double get_target_volume() const noexcept {return target_volume_;}
        float get_surface_tension_energy() const noexcept {return surface_tension_energy_;}
        float get_membrane_elasticity_energy() const noexcept {return membrane_elasticity_energy_;}
        float get_bending_energy() const noexcept {return bending_energy_;}
        float get_pressure_energy() const noexcept {return pressure_energy_;}
        float get_kinetic_energy() const noexcept {return kinetic_energy_;}


        virtual void special_polarization_update(const std::vector<cell_ptr>& cell_lst) noexcept{};


        #if CONTACT_MODEL_INDEX == 1  || CONTACT_MODEL_INDEX == 2
            //Use the cotan-laplacian to get a sense of the curvature at the nodes of the cell
            void compute_node_curvature_and_normals() noexcept;
        #endif


        #if FACE_STORE_CONTACT_ENERGY
            float get_adhesion_energy() const noexcept {
                return std::accumulate(face_lst_.begin(), face_lst_.end(), 0., [](auto acc, const auto& f){
                    return acc + f.get_adhesion_energy();
                });
            }

            float get_repulsion_energy() const noexcept {
                return std::accumulate(face_lst_.begin(), face_lst_.end(), 0., [](auto acc, const auto& f){
                    return acc + f.get_repulsion_energy();
                });
            } 
        #endif




};

#endif