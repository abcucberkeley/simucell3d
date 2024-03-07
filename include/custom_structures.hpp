#ifndef DEF_CUSTOM_STRUCTURES
#define DEF_CUSTOM_STRUCTURES


#include <vector> 
#include <memory>
#include <map>
#include <string>
#include <functional>


#include "vec3.hpp"





//---------------------------------------------------------------------------------------
//Contains the parameters of the simulation (time step, damping coefficient, etc.)
struct global_simulation_parameters{

    //The path to the folder where all the simulation data will be written, it must be an absolute path
    std::string output_folder_path_;

    //The path to the input mesh file, it must be an absolute path
    std::string input_mesh_path_;

    //If equals 0 the cells will not be triangulated at the beginning of the simulation
    bool perform_initial_triangulation_ = 1;

    //If set to true, the local mesh refiner will remove the triangles with high aspect ratio
    //by swapping their longest edge
    bool enable_edge_swap_operation_ = true;

    //The damping coefficient
    double damping_coefficient_;

    //The simulation duration
    double simulation_duration_;

    //The time period between 2 consecutive mesh files written on the disk
    double sampling_period_;

    //The time step of the simulation
    double time_step_;

    //The minimum edge length allowed
    double min_edge_len_;

    //The distance between 2 cells above which cell contact forces are zero
    double contact_cutoff_adhesion_;
    double contact_cutoff_repulsion_;
};
//---------------------------------------------------------------------------------------



//---------------------------------------------------------------------------------------
//Store the parameters of the face types
struct face_type_parameters{
    
    //The name of the face_type
    std::string name_;

    //The id of the face type saved on the mesh, the cells internally use local ids to manage their
    //face types
    short face_type_global_id_ = 0;

    //The surface tension of the face
    double surface_tension_  = 0.;

    //The adherence strength between the faces of adjacent cells
    double adherence_strength_  = 0.;

    //The repulsion strength between the faces of adjacent cells
    double repulsion_strength_;

    //The bending rigidity of the face
    double bending_modulus_ = 0.;

    face_type_parameters() {}
};
//---------------------------------------------------------------------------------------


//---------------------------------------------------------------------------------------
//Store the parameters of the cell types
struct cell_type_parameters{

    //The name of the cell_type
    std::string name_;

    //The id of the cell type which is printed on the meshes
    short global_type_id_ = 0;

    //The mass density of the cell;
    double mass_density_ = 0.;

    //The bulk modulus of the cell cytoplasm (its compressibility)
    double bulk_modulus_ = 1.;

    //The maximal pressure difference between the cell and the exterior space
    double max_pressure_ = 0.;

    //The initial pressure of the cells of the given cell type
    double initial_pressure_ = 0.;

    //The elastic modulus of the cell membrane
    double area_elasticity_modulus_ = 0.;

    //The division volumes are drawn from a normal distribution with the following parameters
    double avg_division_vol_ = 0.;
    double std_division_vol_ = 0.;

    //The growth rates are drawn from a normal distribution with the following parameters
    double avg_growth_rate_ = 0.;
    double std_growth_rate_ = 0.;

    //The minial volume before cell removal
    double min_vol_ = 0.;

    //Used to make sure that the angles of the triangular faces are all approximately 60 deg
    double angle_regularization_factor_ = 0.;

    //The target isoperimetric ratio of the cells
    double target_isoperimetric_ratio_ = 0.;

    //The forces at which the links between the epithelial cells are broken
    double surface_coupling_max_curvature_ = 0.;

    //The different face types of the cell are stored in a vector
    std::vector<face_type_parameters> face_types_;

    //Additional parameters to the cell can be given via this map
    std::map<std::string, double> additional_parameters_;

    //Function used to add a face type to the cell type
    void add_face_type(const face_type_parameters& face_type) noexcept {face_types_.push_back(face_type);}

    //Function used to add an additional parameter to the cell type
    void add_additional_parameter(const std::string& parameter_name, const double parameter_value) noexcept {additional_parameters_[parameter_name] = parameter_value;}

    //Default constructor
    cell_type_parameters() {};
};
//---------------------------------------------------------------------------------------


//---------------------------------------------------------------------------------------
//Use this structure to store the different outputs that can be created during a simulation
//this is used by the python simucell3d_wrapper to pass the simulation outputs to the python interpreter
class cell; //Forward declaration of the cell class
struct simulation_outputs{
    
    //If the return code is 0, the simulation was successful, otherwise it failed
    bool RETURN_CODE_ = 0;

    //If an error was thrown, pass the error message to the python interpreter
    std::string error_message_ = "";

    //All the collected cell statistics will be stored in this string (it's basically the content of a csv file)
    std::string cell_statistics_ = "";

    //Contains pointer of the cells at the end of the simulation
    std::vector<std::shared_ptr<cell>> cell_lst_;

    simulation_outputs() {}
};
//---------------------------------------------------------------------------------------


//---------------------------------------------------------------------------------------
//Just a structure to store the points with their normals
struct oriented_point{
    unsigned id_;
    vec3 position_;
    vec3 normal_;

    //This is used during the mesh division to distinguish the points tha were created by the
    //poisson sampling from the points that were already in the mesh
    bool created_by_poisson_sampling_ = true; 

    //Constructors
    oriented_point(const vec3& point_pos, const vec3& point_normal) noexcept:
        position_(point_pos), normal_(point_normal){}

    oriented_point(const unsigned id, const vec3& point_pos, const vec3& point_normal) noexcept:
        id_(id), position_(point_pos), normal_(point_normal){}

    //Default constructor
    oriented_point():position_(vec3(0., 0., 0.)), normal_(vec3(0., 0., 0.)){}



};
//---------------------------------------------------------------------------------------


//---------------------------------------------------------------------------------------
//This is a very lightweight representation of a mesh. The  mesh stored in this structure does not 
//need to be triangulated. We use this structure to load the cell meshes or to divisde the cell meshes
//since it offers more freedom than the meshes stored in the cell objects.


//The data is stored in the following way in a mesh
//mesh.face_point_ids = {{f0_n0,f0_n1,f0_n2,f0_n3}, {f1_n0,f1_n1,f1_n2},...,{fn_n0,fn_n1,fn_n2}]
//mesh.node_pos_lst =   {x0,y0,z0, x1,y1,z1, ..., xn,yn,zn}   
//Where f0_n0 is id of the first node of the first face, f0_n1 is the id of the second node of the first face etc


struct mesh{

    //Store in a 2D vector the IDs of the nodes constituting the faces of the mesh
    std::vector<std::vector<unsigned>> face_point_ids;


    //Store the coordinates of the nodes of the mesh in this 1D vector
    std::vector<double> node_pos_lst;

    unsigned get_nb_nodes() const {return node_pos_lst.size() / 3;}
    unsigned get_nb_faces() const {return face_point_ids.size();}


    //Return the position of a node
    std::array<double, 3> get_node_pos(const unsigned node_id){
        assert(node_id < get_nb_nodes());

        return {
            node_pos_lst[node_id * 3    ],
            node_pos_lst[node_id * 3 + 1],
            node_pos_lst[node_id * 3 + 2]
        };
    } 
};
//---------------------------------------------------------------------------------------

#endif