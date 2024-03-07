#ifndef DEF_MESH_DATA
#define DEF_MESH_DATA

#include "utils.hpp"

#include "global_configuration.hpp"
#include "mesh_data_mappers.hpp"



/*
This file contains the functions that write the data of the cells in the output meshes. 
*/




//-------------------------------------------------------------------------------------------------------
//Store in this vector some mappers to map the cell area, pressure volume to the output VTK mesh
inline std::vector<cell_data_mapper> cell_data_mapper_lst
{

    //Print the ids of the cells
    cell_data_mapper("cell_id", "int", [](cell_ptr c) -> std::string {
        return format_number(c->get_id(), "%d");
    }),

    //Print the id of the cell type
    cell_data_mapper("cell_type_id", "int", [](cell_ptr c) -> std::string {

        //If the cell does not belong to a cell type, then print -1
        int cell_type_id = (c->get_cell_type() != nullptr) ? c->get_cell_type()->global_type_id_ : -1;
        return format_number(cell_type_id, "%d");
    }),


    //Print the areas of the cells
    cell_data_mapper("cell_area", "float", [](cell_ptr c) -> std::string {
        return format_number(c->get_area(), "%.3e");
    }),

    //Print the volumes of the cells
    cell_data_mapper("cell_volume", "float", [](cell_ptr c) -> std::string {
        return format_number(c->get_volume(), "%.3e");
    }),

    //Print the pressures of the cells
    cell_data_mapper("cell_pressure", "float", [](cell_ptr c) -> std::string {
        return format_number(c->get_pressure(), "%.3e");
    }),

    //Print the cell growth rates
    cell_data_mapper("cell_growth_rate", "float", [](cell_ptr c) -> std::string {
        return format_number(c->get_growth_rate(), "%.3e");
    }),

    //Print the cell division volumes
    cell_data_mapper("cell_division_vol", "float", [](cell_ptr c) -> std::string {
        return std::isfinite(c->get_division_volume()) ? format_number(c->get_division_volume(), "%.3e") : "-1.0";
    }),

    //Print wether or not the cell is static
    cell_data_mapper("cell_is_static", "int", [](cell_ptr c) -> std::string {
        return (c->is_static()) ? "1" : "0";
    })
};
//-------------------------------------------------------------------------------------------------------



//-------------------------------------------------------------------------------------------------------
//All the functions in this vector will be called to write the data of the cells in an output csv file
inline std::vector<cell_data_mapper> file_data_mapper_lst
{

    //Print the ids of the cells
    cell_data_mapper("cell_id", "int", [](cell_ptr c) -> std::string {
        return format_number(c->get_id(), "%d");
    }),

    //Print the id of the cell type
    cell_data_mapper("type_id", "int", [](cell_ptr c) -> std::string {

        //If the cell does not belong to a cell type, then print -1
        int cell_type_id = (c->get_cell_type() != nullptr) ? c->get_cell_type()->global_type_id_ : -1;
        return format_number(cell_type_id, "%d");
    }),


    //Print the areas of the cells
    cell_data_mapper("area", "float", [](cell_ptr c) -> std::string {
        return format_number(c->get_area(), "%.3e");
    }),

    cell_data_mapper("target_area", "float", [](cell_ptr c) -> std::string {

        //Get the type of the cell
        const auto cell_type = c->get_cell_type();
        if(cell_type == nullptr) return format_number(0.0, "%.3e");

        //Get the volume of the cell
        const double cell_volume = c->get_volume();
        return format_number(std::cbrt(cell_type->target_isoperimetric_ratio_ * cell_volume * cell_volume), "%.3e");
    }),

    //Print the volumes of the cells
    cell_data_mapper("volume", "float", [](cell_ptr c) -> std::string {
        return format_number(c->get_volume(), "%.3e");
    }),

    //Print the target volumes of the cells
    cell_data_mapper("target_volume", "float", [](cell_ptr c) -> std::string {
        return format_number(c->get_target_volume(), "%.3e");
    }),

    //Print the pressures of the cells
    cell_data_mapper("pressure", "float", [](cell_ptr c) -> std::string {
        return format_number(c->get_pressure(), "%.3e");
    }),

    //Save the fraction of the cell surface which is in contact with other cells. 
    //Only works if the cell is epithelial and POLARIZATION_MODE_INDEX == 1 = true
    #if POLARIZATION_MODE_INDEX == 1
        cell_data_mapper("cell_contact_area_fraction", "float", [](cell_ptr c) -> std::string {
            return format_number(c->get_contact_area_fraction(), "%.3f");
        }),
    #endif

    #if FACE_STORE_CONTACT_ENERGY
        //Print the adhesion energy of the faces
        cell_data_mapper("adhesion_energy", "float", [](cell_ptr c) -> std::string {
            return format_number(c->get_adhesion_energy(), "%.3e");
        }),

        cell_data_mapper("repulsion_energy", "float", [](cell_ptr c) -> std::string {
            return format_number(c->get_repulsion_energy(), "%.3e");
        }),
    #endif


    cell_data_mapper("kinetic_energy", "float", [](cell_ptr c) -> std::string {
        return format_number(c->get_kinetic_energy(), "%.3e");
    }),

    cell_data_mapper("surface_tension_energy", "float", [](cell_ptr c) -> std::string {
        return format_number(c->get_surface_tension_energy(), "%.3e");
    }),

    cell_data_mapper("membrane_elasticity_energy", "float", [](cell_ptr c) -> std::string {
        return format_number(c->get_membrane_elasticity_energy(), "%.3e");
    }),

    cell_data_mapper("bending_energy", "float", [](cell_ptr c) -> std::string {
        return format_number(c->get_bending_energy(), "%.3e");
    }),


    cell_data_mapper("pressure_energy", "float", [](cell_ptr c) -> std::string {
        return format_number(c->get_pressure_energy(), "%.3e");
    }),


    cell_data_mapper("total_potential_energy", "float", [](cell_ptr c) -> std::string {
        const float total_potential_energy = c->get_surface_tension_energy() + c->get_membrane_elasticity_energy() + c->get_bending_energy() + c->get_pressure_energy();
        return format_number(total_potential_energy, "%.3e");
    }),

};
//---------------------------------------------------------------------------------------------------------




//-------------------------------------------------------------------------------------------------------
//Store in this vector some mappers to map the face attributes to the output VTK mesh
inline std::vector<face_data_mapper> face_data_mapper_lst
{

    //Print the areas of the faces
    face_data_mapper("face_area", "float", [](const face& f) -> std::string {
        return format_number(f.get_area(), "%.3e");
    }),

    //Print the id of the cell owning this face
    face_data_mapper("face_cell_id", "float", [](const face& f) -> std::string {
        int cell_id = f.get_owner_cell() == nullptr ? -1 : f.get_owner_cell()->get_id();
        return format_number(cell_id, "%d");
    }),

    //Print the id of the face type
    face_data_mapper("face_type_id", "int", [](const face& f) -> std::string {
        cell_ptr owner_cell = f.get_owner_cell();

        //If the face does not belonng to any cell, the face type id is -1
        const int face_type_id = (owner_cell != nullptr) ? owner_cell->get_face_type(f.get_local_id()).face_type_global_id_ : -1;
        return format_number(face_type_id, "%d");
    }),


    //Print the id of the face type
    face_data_mapper("face_surface_tension", "float", [](const face& f) -> std::string {
        cell_ptr owner_cell = f.get_owner_cell();
        const double surface_tension = (owner_cell != nullptr) ? owner_cell->get_face_type(f.get_local_id()).surface_tension_ : -1.;
        return format_number(surface_tension, "%.3e");
    }),


    #if FACE_STORE_CONTACT_ENERGY
        //Print the adhesion energy of the faces
        face_data_mapper("face_adhesion_energy", "float", [](const face& f) -> std::string {
            return format_number(f.get_adhesion_energy(), "%.3e");
        }),

        //Print the repulsion energy of the faces
        face_data_mapper("face_repulsion_energy", "float", [](const face& f) -> std::string {
            return format_number(f.get_repulsion_energy(), "%.3e");
        }),
    #endif
};
//-------------------------------------------------------------------------------------------------------




//-------------------------------------------------------------------------------------------------------
//Store in this vector some mappers to map the node attributes to the output VTK mesh
inline std::vector<node_data_mapper> node_data_mapper_lst
{

    #if CONTACT_MODEL_INDEX == 2


    //Print wether or not a node is coupled to a face
    node_data_mapper("node_is_coupled", "int", [](const node& n) -> std::string {
        return (n.is_coupled()) ? "1" : "0";
    }),


    node_data_mapper("mean_curvature", "double", [](const node& n) -> std::string {
        return format_number(n.get_curvature(), "%.3e");
    }),

    
    //Print wether or not a node is coupled to a face
    node_data_mapper("nb_coupled_node", "int", [](const node& n) -> std::string {
        return format_number(n.get_nb_coupled_nodes(), "%d");
    })

    #endif
};
//-------------------------------------------------------------------------------------------------------




#endif