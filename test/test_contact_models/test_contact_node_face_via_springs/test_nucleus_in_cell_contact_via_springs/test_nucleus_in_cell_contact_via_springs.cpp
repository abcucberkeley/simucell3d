#include "test_nucleus_in_cell_contact_via_springs.hpp"


//Check the contact model contact_node_face_via_spring by growing a nucleus inside a cell. If the 
//the nucleus is not pushed out of the cell, then the contact model is working correctly


//---------------------------------------------------------------------------------------------------------
int run_nucleus_in_cell_test() noexcept{

    //Make sure that the correct contact model is used when the test is run
    #if CONTACT_MODEL_INDEX != 0
        std::cout << "ERROR: the CONTACT_MODEL_INDEX macro must be set to 0 to run this test" << std::endl;
        std::cout << "You can modify it in the file include/global_configuration.hpp" << std::endl;
        return 1;
    #endif

    
    //Generate the parameters for the simulation
    //Numerical parameters
    global_simulation_parameters sim_parameters;                                        
    sim_parameters.input_mesh_path_ =               std::string(PROJECT_SOURCE_DIR) + "/test/test_contact_models/test_contact_node_face_via_springs/test_nucleus_in_cell_contact_via_springs/nucleus_in_cell_geometry.vtk";
    sim_parameters.output_folder_path_ =            std::string(CMAKE_BINARY_DIR) + "/test_nucleus_in_cell_contact_via_springs/";
    sim_parameters.damping_coefficient_ =           2.0e-09;
    sim_parameters.simulation_duration_ =           1.0e-03; 
    sim_parameters.sampling_period_ =               5.0e-06;
    sim_parameters.time_step_ =                     1.0e-07;
    sim_parameters.min_edge_len_ =                  7.5e-07;
    sim_parameters.contact_cutoff_adhesion_ =       5.0e-07;
    sim_parameters.contact_cutoff_repulsion_=       5.0e-07;
    sim_parameters.enable_edge_swap_operation_ =    false;
    sim_parameters.perform_initial_triangulation_ = false;

    //The parameters of the apical faces
    face_type_parameters face_type_1;
    face_type_1.name_                = "apical";     
    face_type_1.face_type_global_id_ = 0;
    face_type_1.adherence_strength_  = 0;
    face_type_1.repulsion_strength_  = 2e9;
    face_type_1.surface_tension_     = 1e-3;
    face_type_1.bending_modulus_     = 0;

    //The parameters of the lateral faces
    face_type_parameters face_type_2;
    face_type_2.name_                = "lateral";     
    face_type_2.face_type_global_id_ = 1;
    face_type_2.adherence_strength_  = 0;
    face_type_2.repulsion_strength_  = 2e9;
    face_type_2.surface_tension_     = 5e-4; 
    face_type_2.bending_modulus_     = 0;

    //The parameters global to the epithelial cells
    auto epi_cell_type_ptr = std::make_shared<cell_type_parameters>();
    epi_cell_type_ptr->name_                             = "epithelial";     
    epi_cell_type_ptr->global_type_id_                   = 0; 
    epi_cell_type_ptr->mass_density_                     = 1e3;
    epi_cell_type_ptr->bulk_modulus_                     = 1e4;   
    epi_cell_type_ptr->initial_pressure_                 = 0;
    epi_cell_type_ptr->max_pressure_                     = 2.5e2;
    epi_cell_type_ptr->avg_growth_rate_                  = 0;
    epi_cell_type_ptr->std_growth_rate_                  = 0;  
    epi_cell_type_ptr->target_isoperimetric_ratio_       = 150;
    epi_cell_type_ptr->angle_regularization_factor_      = 0;
    epi_cell_type_ptr->area_elasticity_modulus_          = 0;  
    epi_cell_type_ptr->surface_coupling_max_curvature_   = 8e5;
    epi_cell_type_ptr->avg_division_vol_                 = 1e20;  
    epi_cell_type_ptr->std_division_vol_                 = 0;   
    epi_cell_type_ptr->min_vol_                          = 5e-17;   
    epi_cell_type_ptr->add_face_type(face_type_1);
    epi_cell_type_ptr->add_face_type(face_type_2);

    //Create a face type for the ECM
    face_type_parameters face_type_3;
    face_type_3.name_                = "nucleus_face";     
    face_type_3.face_type_global_id_ = 5;
    face_type_3.adherence_strength_  = 0;
    face_type_3.repulsion_strength_  = 2e9;
    face_type_3.surface_tension_     = 1e-3; 
    face_type_3.bending_modulus_     = 0;

    //The parameters for the ECM cell type
    auto nuc_cell_type_ptr = std::make_shared<cell_type_parameters>();
    nuc_cell_type_ptr->name_                             = "nucleus";     
    nuc_cell_type_ptr->global_type_id_                   = 3; 
    nuc_cell_type_ptr->mass_density_                     = 1e3;
    nuc_cell_type_ptr->bulk_modulus_                     = 1e4;   
    nuc_cell_type_ptr->initial_pressure_                 = 0;
    nuc_cell_type_ptr->max_pressure_                     = 2.5e2;
    nuc_cell_type_ptr->avg_growth_rate_                  = 1e-11;
    nuc_cell_type_ptr->std_growth_rate_                  = 0;  
    nuc_cell_type_ptr->target_isoperimetric_ratio_       = 150;
    nuc_cell_type_ptr->angle_regularization_factor_      = 0;
    nuc_cell_type_ptr->area_elasticity_modulus_          = 0;  
    nuc_cell_type_ptr->surface_coupling_max_curvature_   = 8e5;
    nuc_cell_type_ptr->avg_division_vol_                 = 1e20;  
    nuc_cell_type_ptr->std_division_vol_                 = 0;   
    nuc_cell_type_ptr->min_vol_                          = 5e-17;   
    nuc_cell_type_ptr->add_face_type(face_type_3);



    try{
        //Call the simulation initializer
        simulation_initializer initializer(sim_parameters, {epi_cell_type_ptr, nuc_cell_type_ptr, nuc_cell_type_ptr, nuc_cell_type_ptr , nuc_cell_type_ptr}, /*verbose = */ true);

        //Initialize the solver, we use a custom solver where we add a force in opposite direction to the 2 cells
        //at each iteration
        solver base_solver(
            /*sim_parameters = */               initializer.get_simulation_parameters(), 
            /*cell_lst = */                     initializer.get_cell_lst(), 
            /*nb_threads = */                   1,
            /*write_cell_stats_in_string =*/    true, 
            /*verbose = */                      true
        );

        //Run the simulation
        base_solver.run();
    }

    //Catch and print any exception
    catch(std::exception const& e){
        std::cout << e.what() <<std::endl;
        return 1;
    }

    return 0;
}
//---------------------------------------------------------------------------------------------------------






//---------------------------------------------------------------------------------------------------------
int main (int argc, char** argv){

    return run_nucleus_in_cell_test();

}
//---------------------------------------------------------------------------------------------------------
