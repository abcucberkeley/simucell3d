#include "test_cell_separation.hpp"



/*
Initialize the simulation with a cell doublet geometry and then pull on the cells in opposite directions
to separate them. This test is used to check that the contact model contact_node_node_via_coupling allows two adjacent 
cells to be separated.
*/


//---------------------------------------------------------------------------------------------------------
int run_separation_test() noexcept{

    //Make sure that the correct contact model is used when the test is run
    #if CONTACT_MODEL_INDEX != 1
        std::cout << "ERROR: the CONTACT_MODEL_INDEX macro must be set to 1 to run this test" << std::endl;
        std::cout << "You can modify it in the file include/global_configuration.hpp" << std::endl;
        return 1;
    #endif

    
    //Generate the parameters for the simulation
    //Numerical parameters
    global_simulation_parameters sim_parameters;
    sim_parameters.input_mesh_path_ =               std::string(PROJECT_SOURCE_DIR) + "/test/test_contact_models/test_contact_node_node_via_coupling/test_cell_separation/cell_doublet.vtk";
    sim_parameters.output_folder_path_ =            std::string(CMAKE_BINARY_DIR) + "/test_cell_separation/";
    sim_parameters.damping_coefficient_ =           2.0e-09;
    sim_parameters.simulation_duration_ =           1.0e-02; 
    sim_parameters.sampling_period_ =               2.0e-06;
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
    auto cell_type_ptr = std::make_shared<cell_type_parameters>();
    cell_type_ptr->name_                             = "epithelial";     
    cell_type_ptr->global_type_id_                   = 0; 
    cell_type_ptr->mass_density_                     = 1e3;
    cell_type_ptr->bulk_modulus_                     = 1e4;   
    cell_type_ptr->initial_pressure_                 = 0;
    cell_type_ptr->max_pressure_                     = 1e20;
    cell_type_ptr->avg_growth_rate_                  = 0;
    cell_type_ptr->std_growth_rate_                  = 0;  
    cell_type_ptr->target_isoperimetric_ratio_       = 150;
    cell_type_ptr->angle_regularization_factor_      = 0;
    cell_type_ptr->area_elasticity_modulus_          = 0;  
    cell_type_ptr->surface_coupling_max_curvature_   = 8e5;
    cell_type_ptr->avg_division_vol_                 = 1e20;  
    cell_type_ptr->std_division_vol_                 = 0;   
    cell_type_ptr->min_vol_                          = 5e-17;   
    cell_type_ptr->add_face_type(face_type_1);
    cell_type_ptr->add_face_type(face_type_2);


    try{
        //Call the simulation initializer
        simulation_initializer initializer(sim_parameters, {cell_type_ptr}, /*verbose = */ true);

        //Initialize the solver, we use a custom solver where we add a force in opposite direction to the 2 cells
        //at each iteration
        solver_separation_test solver(
            /*sim_parameters = */               initializer.get_simulation_parameters(), 
            /*cell_lst = */                     initializer.get_cell_lst(), 
            /*nb_threads = */                   1,
            /*write_cell_stats_in_string =*/    true, 
            /*verbose = */                      true
        );

        //Run the simulation
        solver.run();
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

    return run_separation_test();

}
//---------------------------------------------------------------------------------------------------------
