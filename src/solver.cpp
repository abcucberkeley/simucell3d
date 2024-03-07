#include "solver.hpp"


//--------------------------------------------------------------------------------------------------------------------------
//Constructor
solver::solver(
    const global_simulation_parameters& sim_parameters, 
    const std::vector<cell_ptr>& cell_lst,
    int nb_threads, 
    bool write_cell_stats_in_string, 
    bool verbose 
) noexcept(false): sim_parameters_(sim_parameters), cell_lst_(cell_lst), verbose_(verbose){

    assert(cell_lst_.size() > 0);

    //Remove the output folder if it already exists
    std::filesystem::remove_all(sim_parameters_.output_folder_path_);

    //Create the folders where the data will be written
    bool t1 = std::filesystem::create_directories(sim_parameters_.output_folder_path_);
    bool t2 = std::filesystem::create_directories(sim_parameters_.output_folder_path_ + "/cell_data");
    bool t3 = std::filesystem::create_directories(sim_parameters_.output_folder_path_ + "/face_data");
    if(!t2) throw intialization_exception("The output folder could not be created: " + sim_parameters_.output_folder_path_);
    if(!t2) throw intialization_exception("The output folder could not be created: " + sim_parameters_.output_folder_path_ + "/cell_data");
    if(!t3) throw intialization_exception("The output folder could not be created: " + sim_parameters_.output_folder_path_ + "/face_data");

    //Create a folder where the results of the automatic polarizer will be stored
    #if POLARIZATION_MODE_INDEX == 2
        bool t4 = std::filesystem::create_directories(sim_parameters_.output_folder_path_ + "/region_data");
        if(!t4) throw intialization_exception("The output folder could not be created: " + sim_parameters_.output_folder_path_ + "/region_data");
    #endif

    if(verbose_) std::cout << "The output folder is: " << std::filesystem::absolute(sim_parameters_.output_folder_path_)<< std::endl;

    //Set the ids of all the cells
    for(size_t i = 0; i< cell_lst_.size(); i++){
        cell_lst_[i]->set_id(max_cell_id_++);
        cell_lst_[i]->set_local_id(cell_lst_[i]->get_id());
    }

    //Initialize the local mesh refiner
    lmr_ptr_ = std::make_unique<local_mesh_refiner>(sim_parameters_.min_edge_len_, sim_parameters_.min_edge_len_ *  3., sim_parameters_.enable_edge_swap_operation_);

    //Initialize the time integration scheme
    time_integrator_ptr_ = std::make_unique<time_integration_scheme>(sim_parameters_, verbose);


    #if POLARIZATION_MODE_INDEX == 2
        cell_surface_polarizer_ptr_ = std::make_unique<automatic_polarizer>(sim_parameters_.min_edge_len_ *  3.);
    #endif

    //Initialize the contact model used to compute the contact forces between the cells
    #if CONTACT_MODEL_INDEX == 0
        contact_model_ptr_ = std::make_unique<contact_node_face_via_spring>(sim_parameters_);
    #elif CONTACT_MODEL_INDEX == 1
        contact_model_ptr_ = std::make_unique<contact_node_node_via_coupling>(sim_parameters_);
    #elif CONTACT_MODEL_INDEX == 2
        contact_model_ptr_ = std::make_unique<contact_face_face_via_coupling>(sim_parameters_);
    #else
        throw intialization_exception("The contact model index is not valid");
    #endif

    //Initialize the file writer. The simulation statistics can be written in a file or in a string
    if(write_cell_stats_in_string){
        statistic_writer_ptr_ = std::make_unique<string_statistics_writer>();
    }
    else{
        statistic_writer_ptr_ = std::make_unique<csv_file_statistics_writer>(sim_parameters_.output_folder_path_ + "/simulation_statistics.csv");
    }


    //Set the initial pressure of all the cells
    for(cell_ptr c: cell_lst_){

        //Get the type of the cell
        const cell_type_param_ptr cell_type_ = c->get_cell_type();

        //Compute the target volume of the cell based on its initial pressure
        const double target_volume_ = c->get_volume() * std::exp(cell_type_->initial_pressure_ / cell_type_->bulk_modulus_);

        //Set the target volume of the cell
        c->set_target_volume(target_volume_);
        c->update_pressure();
    }

    //By default the number of threads corresponds to all the available cores
    nb_threads_ = (nb_threads <= 0 ) ? omp_get_num_procs() : nb_threads;
    omp_set_num_threads(nb_threads_);
}
//--------------------------------------------------------------------------------------------------------------------------




//--------------------------------------------------------------------------------------------------------------------------
//The method that contains the main loop of the program
void solver::run() noexcept(false){

    //The main loop of the program
    while(time_integrator_ptr_->get_simulation_time() < sim_parameters_.simulation_duration_ && cell_lst_.size() > 0){
        run_iteration();
    }

    statistic_writer_ptr_->write_data(iteration_, time_integrator_ptr_->get_simulation_time(), cell_lst_); //noexcept(false)

    //Rebase all the cells before stopping the main loop
    std::for_each(cell_lst_.begin(), cell_lst_.end(), [](cell_ptr c){c->rebase();});
}
//--------------------------------------------------------------------------------------------------------------------------



//--------------------------------------------------------------------------------------------------------------------------
//Run one iteration of the program
void solver::run_iteration() noexcept(false){

    //Save the mesh if the time is right
    if(!time_integrator_ptr_->is_step_tmp()) save_mesh(); //noexcept(false)

    //Every 5 iterations check if some of the cells need to be divided
    if(!time_integrator_ptr_->is_step_tmp() && iteration_ % 5 == 0) cell_divider::run(cell_lst_, sim_parameters_.min_edge_len_, *lmr_ptr_, max_cell_id_, verbose_);

    #pragma omp parallel for
    for(size_t i = 0; i< cell_lst_.size(); i++){cell_lst_[i]->update_face_types();} 

    //Refine the mesh of the cells in parallel
    lmr_ptr_->refine_meshes(cell_lst_); //noexcept(false)

    //Compute the contact forces between the cells
    contact_model_ptr_->run(cell_lst_); //noexcept

    //If enabled, use the automatic polarizer to polarize the faces of the cells, 
    #if POLARIZATION_MODE_INDEX == 2
        cell_surface_polarizer_ptr_->polarize_faces(cell_lst_); //noexcept
    #endif

    //Polarize the faces based on their contacts
    #pragma omp parallel for
    for(size_t i = 0; i< cell_lst_.size(); i++){cell_lst_[i]->special_polarization_update(cell_lst_);}

    //Update the internal forces of the cells in parallel
    #pragma omp parallel for
    for(size_t i = 0; i< cell_lst_.size(); i++){cell_lst_[i]->apply_internal_forces(sim_parameters_.time_step_);}

    //Update the positions of the nodes in parallel
    time_integrator_ptr_->update_nodes_positions(cell_lst_); //noexcept

    //Save the cell properties every 50 iterations in a csv file
    if(iteration_ % 50 == 0) statistic_writer_ptr_->write_data(iteration_, time_integrator_ptr_->get_simulation_time(), cell_lst_); //noexcept(false)

    //Remove the cells that have a volume smaller than the minimum volume
    cell_lst_.erase(
        std::remove_if(cell_lst_.begin(), cell_lst_.end(), [](const cell_ptr& c){
            
            //Clear the data of the cell if it is below the minimum volume
            if(c->is_below_min_vol()) c->clear_data();

            //Return true if the cell is below the minimum volume
            return c->is_below_min_vol();

        }), 
        cell_lst_.end()
    );

    //Print the progression of the program
    if(verbose_ && !time_integrator_ptr_->is_step_tmp()) printf("Progression %d%%, iteration: %d, file number %d, nb cells %d\n", 
    static_cast<unsigned short>(100. * time_integrator_ptr_->get_simulation_time() / sim_parameters_.simulation_duration_),  iteration_, file_number_, static_cast<int>(cell_lst_.size()));

    //Update the iteration number
    iteration_++;
}     
//--------------------------------------------------------------------------------------------------------------------------








//--------------------------------------------------------------------------------------------------------------------------
//Write the mesh of the tissue in a VTK file
void solver::save_mesh() noexcept(false){

    unsigned new_file_nb = static_cast<unsigned>(std::floor(time_integrator_ptr_->get_simulation_time() / sim_parameters_.sampling_period_) + 1);
    if(new_file_nb != file_number_){
        file_number_ = new_file_nb;

        //Save the meshes of the cells in VTK format
        const std::string cell_mesh_path = sim_parameters_.output_folder_path_ + "/cell_data/result_" +std::to_string(file_number_)+ ".vtk";
        const std::string face_mesh_path = sim_parameters_.output_folder_path_ + "/face_data/result_" +std::to_string(file_number_)+ ".vtk";
        mesh_writer::write(cell_mesh_path, face_mesh_path, cell_lst_);
    
        //Save the discretization of the space by the automatic polarizer
        #if POLARIZATION_MODE_INDEX == 2
            if(iteration_ > 1) automatic_polarization_writer::write(sim_parameters_.output_folder_path_ + "/region_data/result_" +std::to_string(file_number_)+ ".vtk", *(cell_surface_polarizer_ptr_));
        #endif

    }
}
//--------------------------------------------------------------------------------------------------------------------------






//--------------------------------------------------------------------------------------------------------------------------
//If the simulation statistics have been written in a string, then 
//returns the string. Do not call this function if you have 
//constructed the solver with write_cell_stats_in_string = false
std::string solver::get_simulation_statistics() const noexcept(false){

    //Downcast the statistic_writer_ to a string_statistics_writer
    const string_statistics_writer* string_writer = dynamic_cast<const string_statistics_writer*>(statistic_writer_ptr_.get());
    if(string_writer == nullptr) throw std::runtime_error("The simulation statistics have not been written in a string");
    return string_writer->get_string();
}
//--------------------------------------------------------------------------------------------------------------------------
