#include "simucell3d_wrapper.hpp"

namespace py = pybind11;




/*
    This class is a python wrapper. It can be called from python and it will run the simulation.
    Then the results of the simulation can be extracted from this class and analyzed in python.
    The initialization and the execution of the simulation are both done in the constructor.
*/


//--------------------------------------------------------------------------------------------------------------------------
//Custom constructor, all the code is run from here
simucell3d_wrapper::simucell3d_wrapper(
            const global_simulation_parameters& sim_parameters,
            const std::vector<cell_type_parameters>& cell_type_param_lst,
            const bool verbose, //If set to true, the progress of the simulation will be printed
            const int nb_threads
        ) noexcept(false){

    //Transform the cell_type_param to shared pointers
    std::vector<std::shared_ptr<cell_type_parameters>> cell_type_param_lst_shared_ptr;
    cell_type_param_lst_shared_ptr.reserve(cell_type_param_lst.size());
    for(const auto& cell_type_param : cell_type_param_lst){
        cell_type_param_lst_shared_ptr.push_back(std::make_shared<cell_type_parameters>(cell_type_param));
    }

    //Catch any exception and return error msg
    try{
        //Load all the parameters and geometrical information of the tissue
        simulation_initializer sim_init(sim_parameters,  cell_type_param_lst_shared_ptr, verbose);

        //Initialize the solver, we use a derived class of the solver that can be stopped with a ctrl-c command
        solver_python_wrapper solver_(
            sim_init.get_simulation_parameters(), 
            sim_init.get_cell_lst(), 
            nb_threads,
            true, //Write the cell staistics to a string
            verbose
        );

        //Run the simulation
        solver_.run();

        //If the simulation was successful, return the statistics collected about the cells
        simulation_outputs_.RETURN_CODE_ = 0;
        simulation_outputs_.cell_statistics_ = solver_.get_simulation_statistics();
        simulation_outputs_.cell_lst_ = solver_.get_cell_lst();
    }
    
    //Catch and print any exception
    catch(std::exception const& e){
        //If the simulation failed
        simulation_outputs_.RETURN_CODE_ = 1;
        simulation_outputs_.error_message_ = e.what();
    }
}






//----------------------------------------------------------------------------------------
PYBIND11_MODULE(simucell3d_python_wrapper, m) {

    //The stucture used to store the parameters of the simulation
    py::class_<global_simulation_parameters>(m, "global_simulation_parameters")
        .def(py::init<>())
        .def_readwrite("output_folder_path_",           &global_simulation_parameters::output_folder_path_)
        .def_readwrite("input_mesh_path_",              &global_simulation_parameters::input_mesh_path_)
        .def_readwrite("damping_coefficient_",          &global_simulation_parameters::damping_coefficient_)
        .def_readwrite("simulation_duration_",          &global_simulation_parameters::simulation_duration_)
        .def_readwrite("sampling_period_",              &global_simulation_parameters::sampling_period_)
        .def_readwrite("perform_initial_triangulation_",            &global_simulation_parameters::perform_initial_triangulation_)
        .def_readwrite("time_step_",                    &global_simulation_parameters::time_step_)
        .def_readwrite("min_edge_len_",                 &global_simulation_parameters::min_edge_len_)
        .def_readwrite("contact_cutoff_adhesion_",      &global_simulation_parameters::contact_cutoff_adhesion_)
        .def_readwrite("contact_cutoff_repulsion_",     &global_simulation_parameters::contact_cutoff_repulsion_)
        .def_readwrite("enable_edge_swap_operation_",   &global_simulation_parameters::enable_edge_swap_operation_);


    //The structure used to store the parameters of the cells
    py::class_<cell_type_parameters>(m, "cell_type_parameters")
        .def(py::init<>())
        .def_readwrite("name_",                             &cell_type_parameters::name_)
        .def_readwrite("global_type_id_",                   &cell_type_parameters::global_type_id_)
        .def_readwrite("mass_density_",                     &cell_type_parameters::mass_density_)
        .def_readwrite("bulk_modulus_",                     &cell_type_parameters::bulk_modulus_)
        .def_readwrite("max_pressure_",                     &cell_type_parameters::max_pressure_)
        .def_readwrite("initial_pressure_",                 &cell_type_parameters::initial_pressure_)
        .def_readwrite("area_elasticity_modulus_",          &cell_type_parameters::area_elasticity_modulus_)
        .def_readwrite("avg_division_vol_",                 &cell_type_parameters::avg_division_vol_)
        .def_readwrite("std_division_vol_",                 &cell_type_parameters::std_division_vol_)
        .def_readwrite("avg_growth_rate_",                  &cell_type_parameters::avg_growth_rate_)
        .def_readwrite("std_growth_rate_",                  &cell_type_parameters::std_growth_rate_)
        .def_readwrite("angle_regularization_factor_",      &cell_type_parameters::angle_regularization_factor_)
        .def_readwrite("target_isoperimetric_ratio_",       &cell_type_parameters::target_isoperimetric_ratio_)
        .def_readwrite("min_vol_",                          &cell_type_parameters::min_vol_)
        .def_readwrite("surface_coupling_max_curvature_",   &cell_type_parameters::surface_coupling_max_curvature_)
        .def("add_additional_parameter",                    &cell_type_parameters::add_additional_parameter)
        .def("add_face_type",                               &cell_type_parameters::add_face_type);


    //The structure used to store the parameters of the faces
    py::class_<face_type_parameters>(m, "face_type_parameters")
        .def(py::init<>())
        .def_readwrite("name_",                         &face_type_parameters::name_)
        .def_readwrite("face_type_global_id_",          &face_type_parameters::face_type_global_id_)
        .def_readwrite("surface_tension_",              &face_type_parameters::surface_tension_)
        .def_readwrite("adherence_strength_",           &face_type_parameters::adherence_strength_)
        .def_readwrite("repulsion_strength_",           &face_type_parameters::repulsion_strength_)
        .def_readwrite("bending_modulus_",              &face_type_parameters::bending_modulus_);

    //The structure used to store the outputs of the simulation
    py::class_<simulation_outputs>(m, "simulation_outputs")
        .def(py::init<>())
        .def_readwrite("RETURN_CODE_",                  &simulation_outputs::RETURN_CODE_)
        .def_readwrite("error_message_",                &simulation_outputs::error_message_)
        .def_readwrite("cell_statistics_",              &simulation_outputs::cell_statistics_)
        .def_readwrite("cell_lst_",                     &simulation_outputs::cell_lst_);

    //This class makes the cell object accessible from python
    py::class_<cell, std::shared_ptr<cell>>(m, "cell")
        .def("get_id",              &cell::get_id)
        .def("get_cell_type_id",    &cell::get_cell_type_id)
        .def("get_growth_rate",     &cell::get_growth_rate)
        .def("get_area",            &cell::get_area)
        .def("get_volume",          &cell::get_volume)
        .def("get_pressure",        &cell::get_pressure)
        .def("get_target_volume",   &cell::get_target_volume)
        .def("get_node_coord_lst",  &cell::get_node_coord_lst)
        .def("get_face_point_ids",  &cell::get_face_point_ids);

    //The class used to make the c++ code accessible from python
    py::class_<simucell3d_wrapper>(m, "simucell3d_wrapper")
        .def(py::init<const global_simulation_parameters&, const std::vector<cell_type_parameters>&, bool, int>())
        .def("get_simulation_outputs",                  &simucell3d_wrapper::get_simulation_outputs);

}



