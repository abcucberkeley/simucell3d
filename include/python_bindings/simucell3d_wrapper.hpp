#ifndef DEF_SIMUCELL3D_WRAPPER
#define DEF_SIMUCELL3D_WRAPPER


#include <string>
#include <cassert>
#include <vector>
#include <filesystem>
#include <omp.h>

#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

/*
    This class is a python wrapper. It can be called from python and it will run the simulation.
    Then the results of the simulation can be extracted from this class and analyzed in python.
*/


#include "custom_structures.hpp"
#include "custom_exception.hpp"
#include "solver.hpp"
#include "solver_python_wrapper.hpp"
#include "simulation_initializer.hpp"

class simucell3d_wrapper{
    private:
        //Return the state of the simulation in this structure, which can be used in python
        simulation_outputs simulation_outputs_;


    public:
        simucell3d_wrapper() = default;                          //default constructor
        simucell3d_wrapper(const simucell3d_wrapper& v) = delete;            //copy constructor
        simucell3d_wrapper(simucell3d_wrapper&& v) = delete;                 //move constructor
        simucell3d_wrapper& operator=(const simucell3d_wrapper& v) = delete; //copy assignment operator
        simucell3d_wrapper& operator=(simucell3d_wrapper&& v) = default;     //move assignment operator 

        //Custom constructor, all the code is run from here
        simucell3d_wrapper(
            const global_simulation_parameters& sim_parameters,
            const std::vector<cell_type_parameters>& cell_type_param_lst,
            const bool verbose, //If set to true, the progress of the simulation will be printed
            const int nb_threads //Number of threads to use in the simulatiions
        ) noexcept(false);

        //Return if the simulation was successful, and returns the statistics collected about the cells
        simulation_outputs get_simulation_outputs() const noexcept(false){return simulation_outputs_;}
};

#endif