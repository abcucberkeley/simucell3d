#ifndef DEF_TEST_CELL_SEPARATION
#define DEF_TEST_CELL_SEPARATION


/*
In this test, we check that the separation of two cells is possible when the contacts between the cells 
are resolved with the node-face coupling method.
*/
#include <string>
#include <cassert>
#include <vector>
#include <filesystem>
#include <omp.h>


#include "simulation_initializer.hpp"
#include "solver.hpp"
#include "custom_structures.hpp"


//--------------------------------------------------------------------------------------------------------------------------
class solver_separation_test: public solver{

    public:

        solver_separation_test() = default;                                          //default constructor
        solver_separation_test(const solver_separation_test& v) = delete;            //copy constructor
        solver_separation_test(solver_separation_test&& v) = delete;                 //move constructor
        solver_separation_test& operator=(const solver_separation_test& v) = delete; //copy assignment operator
        solver_separation_test& operator=(solver_separation_test&& v) = default;     //move assignment operator 
      
        //Constructor, just call the base class constructor
        solver_separation_test(
            const global_simulation_parameters& sim_parameters, 
            const std::vector<cell_ptr>& cell_lst,
            int nb_threads = -1, //If set to -1, SimuCell3D will use the available cores
            bool write_cell_stats_in_string = false,  //If true, the program will write the cell statistics in a string instead of a file   
            bool verbose = true //If true, the program will print some information about the simulation progress
        )  noexcept(false) : solver(sim_parameters, cell_lst, nb_threads, write_cell_stats_in_string, verbose){};

        
        //Add a force in opposite direction to the 2 cells
        void run_iteration() noexcept(false) override{

            //Make sure that there are 2 cells in the simulation
            if(cell_lst_.size() != 2) throw std::runtime_error("ERROR: the solver_separation_test class can only be used with 2 cells");

            //Get the 2 cells
            cell_ptr cell_1 = cell_lst_[0];
            cell_ptr cell_2 = cell_lst_[1];

            //Apply a force in opposite direction to the 2 cells
            const vec3 force(1e-10, 0.0, 0.0);

            cell_1->add_force(force * (-1.));
            cell_2->add_force(force);

            //Run the rest of the base class run_iteration method
            solver::run_iteration();
        };


};
//--------------------------------------------------------------------------------------------------------------------------




#endif









