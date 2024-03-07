#ifndef DEF_SOLVER
#define DEF_SOLVER


#include <string>
#include <cassert>
#include <vector>
#include <filesystem>
#include <omp.h>

/*
    This class contains the solver of the simulation. First the solver initialize all the modules used to run 
    the simulation (local_mesh_refiner, cell_contacts, time_integration etc). Then it runs the main loop of the 
    program.
*/


#include "custom_structures.hpp"
#include "custom_exception.hpp"
#include "local_mesh_refiner.hpp"
#include "time_integration.hpp"
#include "mesh_writer.hpp"
#include "cell_divider.hpp"

#include "automatic_polarizer.hpp"
#include "automatic_polarization_writer.hpp"

#include "contact_model_abstract.hpp"
#include "statistics_writer.hpp"
#include "cell.hpp"


//Load the correct contact model depending on the value of the macro CONTACT_MODEL_INDEX
#include "contact_node_node_via_coupling.hpp"
#include "contact_node_face_via_spring.hpp"
#include "contact_face_face_via_coupling.hpp"


//--------------------------------------------------------------------------------------------------------------------------
class solver{
    
    protected:
        //The parameters used to run the simulation, they come from the input file
        //xml parameter file.
        global_simulation_parameters sim_parameters_;

        //The local mesh refiner is used to make sure that all the edges of the cell meshes have lengths
        std::unique_ptr<local_mesh_refiner> lmr_ptr_;

        //There are two types of statistic writers: one that writes the statistics in a file and one that writes in a string
        std::unique_ptr<abstract_statistics_writer> statistic_writer_ptr_;

        //There are also 2 types of contact models: one for face-face contacts and one for node-face contacts
        std::unique_ptr<contact_model_abstract> contact_model_ptr_;

        //There are also 2 types of time integrators: one for forward euler and one for semi-implicit euler
        std::unique_ptr<time_integration_scheme> time_integrator_ptr_;

        #if POLARIZATION_MODE_INDEX == 2
            std::unique_ptr<automatic_polarizer> cell_surface_polarizer_ptr_;
        #endif

        //The list of cells that will be simulated. Their geometries come from the 
        //input mesh file.
        std::vector<cell_ptr> cell_lst_;

        //The number of iterations of the main loop
        unsigned int iteration_ = 0;
        unsigned int file_number_ = 0;

        //Use this number to assign a unique id to each cell
        unsigned max_cell_id_ = 0;

        //The time of the simulation
        double simulation_time_ = 0.;

        //Number of threads used to run the computation in parallel
        unsigned short int nb_threads_;

        //If true, the program will print some information about the simulation progress
        bool verbose_;


    public:
        solver() = default;                          //default constructor
        solver(const solver& v) = delete;            //copy constructor
        solver(solver&& v) = delete;                 //move constructor
        solver& operator=(const solver& v) = delete; //copy assignment operator
        solver& operator=(solver&& v) = default;     //move assignment operator 
      
        //Constructor
        solver(
            const global_simulation_parameters& sim_parameters, 
            const std::vector<cell_ptr>& cell_lst,
            int nb_threads = -1, //If set to -1, SimuCell3D will use the available cores
            bool write_cell_stats_in_string = false,  //If true, the program will write the cell statistics in a string instead of a file   
            bool verbose = true //If true, the program will print some information about the simulation progress
        ) noexcept(false);

        //The method that contains the main loop of the program
        void run() noexcept(false);

        //Run one iteration of the program. This method can be overloaded and thus allow
        //to easily create new derived solvers
        virtual void run_iteration() noexcept(false);

        //Write the output file
        void save_mesh() noexcept(false);

        const std::vector<cell_ptr>& get_cell_lst() const noexcept{return cell_lst_;}
        const global_simulation_parameters& get_sim_parameters() const noexcept{return sim_parameters_;}


        //If the simulation statistics have been written in a string, then 
        //returns the string. Do not call this function if you have 
        //constructed the solver with write_cell_stats_in_string_ = false
        std::string get_simulation_statistics() const noexcept(false);

};
//--------------------------------------------------------------------------------------------------------------------------


#endif