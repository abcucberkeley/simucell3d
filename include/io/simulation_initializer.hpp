#ifndef DEF_SIMULATION_INITIALIZER
#define DEF_SIMULATION_INITIALIZER


#include <string>
#include <iostream>
#include <cassert>
#include <vector>


#include "custom_structures.hpp"
#include "custom_exception.hpp"
#include "parameter_reader.hpp"
#include "mesh_reader.hpp"
#include "mesh_writer.hpp"
#include "initial_triangulation.hpp"


#include "cell.hpp"
#include "epithelial_cell.hpp"
#include "static_cell.hpp"
#include "nucleus_cell.hpp"
#include "lumen_cell.hpp"
#include "ecm_cell.hpp"



/*
   This class contains all the methods that takes care of the initialization of the simulation.
   It first reads the parameters from the XML file, then it creates the cells and finally it triangulates the cell surfaces.
   All the simulation parameters are stored in structures that can then  be used by the solver class to run the simulation.
*/




class simulation_initializer{
   private:

      //The parameters that should be passed to the simulation_solver
      global_simulation_parameters sim_parameters_; 

      //The list of cells
      std::vector<cell_ptr> cell_lst_; 

      //If set to true all the cells will get the attribute of the 
      //epithelial cell type. However, the different cells can still possess 
      //different parameters
      bool make_all_cell_epithelial_ = false;

      //If set to false, the cells will not be triangulated
      bool perform_initial_triangulation_ = true;

      bool verbose_;
      

   public:
      simulation_initializer() = delete;                               //default constructor
      simulation_initializer(const simulation_initializer& v) = delete;           //copy constructor
      simulation_initializer(simulation_initializer&& v) = delete;                //move constructor
      simulation_initializer& operator=(const simulation_initializer& v) = delete;//copy assignment operator
      simulation_initializer& operator=(simulation_initializer&& v) = delete;     //move assignment operator 
      

      //Initialize the simulation by directly providing the structures containing the simulation and cell parameters
      simulation_initializer(
         const global_simulation_parameters& sim_parameters,
         const std::vector<cell_type_param_ptr>& cell_type_param_lst, 
         const bool verbose = true
      ) noexcept(false);

      //Initialize the simulation with an XML parameter file
      simulation_initializer(
         const std::string& parameter_file_path,
         const bool verbose = true
      ) noexcept(false);

      //Starts the initialization process
      void run(const std::vector<cell_type_param_ptr>& cell_type_param_lst) noexcept(false); 

      //Return the list of cells
      std::vector<cell_ptr> get_cell_lst() const noexcept {return cell_lst_;}

      //Return the simulation parameters
      global_simulation_parameters get_simulation_parameters() const noexcept {return sim_parameters_;}

      //Tries five times to triangulate the cell surface
      cell_ptr triangulate_surface(
        const mesh& cell_mesh, 
        const size_t cell_id,
        cell_type_param_ptr cell_type
      ) noexcept(false);

};

#endif