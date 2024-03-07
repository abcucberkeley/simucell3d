#ifndef DEF_TIME_INTEGRATION
#define DEF_TIME_INTEGRATION


#include <cassert>

#include "global_configuration.hpp"

#include "node.hpp"
#include "cell.hpp"


/*
SimuCell3D possesses three methods to solve the equations of motions:
    - DYNAMIC_MODEL_INDEX = 0 : Full equations of motions solved with semi-implicit euler  (1st order integration scheme)
    - DYNAMIC_MODEL_INDEX = 1 : Overdamped equations of motions solved with forward euler  (1st order integration scheme)
    - DYNAMIC_MODEL_INDEX = 2 : Overdamped equations of motions solved with improved euler (2nd order integration scheme)
//The DYNAMIC_MODEL_INDEX can be changed in the ./include/global_configuration.hpp file.

For more information, please refer to the paper associated to SimuCell3D.
*/ 


//------------------------------------------------------------------------------------------------------------------------
//This is an abstract class that defines the interface of the time integration schemes
class time_integration_scheme{

   protected:
      double dt_; //Time step
      double damping_coeff_; //Damping coefficient
      double simulation_time_ = 0.0; //Simulation time
      bool verbose_;

      //True if the time step is a temporary one. Temporary time steps are used by the improved Euler integration scheme
      bool tmp_step_ = false; 

   public:
      time_integration_scheme() = default;                                           //default constructor
      time_integration_scheme(const time_integration_scheme& v) = delete;            //copy constructor
      time_integration_scheme(time_integration_scheme&& v) = delete;                 //move constructor
      time_integration_scheme& operator=(const time_integration_scheme& v) = delete; //copy assignment operator
      time_integration_scheme& operator=(time_integration_scheme&& v) = default;     //move assignment operator

      //Trivial constructor
      time_integration_scheme(const global_simulation_parameters& sim_parameters, const bool verbose) noexcept: 
         dt_(sim_parameters.time_step_), 
         damping_coeff_(sim_parameters.damping_coefficient_),
         verbose_(verbose)
      {};

      double get_simulation_time() const noexcept {return simulation_time_;}

      //Return wether or not the time step is a temporary one (temporary time steps are used by the improved Euler integration scheme)
      bool is_step_tmp() const noexcept {return tmp_step_;}

      void update_nodes_positions(const std::vector<cell_ptr>& cell_lst);



};
//------------------------------------------------------------------------------------------------------------------------





#endif