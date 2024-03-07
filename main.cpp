
#include <iostream>

#include "simulation_initializer.hpp"
#include "solver.hpp"
#include "mesh_writer.hpp"




int main (int argc, char** argv){

    //Get the path to the input parameter file
	std::string path_parameter_file;
    if (argc == 2){
        path_parameter_file = argv[1];
    }
    else{
        std::cout << "No input file given in the command line" << std::endl;
        std::cout << "Usage example: ./main path/to/parameter_file.xml" << std::endl;

        return 1;
    }

    //Catch any exception and print error msg
    solver solver_;
    try{
        //Load all the parameters and geometrical information of the tissue
        simulation_initializer sim_init(path_parameter_file);

        //Initialize the solver
        solver_ =  solver(sim_init.get_simulation_parameters(), sim_init.get_cell_lst());

        //Run the simulation
        solver_.run();
    }
    
    //Catch and print any exception
    catch(std::exception const& e){
        std::cerr << e.what() << std::endl;

        return 1;
    }

    return 0;
}