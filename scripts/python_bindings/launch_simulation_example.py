import sys
from os import path, getcwd, _exit
import pandas as pd
from io import StringIO

"""
This python script lauches a simucell3d simulation and retrieves the cell statistics data.
"""


#--------------------------------------------------------------------------------------------------------
#Get the path to the directory where this script is stored 
path_to_script_folder = path.dirname(path.realpath(__file__))

#Get the path to the simucell3d python library
path_to_build = path.realpath(path.join(path_to_script_folder, "..", "..", "build"))
if not path.exists(path_to_build): 
    raise SystemError("\n\nThe given path to the build folder: {} does not exist.\n".format(path_to_build))

#Get the path to the simucell3d python library
path_to_simucell3d_python_library = path.join(path_to_build, "bin", "python_bindings")

if not path.exists(path_to_simucell3d_python_library):
    raise SystemError(
            "\n\nThe path to the simucell3d python library does not exist.\n"+
            "Run the following commands:\n"+
            "cd {}:\n".format(path_to_build) +
            "cmake -DENABLE_PYTHON_BINDINGS=TRUE -DPYTHON_EXECUTABLE=$(which python3) -DCMAKE_BUILD_TYPE=Release .. && make -j6\n"
        )

#Get the path to the library
sys.path.insert(0, path_to_simucell3d_python_library)


#Import the simucell3d python library
import simucell3d_python_wrapper
#---------------------------------------------------------------------------------------------------------------------------------------------------------







#--------------------------------------------------------------------------------------------------------
def get_cell_statistics_data(simulation_outputs):
    """
    This function retrieves the cell statistics data from the simulation_wrapper object.

    Parameters:
    -----------
    simulation_outputs: simucell3d.simulation_outputs
        The sstructure storing the results of the simulation.

    Returns:
    --------
    cell_statistics_data_df: pd.DataFrame
        A Dataframe containing the cell statistics data.
    """


    #Wrap the string data in StringIO function
    string_buffer = StringIO(simulation_outputs.cell_statistics_)
    cell_statistics_data_df = pd.read_csv(string_buffer, sep =",")
    return cell_statistics_data_df
#--------------------------------------------------------------------------------------------------------






#--------------------------------------------------------------------------------------------------------
def launch_simulation(
        global_simulation_parameters, 
        cell_type_lst, 
        verbose=False,
        nb_threads = -1
    ):
    """
    Launch the simulation with the given parameters.

    Parameters:
    -----------
    global_simulation_parameters: simucell3d.global_simulation_parameters
        The structure storing the global simulation parameters.

    cell_type_lst: list
        A list of simucell3d.cell_type_parameters objects.

    perform_initial_triangulation: bool
        If set to True SimuCell3D will try to triangulate the cells at the beginning of the simulation.

    verbose: bool
        If True, print the simulation progress.

    nb_threads: int
        The number of threads used to run a simulation. If set to -1 it will use all the available cores

    Returns:
    --------
    simulation_outputs: simucell3d.simulation_outputs
        The structure storing the results of the simulation.
    """


    #Run the simulation
    simulation_wrapper = simucell3d_python_wrapper.simucell3d_wrapper(
        global_simulation_parameters, 
        cell_type_lst,
        verbose, 
        nb_threads
    )


    #Get the ouput of the simulation
    return simulation_wrapper.get_simulation_outputs()
#--------------------------------------------------------------------------------------------------------






#--------------------------------------------------------------------------------------------------------
if __name__ == "__main__":


    #--------------------------------------------------------------------------------------------------------
    #Create a global_simulation_parameters object to store the simulation parameters
    global_simulation_parameters = simucell3d_python_wrapper.global_simulation_parameters()
    global_simulation_parameters.output_folder_path_ =      path.join(path_to_build, "python_sim_example")
    global_simulation_parameters.input_mesh_path_ =         "data/input_meshes/big_sphere.vtk"
    global_simulation_parameters.damping_coefficient_ =     2e-9
    global_simulation_parameters.simulation_duration_ =     1e-3 
    global_simulation_parameters.sampling_period_ =         1e-5 
    global_simulation_parameters.time_step_ =               1e-07
    global_simulation_parameters.min_edge_len_ =            1e-6
    global_simulation_parameters.contact_cutoff_adhesion_ = 7.5e-7
    global_simulation_parameters.contact_cutoff_repulsion_= 7.5e-7
    global_simulation_parameters.enable_edge_swap_operation_ =  0
    global_simulation_parameters.perform_initial_triangulation_ = 0
    #--------------------------------------------------------------------------------------------------------


    #--------------------------------------------------------------------------------------------------------
    #The parameters of the apical faces
    face_type_1 = simucell3d_python_wrapper.face_type_parameters()
    face_type_1.name_                = "apical"     
    face_type_1.face_type_global_id_ = 0
    face_type_1.adherence_strength_  = 0
    face_type_1.repulsion_strength_  = 1e9
    face_type_1.surface_tension_     = 1e-3
    face_type_1.bending_modulus_     = 0

    #The parameters of the lateral faces
    face_type_2 = simucell3d_python_wrapper.face_type_parameters()
    face_type_2.name_                = "lateral"     
    face_type_2.face_type_global_id_ = 1
    face_type_2.adherence_strength_  = 0
    face_type_2.repulsion_strength_  = 1e9
    face_type_2.surface_tension_     = 5e-4
    face_type_2.bending_modulus_     = 0

    #The other parameters of the cell
    cell_type = simucell3d_python_wrapper.cell_type_parameters()
    cell_type.name_                             = "epithelial"     
    cell_type.global_type_id_                   = 0 
    cell_type.mass_density_                     = 1e3
    cell_type.bulk_modulus_                     = 1e4   
    cell_type.initial_pressure_                 = 0
    cell_type.max_pressure_                     = 1e20
    cell_type.avg_growth_rate_                  = 2e-11
    cell_type.std_growth_rate_                  = 0  
    cell_type.target_isoperimetric_ratio_       = 150
    cell_type.angle_regularization_factor_      = 0
    cell_type.area_elasticity_modulus_          = 0  
    cell_type.surface_coupling_max_curvature_   = 1e7
    cell_type.avg_division_vol_                 = 1e-14  
    cell_type.std_division_vol_                 = 0   
    cell_type.min_vol_                          = 3.7e-17   
    cell_type.add_face_type(face_type_1)
    cell_type.add_face_type(face_type_2)
    cell_type_lst = [cell_type]
    #--------------------------------------------------------------------------------------------------------


    #Launch the simulation
    simulation_outputs = launch_simulation(
        global_simulation_parameters, 
        cell_type_lst,
        verbose=True,
        nb_threads=-1                    
    )

    if simulation_outputs.RETURN_CODE_ == 0:
        #Get the cell statistics data
        cell_statistics_data_df = get_cell_statistics_data(simulation_outputs)

        #Save the cell statistics data
        cell_statistics_data_df.to_csv(path.join(path_to_build, "python_sim_example", "cell_statistics_data.csv"), header = True, index = False)

        print(cell_statistics_data_df)
    
    else:
        print("The simulation did not terminate")
        print(simulation_outputs.error_message_)


