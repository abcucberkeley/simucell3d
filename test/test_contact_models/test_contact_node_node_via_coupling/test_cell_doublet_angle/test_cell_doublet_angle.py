from os import path
import sys



"""
Simulate a cell doublet getometry, where the two cells have the same apical surface tension and lateral surface tensions.
The contact angle between the cells should follow arccos(lateral_tension / apical_tension).
"""

#--------------------------------------------------------------------------------------------------------
#Get the path to the directory where this script is stored 
path_to_script_folder = path.dirname(path.realpath(__file__))

#Get the path to the simucell3d python library
path_to_build = path.join(path_to_script_folder, "..", "..", "..", "..", "build")
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




#---------------------------------------------------------------------------------------------------------------------------------------------------------
def generate_cell_types_embryo(parameter_set, nb_cells): 
    """
    Each cell will have a different lateral surface tension value. This function 
    returns the simucell3d.cell_type_parameters objects of each cell
    

    Parameters:
    -----------

    parameter_set: list of floats
        List of the parameters used for each cell

    nb_cells: int
        Number of cells in the tissue

    Returns:
    --------
    cell_types: list of simucell3d.cell_type_parameters objects
        List of the simucell3d.cell_type_parameters objects of each cell
    """
    
    cell_type_lst = []

    for cell_id in range(nb_cells):

        #The parameters of the apical faces
        face_type_1 = simucell3d_python_wrapper.face_type_parameters()
        face_type_1.name_                = "apical"     
        face_type_1.face_type_global_id_ = 0
        face_type_1.adherence_strength_  = 0
        face_type_1.repulsion_strength_  = 2e9
        face_type_1.surface_tension_     = 1e-3
        face_type_1.bending_modulus_     = 0

        #The parameters of the lateral faces
        face_type_2 = simucell3d_python_wrapper.face_type_parameters()
        face_type_2.name_                = "lateral"     
        face_type_2.face_type_global_id_ = 1
        face_type_2.adherence_strength_  = 0
        face_type_2.repulsion_strength_  = 2e9
        face_type_2.surface_tension_     = parameter_set[cell_id] 
        face_type_2.bending_modulus_     = 0

        #The other parameters of the cell
        cell_type = simucell3d_python_wrapper.cell_type_parameters()
        cell_type = simucell3d_python_wrapper.cell_type_parameters()
        cell_type.name_                             = "epithelial"     
        cell_type.global_type_id_                   = 0 
        cell_type.mass_density_                     = 1e3
        cell_type.bulk_modulus_                     = 2.5e3   
        cell_type.initial_pressure_                 = 0
        cell_type.max_pressure_                     = 1e20
        cell_type.avg_growth_rate_                  = 0
        cell_type.std_growth_rate_                  = 0  
        cell_type.target_isoperimetric_ratio_       = 150
        cell_type.angle_regularization_factor_      = 1e-16
        cell_type.area_elasticity_modulus_          = 0
        cell_type.surface_coupling_max_curvature_   = 1e20
        cell_type.avg_division_vol_                 = 1e20  
        cell_type.std_division_vol_                 = 0   
        cell_type.min_vol_                          = 5e-17   
        cell_type.add_face_type(face_type_1)
        cell_type.add_face_type(face_type_2)
        cell_type_lst.append(cell_type)

    return cell_type_lst
#---------------------------------------------------------------------------------------------------------------------------------------------------------









#---------------------------------------------------------------------------------------------------------------------------------------------------------
def cell_doublet_test(
        parameter_set, 
        output_folder_path, 
        simulation_id, 
        simulation_duration, 
        failure_loss = 10.0,
        verbose = False
    ):
    """
    Launch a simulation with the variable parameters given in input and return the loss function value.
    The loss is calculated by measuring the deviation between the initial cell shapes and the final cell shapes.

    Parameters:
    -----------

    lateral_surface_tension_lst: list
        List containing the lateral surface tension of each cell

    simulation_id: int
        The id of the simulation

    output_folder_path: str
        The folder where the simulation files will be written

    ground_truth_cell_mesh_lst: list of trimesh.Trimesh objects
        The shape of the cells at end of the simulations if the ground truth parameters were used

    failure_loss: float
        The loss value to return if the simulation fails
    
    """

    nb_cells = 2

    #Create a global_simulation_parameters object to store the simulation parameters
    global_simulation_parameters = simucell3d_python_wrapper.global_simulation_parameters()
    global_simulation_parameters.output_folder_path_ =              path.join(output_folder_path, "simulation_{}".format(simulation_id))
    global_simulation_parameters.input_mesh_path_ =                 "data/input_meshes/2_big_cubes.vtk"
    global_simulation_parameters.damping_coefficient_ =             2e-09
    global_simulation_parameters.simulation_duration_ =             5e-03 
    global_simulation_parameters.sampling_period_ =                 5e-06
    global_simulation_parameters.time_step_ =                       1e-07
    global_simulation_parameters.min_edge_len_ =                    7.5e-7
    global_simulation_parameters.contact_cutoff_adhesion_ =         5e-07
    global_simulation_parameters.contact_cutoff_repulsion_=         5e-07
    global_simulation_parameters.enable_edge_swap_operation_ =      0
    global_simulation_parameters.perform_initial_triangulation_ =   1


    #Generate the cell type of all the cells. Give the same parameter set to all the cells
    cell_type_lst = generate_cell_types_embryo(parameter_set, nb_cells)

    #Launch the simulation
    simulation_outputs = launch_simulation(
        global_simulation_parameters, 
        cell_type_lst,
        verbose=True,
        nb_threads=-1                    
    )

    if simulation_outputs.RETURN_CODE_ != 0: print(simulation_outputs.error_message_)
#---------------------------------------------------------------------------------------------------------------------------------------------------------







#--------------------------------------------------------------------------------------------------------
if __name__ == "__main__":

    #To test the effect of all the loss functions on the parameter estimation process run the following command:
    #sbatch -n 1 --cpus-per-task=6 --time=48:00:00 --job-name="loss_function_embryo" --array=1-250 --mem-per-cpu=1000 --wrap="python3 parameter_estimation_embryo.py \$SLURM_ARRAY_TASK_ID"

    #The ground truth apical surface tensions used to generate the embryo geometry
    ground_truth_lateral_surface_tension_lst = [5e-4, 5e-4]
    #ground_truth_lateral_surface_tension_lst = [1e-3, 1e-3]

    path_to_output_folder = path.join(path_to_build, "cell_doublet_angle_output")

    eval(
        f"""cell_doublet_test(
            parameter_set = ground_truth_lateral_surface_tension_lst,
            output_folder_path = r"{path_to_output_folder}",
            simulation_id = {1},
            simulation_duration = 1e-2,
            verbose = True
        )"""
    )


    
