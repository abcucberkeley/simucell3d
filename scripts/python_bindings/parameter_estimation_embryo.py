import sys
from os import path, getcwd, _exit, environ, cpu_count, mkdir
import pandas as pd
#from joblib import Parallel, delayed
import trimesh as tm
from trimesh.voxel import creation
import vtk
import numpy as np
np.random.seed(123)
import matplotlib.pyplot as plt
from skopt import  Optimizer
from skopt.utils import cook_initial_point_generator
import time
from sklearn.linear_model import LinearRegression
from sklearn.metrics import  r2_score, mean_squared_error


from timeout_decorator import *




"""
Parameter estimation pipeline to estimate the differential tension of the cells of a 16 mouse embryo.
Each cell will have a different tension value.
"""

#--------------------------------------------------------------------------------------------------------
#Get the path to the directory where this script is stored 
path_to_script_folder = path.dirname(path.realpath(__file__))

#Get the path to the simucell3d python library
path_to_build = path.join(path_to_script_folder, "..", "..", "build")

if not path.exists(path_to_build):
    raise SystemError(
            "\n\nThe path to the build folder does not exist.\n"
        )

path_to_simucell3d_python_library = path.join(path_to_build, "bin", "python_bindings")

if not path.exists(path_to_simucell3d_python_library):
    raise SystemError(
            "\n\nThe path to the simucell3d python library does not exist.\n"+
            "Make sure that you have compiled the code with python support enabled:\n"+
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
def convert_cell_meshes(simulation_outputs):
    """
    Convert the cell meshes at the end of the simulation to trimesh objects.

    Parameters:
    -----------
    simulation_outputs: simucell3d.simulation_outputs
        The structure storing the results of the simulation.

    Returns:
    --------

    cell_meshes: list
        A list of trimesh objects representing the cell meshes at the end of the simulation.
    """
    
    #Get the cell_lst from the simulation_wrapper object
    cell_lst = simulation_outputs.cell_lst_

    #Create a list to store the cell meshes
    cell_meshes = []

    #Iterate over the cells and convert the meshes to trimesh objects
    for cell in cell_lst:
        tm_cell_mesh = tm.Trimesh(vertices=cell.get_node_coord_lst(), faces=cell.get_face_point_ids(), process=False)
        cell_meshes.append(tm_cell_mesh)
    return cell_meshes
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
        face_type_1.repulsion_strength_  = 1e9
        face_type_1.surface_tension_     = 1e-3
        face_type_1.bending_modulus_     = 0

        #The parameters of the lateral faces
        face_type_2 = simucell3d_python_wrapper.face_type_parameters()
        face_type_2.name_                = "lateral"     
        face_type_2.face_type_global_id_ = 1
        face_type_2.adherence_strength_  = 0
        face_type_2.repulsion_strength_  = 1e9
        face_type_2.surface_tension_     = parameter_set[cell_id] 
        face_type_2.bending_modulus_     = 0

        #The other parameters of the cell
        cell_type = simucell3d_python_wrapper.cell_type_parameters()
        cell_type.name_                             = "epithelial"     
        cell_type.global_type_id_                   = 0 
        cell_type.mass_density_                     = 1e3
        cell_type.bulk_modulus_                     = 1e4   
        cell_type.initial_pressure_                 = parameter_set[nb_cells + cell_id]
        cell_type.max_pressure_                     = 1e20
        cell_type.avg_growth_rate_                  = 0
        cell_type.std_growth_rate_                  = 0  
        cell_type.target_isoperimetric_ratio_       = 150
        cell_type.angle_regularization_factor_      = 0
        cell_type.area_elasticity_modulus_          = 0  
        cell_type.surface_coupling_max_curvature_   = 5e7
        cell_type.avg_division_vol_                 = 1e20  
        cell_type.std_division_vol_                 = 0   
        cell_type.min_vol_                          = 5e-17   
        cell_type.add_face_type(face_type_1)
        cell_type.add_face_type(face_type_2)
        cell_type_lst.append(cell_type)

    return cell_type_lst
#---------------------------------------------------------------------------------------------------------------------------------------------------------


#---------------------------------------------------------------------------------------------------------------------------------------------------------
def convert_vtk_cell_to_trimesh(vtk_cell):
    """
    Convert a vtk cell to a trimesh cell

    Parameters:
    -----------
    vtk_cell: vtk.vtkCell
        The vtk cell to convert

    Returns:
    --------
    trimesh_cell: trimesh.Trimesh
        The trimesh cell
    """


    #Get the point id of all the points of the cell
    global_point_id_lst = [vtk_cell.GetPointId(i) for i in range(vtk_cell.GetNumberOfPoints())]

    #Get the coordinates of all the points of the cell
    point_coordinate_lst = [vtk_cell.GetPoints().GetPoint(i) for i in range(vtk_cell.GetNumberOfPoints())]

    #Create a map between the global point id and the local point ids
    global_to_local_point_id_map = {global_point_id_lst[i] : i for i in range(len(global_point_id_lst))}

    faces_lst = []

    #Loop over the faces of the cell
    for i in range(vtk_cell.GetNumberOfFaces()):
            
            #Get the face
            face = vtk_cell.GetFace(i)
    
            #Get the number of points of the face
            num_points = face.GetNumberOfPoints()
    
            #Create a list to store the points of the face
            face_points_lst = []
    
            #Loop over the points of the face
            for j in range(num_points):
    
                #Get the global id of the point
                global_point_id = face.GetPointId(j)

                #Convert the global id to the local id
                local_point_id = global_to_local_point_id_map[global_point_id]
    
                #Append the point to the list
                face_points_lst.append(local_point_id)
    
            #Append the face to the list
            faces_lst.append(face_points_lst)


    #Create the trimesh cell
    trimesh_cell = tm.Trimesh(point_coordinate_lst, faces_lst)

    return trimesh_cell
#---------------------------------------------------------------------------------------------------------------------------------------------------------


#---------------------------------------------------------------------------------------------------------------------------------------------------------
def get_cell_geometries_from_file(path_to_input_mesh):
    """
    Load the initial mesh file and return each cell as a trimesh object in a list

    Parameters:
    -----------
    path_to_input_mesh: string
        The path to the input mesh file

    Returns:
    --------
    cell_lst: list of trimesh.Trimesh objects
        The list of the cells
    """

    reader = vtk.vtkUnstructuredGridReader()
    reader.SetFileName(path_to_input_mesh)
    reader.Update()
    ugrid = reader.GetOutput()

    cell_mesh_lst =  []
    for i in range(ugrid.GetNumberOfCells()):
        vtk_cell = ugrid.GetCell(i)
        trimesh_cell_mesh = convert_vtk_cell_to_trimesh(vtk_cell)
        cell_mesh_lst.append(trimesh_cell_mesh)

    return cell_mesh_lst
#---------------------------------------------------------------------------------------------------------------------------------------------------------


#---------------------------------------------------------------------------------------------------------------------------------------------------------
def compute_iou(mesh_1, mesh_2, voxel_size=1e-6, cell_id = None):
    """
    Voxelize the two meshes and compute the intersection over union between them

    Parameters:
    -----------
    mesh_1: trimesh.Trimesh
        The first mesh

    mesh_2: trimesh.Trimesh
        The second mesh

    voxel_size: float
        The size of the voxel to use for the voxelization

    Returns:
    --------
    iou: float
        The intersection over union between the two meshes
    """


    #Get the bounds of the cell initial and final shape
    initial_cell_mesh_bounds = mesh_1.bounds
    final_cell_mesh_bounds =  mesh_2.bounds

    #Get the minimum and maximum bounds of the cell
    min_x = min(initial_cell_mesh_bounds[0,0], final_cell_mesh_bounds[0, 0])
    min_y = min(initial_cell_mesh_bounds[0,1], final_cell_mesh_bounds[0, 1])
    min_z = min(initial_cell_mesh_bounds[0,2], final_cell_mesh_bounds[0, 2])

    max_x = max(initial_cell_mesh_bounds[1,0], final_cell_mesh_bounds[1,0])
    max_y = max(initial_cell_mesh_bounds[1,1], final_cell_mesh_bounds[1,1])
    max_z = max(initial_cell_mesh_bounds[1,2], final_cell_mesh_bounds[1,2])

    #Get the global bounds of the two meshes
    global_bounds = np.array([[min_x, min_y, min_z],[max_x, max_y, max_z]])

    #Calculate the centroid of the global bounds
    global_centroid = np.mean(global_bounds, axis=0)

    #Get the longest axis of the global bounds
    global_longest_axis_length = np.max(global_bounds[1,:] - global_bounds[0,:])

    radius = int(global_longest_axis_length/ voxel_size)

    voxel_grid_1 = creation.local_voxelize(
        mesh_1, 
        point = global_centroid, 
        pitch = voxel_size, 
        radius = radius, 
        fill=True
    )
    
    voxel_grid_2 = creation.local_voxelize(
        mesh_2, 
        point = global_centroid, 
        pitch = voxel_size, 
        radius = radius, 
        fill=True
    )

    #That can happen in case of numerical instabilities during the simulation
    if voxel_grid_1 == None or voxel_grid_2 == None: return 0.

    voxel_grid_1 = voxel_grid_1.matrix.astype(bool)
    voxel_grid_2 = voxel_grid_2.matrix.astype(bool)

    #Compute the intersection over union between the two voxel grids
    iou = np.sum(np.logical_and(voxel_grid_1, voxel_grid_2)) / np.sum(np.logical_or(voxel_grid_1, voxel_grid_2))
    return iou
#---------------------------------------------------------------------------------------------------------------------------------------------------------


#---------------------------------------------------------------------------------------------------------------------------------------------------------
def call_back_function(
        nb_calls_lst, 
        loss_lst, 
        avg_relative_error_lst,
        min_loss_lst, 
        min_avg_relative_error_lst,
        time_per_iteration_lst, 
        output_path,
        title = "None"
    ):
    """
    Call back function used to print the progress of the minimization

    Parameters:
    ----------

    nb_calls_lst: list of int
        Number of time the loss function has been called

    loss_lst: list of floats
        The loss value for each call of the loss function

    avg_relative_error_lst: list of floats
        The average deviation from the ground truth for each call of the loss function

    min_loss_lst: list of floats
        The minimum loss after nb update of the surrogate model

    time_per_iteration_lst: list of float
        Time in minute that it takes to update the model

    output_path: str
        Path to where the graph should be stored
    
    Retunrs:
    -------

    None
    
    """

    #Create a matplotlib plot to show the progress
    fig = plt.figure(figsize=(12, 8))

    ax1 = fig.add_subplot(231)
    ax1.set_xlabel("nb calls")
    ax1.set_ylabel("min loss value recorded")
    ax1.plot(nb_calls_lst, min_loss_lst,'-o',  c = "b")
    ax1.set_xlim(left = nb_calls_lst[0])

    ax2 = fig.add_subplot(232)
    ax2.set_xlabel("nb calls")
    ax2.set_ylabel("average relative error")
    ax2.plot(nb_calls_lst, min_avg_relative_error_lst,'-o',  c = "b")
    ax2.set_xlim(left = nb_calls_lst[0])

    #Plot the correlation between the loss value and the distance in parameter space
    ax3 = fig.add_subplot(233)
    ax3.scatter(avg_relative_error_lst, loss_lst,  c = "k")
    ax3.set_xlabel("average relative error")
    ax3.set_ylabel("loss value")

    #If there are more than 3 points we compute a linear regression
     
    if len(nb_calls_lst) > 3:

        avg_relative_error_ar = np.array(avg_relative_error_lst).reshape(-1, 1)
        loss_ar = np.array(loss_lst).reshape(-1, 1)

        #Get the rows corresponding to failed simulations
        failed_simulation_idx = np.where(loss_ar == 10.0)[0]

        #Remove the failed simulations from the regression
        avg_relative_error_ar = np.delete(avg_relative_error_ar, failed_simulation_idx, axis = 0)
        loss_ar = np.delete(loss_ar, failed_simulation_idx, axis = 0)
        
        reg = LinearRegression().fit(avg_relative_error_ar, loss_ar)
        reg_x = np.linspace(np.min(avg_relative_error_ar), np.max(avg_relative_error_ar), 100)
        reg_y = reg.predict(reg_x.reshape(-1, 1))

        ax3.plot(reg_x, reg_y, color = "r", linestyle = "--")

        #Print the equation of the line on the plot
        a = reg.coef_[0][0]
        b = reg.intercept_[0]

        r2 = r2_score(loss_ar, reg.predict(avg_relative_error_ar))
        mse = mean_squared_error(loss_ar, reg.predict(avg_relative_error_ar))

        ax3.text(0.5, 0.85, "r^2:   {:.2f},  mse: {:.1e}".format(r2, mse), horizontalalignment='center', verticalalignment='center', transform = plt.gca().transAxes)


    ax4 = fig.add_subplot(234)
    ax4.set_xlabel("nb calls")
    ax4.set_ylabel("time per iteration (min)")
    ax4.plot(nb_calls_lst, time_per_iteration_lst,'-o',  c = "r")
    
    #Convert the time to hours
    time_per_iteration_h_lst = list(map(lambda x : x / 60., time_per_iteration_lst))

    #Cumulate the time per iteration to get the total running time
    total_running_time = np.cumsum(time_per_iteration_h_lst).tolist()

    ax5 = fig.add_subplot(235)
    ax5.set_xlabel("nb calls")
    ax5.set_ylabel("Total time (h)")
    ax5.plot(nb_calls_lst, total_running_time,'-o',  c = "r")
    if title != None: fig.suptitle(title)


    fig.tight_layout()
    plt.savefig(output_path)
    plt.close()
#---------------------------------------------------------------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------------------------------------------------------------
#After 10 minutes the simulation is considered as failed
#@timeout(seconds = 10*60, default_value=10.0)
def loss_function_1(
        parameter_set, 
        ground_truth_cell_mesh_lst,
        output_folder_path, 
        job_id, 
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

    nb_cells = ground_truth_cell_mesh_lst.__len__()

    simulation_duration = 2.5e-4

    #Create a global_simulation_parameters object to store the simulation parameters
    global_simulation_parameters = simucell3d_python_wrapper.global_simulation_parameters()
    global_simulation_parameters.output_folder_path_ =      path.join(output_folder_path, "simulation_{}".format(job_id))
    global_simulation_parameters.input_mesh_path_ =         "data/input_meshes/embryo_5.vtk"
    global_simulation_parameters.damping_coefficient_ =     2e-9
    global_simulation_parameters.simulation_duration_ =     simulation_duration 
    global_simulation_parameters.sampling_period_ =         2.5e-6
    global_simulation_parameters.time_step_ =               1e-07
    global_simulation_parameters.min_edge_len_ =            7.5e-07
    global_simulation_parameters.contact_cutoff_adhesion_ = 5e-07
    global_simulation_parameters.contact_cutoff_repulsion_= 5e-07
    global_simulation_parameters.enable_edge_swap_operation_ =  0
    global_simulation_parameters.perform_initial_triangulation_ = 0

    #Generate the cell type of all the cells. Give the same parameter set to all the cells
    cell_type_lst = generate_cell_types_embryo(parameter_set, nb_cells)

    #Launch the simulation
    simulation_outputs = launch_simulation(
        global_simulation_parameters, 
        cell_type_lst,
        verbose=verbose,
        nb_threads=-1                    
    )

    #If the simulation was unstable or some cells shrinked below the min volume, return a high loss
    if (simulation_outputs.RETURN_CODE_ != 0) or (simulation_outputs.cell_lst_.__len__() != ground_truth_cell_mesh_lst.__len__()): 
        loss = failure_loss
        print(simulation_outputs.error_message_)
    else:

        #Convert the final cell meshes to trimesh objects
        final_cell_mesh_lst = convert_cell_meshes(simulation_outputs)

        #Calculate the intersection over union between the initial and final cell shapes
        loss = 0.
        assert len(ground_truth_cell_mesh_lst) == len(final_cell_mesh_lst)
        for cell_id  in range(len(ground_truth_cell_mesh_lst)):
            iou = compute_iou(ground_truth_cell_mesh_lst[cell_id], final_cell_mesh_lst[cell_id], cell_id = cell_id)
            loss += 1. - iou
    return loss
#---------------------------------------------------------------------------------------------------------------------------------------------------------



#---------------------------------------------------------------------------------------------------------------------------------------------------------
def launch_loss_function_minimization(
        job_id,
        loss_function,
        path_to_screen_output_folder, 
        ground_truth_cell_mesh_lst,
        ground_truth_parameter_set_lst,
        nb_surrogate_update = 1000, 
        aquisition_function = "LCB",
        xi = 0.01,
        kappa = 1.96
    ):
    """Function that contains all the minimization process
    
    Parameters:
    -----------
    path_to_screen_output_folder: str
        The path to the folder where the screen output will be written

    simulation_duration: float
        The duration of each simulation

    nb_surrogate_update: int
        The number of time the surrogate model will be updated

    aquisition_function: str
        The aquisition function to use for the bayesian optimization (EI, PI, LCB)

    xi : float
        The xi parameter of the aquisition function

    kappa : float
        The kappa parameter of the aquisition function

    nb_threads: int
        The number of threads to used to run simulations simultaneously. Each simulation is run in a separate thread
    """

    if not path.exists(path_to_screen_output_folder): mkdir(path_to_screen_output_folder)
    nb_cells = len(ground_truth_cell_mesh_lst)

    #The initial guess of the cell parameters

    #Dataframe storing a summary of the minimization process


    summary_df = pd.DataFrame(columns = [f"tension_{i}" for i in range(nb_cells)] + [f"pressure_{i}" for i in range(nb_cells)] + ["loss"])
    summary_df.to_csv(path.join(path_to_screen_output_folder, "loss_function_summary.csv"), header = True, index = True)

    #The minimum loss recorded 
    min_loss = float("inf")

    #Also store the minimum distance to the optimal parameter set
    ground_truth_parameter_set_ar = np.array(ground_truth_parameter_set_lst)

    min_avg_relative_error = float("inf")

    #The id of the simulation where the loss has been minimized the most
    best_simulation = None

    #Keep track of the min loss at each update of the surrogate model
    min_loss_lst = [] 
    nb_calls_lst = []
    min_avg_relative_error_lst = []

    #Keep track of the loss and the average deviation from the ground truth of each tested parameter set
    loss_lst = []
    avg_relative_error_lst = []

    #Save how much time it takes to update the surrogate model at each iteration
    time_per_iteration_lst = [] 

    #The parameter space that will be explored by the minimizer
    parameter_bounds_lst = [(2.5e-4, 1e-3) for _ in range(nb_cells)] + [(100, 200) for _ in range(nb_cells)] 

    #Create the optimizer
    lhs_maximin = cook_initial_point_generator("lhs", criterion="maximin")
    opt = Optimizer(
        base_estimator = "GP",                                                                  # Use a Gaussian process as surrogate model
        dimensions = parameter_bounds_lst,                                                      # The bounds of each parameter
        acq_func=aquisition_function,                                                           # The acquisition function (Expected improvement)
        n_jobs = -1,                                                                            # The number of jobs used to minimize the surrogate model (-1 == all process available)
        acq_optimizer = "lbfgs",                                                                # Use a minimizer to minimize the surrogate model (Gaussian process) instead of sampling
        initial_point_generator = lhs_maximin, #"lhs",                                          # Method to draw the initial samples                  
        n_initial_points= 50,                                                                   # The number of initialization points                       
        random_state=job_id,                                                                    # Seed                                                
        acq_func_kwargs = {"xi": xi, "kappa": kappa},                                           # The exploration versus exploitation parameter                                                                                                                                                                                                                                                           # Seed                                                                                                   
        acq_optimizer_kwargs = {"maxiter": 200, "n_restarts_optimizer" : 30}                    # The maximum number of iteration of the minimizer (lbfgs)                                                                                                                                                                                                                                                                                                                    # Seed
    )    

    for iteration in range(nb_surrogate_update):

        start = time.time()

       #We ask one parameter for each cell
        x_lst = opt.ask()
     
        #We get in return one loss for each cell
        y = loss_function(
            parameter_set = x_lst, 
            ground_truth_cell_mesh_lst =  ground_truth_cell_mesh_lst,
            output_folder_path = path_to_screen_output_folder, 
            simulation_id = iteration
        )

        #Store the parameter values tested and their associated loss
        x_ar = np.array(x_lst)
        y_ar = np.array([y])

        #Store the parameter values tested and their associated loss in a CSV file
        summary_iteration_ar = np.hstack([x_ar, y_ar])
        summary_iteration_ar = np.expand_dims(summary_iteration_ar, axis = 0)
        summary_iteration_df = pd.DataFrame(data = summary_iteration_ar, columns= summary_df.columns, index = [iteration])
        summary_iteration_df.to_csv(path.join(path_to_screen_output_folder, "loss_function_summary.csv"), mode = "a", header = False, index = True)

        #Compute the average relative deviation from the ground truth
        avg_relative_error = np.mean([abs(x_ar[i] - ground_truth_parameter_set_ar[i]) / ground_truth_parameter_set_ar[i] for i in range(len(x_ar))])

        #Store the loss and the average deviation from the ground truth
        loss_lst.append(y)
        avg_relative_error_lst.append(avg_relative_error)

        #If the loss recoreded at this iteration is the minimum one 
        if y < min_loss: 

            #Keep track of the smallest lost and its devriation from the ground truth
            min_loss = y
            best_simulation = iteration
            min_avg_relative_error = avg_relative_error

        #Update the optimizer
        opt.tell(x_lst, y) 

        end = time.time()
        time_per_iteration_lst.append((end - start) / 60.)

        #Call the call back function to print the progress
        min_loss_lst.append(min_loss)
        min_avg_relative_error_lst.append(min_avg_relative_error)
        nb_calls_lst.append((iteration+1))
        
        #Create a plot of the minimization progress
        call_back_function(
            nb_calls_lst, 
            loss_lst, 
            avg_relative_error_lst,
            min_loss_lst, 
            min_avg_relative_error_lst,
            time_per_iteration_lst, 
            output_path = path.join(getcwd(), f"minimization_progress_{job_id}.png"),
            title = f"Minimization progress {job_id}"
        )
        
        print(f"Nb calls: {(iteration+1)}, Min loss: {min_loss}, Avg deviation from ground truth: {min_avg_relative_error}, Best simulation: {best_simulation}")
#---------------------------------------------------------------------------------------------------------------------------------------------------------









#--------------------------------------------------------------------------------------------------------
if __name__ == "__main__":


    #To test the effect of all the loss functions on the parameter estimation process run the following command:
    #sbatch -n 1 --cpus-per-task=6 --time=48:00:00 --job-name="loss_function_embryo" --array=1-10 --mem-per-cpu=1000 --wrap="python3 parameter_estimation_embryo.py \$SLURM_ARRAY_TASK_ID"

    #To just test the first loss function run the following command:
    #sbatch -n 1 --cpus-per-task=6 --time=48:00:00 --job-name="loss_function_embryo" --mem-per-cpu=1000 --wrap="python3 parameter_estimation_embryo.py 1"


    #Get the id of the simulation
    sim_id = 0
    
    ground_truth_parameter = [5e-4 for i in range(16)] + [116, 113.2, 118.5, 112.5, 118.4, 125.3, 117.6, 110.8, 121.8, 118.1, 119.6, 115.5, 116.9, 108.8, 122.4, 124.6] 

    #Get what should be the final geometry of the cells if the parameter set used is the ground truth
    ground_truth_cell_mesh_lst =  get_cell_geometries_from_file(path.join(path_to_build, "..", "data/input_meshes/embryo_5.vtk"))

    loss_function_1(
        ground_truth_parameter, 
        ground_truth_cell_mesh_lst,
        path.join(path_to_build, "parameter_estimation_embryo"), 
        sim_id, 
        failure_loss = 10.0,
        verbose = True
    )
        
    _exit(0)

    eval(
        f"""launch_loss_function_minimization(
            job_id = {sim_id},
            loss_function = loss_function_{sim_id},
            path_to_screen_output_folder = r"{path.join(path_to_build, "parameter_estimation_embryo")}",
            ground_truth_cell_mesh_lst = ground_truth_cell_mesh_lst,
            ground_truth_parameter_set_lst = ground_truth_parameter,
            nb_surrogate_update = 1000, 
            aquisition_function = "LCB",
            xi = 0.01,
            kappa = 1.96
        )"""
    )


    #launch_loss_function_minimization(
    #    0,
    #    path.join(path_to_build, "loss_function_embryo"), 
    #    ground_truth_cell_mesh_lst,
    #    ground_truth_lateral_surface_tension_lst,
    #    simulation_duration = 2.5e-4,
    #    nb_surrogate_update = 1000, 
    #    aquisition_function = "LCB",
    #    xi = 0.01,
    #    kappa = 1.96
    #)

    #loss_value = loss_function_embryo(
    #    parameter_set = ground_truth_lateral_surface_tension_lst,
    #    ground_truth_cell_mesh_lst = ground_truth_cell_mesh_lst,
    #    output_folder_path = path.join(path_to_build, "loss_function_embryo_5"),
    #    simulation_id = 1,
    #    simulation_duration = 1e-2,
    #    verbose = True
    #)
    #_exit(0)

    #sim_id = 0  
    ##Create some random parameters that deviates from the ground truth by the given percentage
    #for expected_devation_from_ground_truth_percentage in  [0., 0.1, 0.25, 0.4]:
#
    #    #Generate a random parameter
    #    random_parameter_set = [
    #        ground_truth_parameter[i] + ground_truth_parameter[i] * expected_devation_from_ground_truth_percentage * np.random.uniform(-3, 3) for i in range(len(ground_truth_parameter))
    #    ]
    #    
    #    #Calculate the true devriation from the ground truth
    #    true_deviation_from_ground_truth_percentage = np.mean([abs(random_parameter_set[i] - ground_truth_parameter[i]) / ground_truth_parameter[i] for i in range(len(random_parameter_set))])
    #    print("\n")
    #    for sim_duration in [1e-5, 1e-4, 2.5e-4, 5e-4]:
#
    #        #Compute the time in minutes needed to run the simulation
    #        start = time.time()
#
    #        loss_value = loss_function_embryo_0(
    #            parameter_set = random_parameter_set,
    #            ground_truth_cell_mesh_lst = ground_truth_cell_mesh_lst,
    #            output_folder_path = path.join(path_to_build, "loss_function_embryo"),
    #            simulation_id = sim_id,
    #            simulation_duration = sim_duration,
    #            verbose = False
    #        )
#
    #        end = time.time()
    #        duration = (end - start) / 60.            
    #        print("sim id {}, deviation {:.2f}, sim duration {:.1e}, computation duration {:.2f} min, loss {:.2f}".format(sim_id, true_deviation_from_ground_truth_percentage, sim_duration, duration, loss_value))
    #        sim_id += 1