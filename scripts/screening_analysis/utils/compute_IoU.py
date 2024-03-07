import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import trimesh as tm
from trimesh.voxel import creation
from sys import argv
from filelock import FileLock
from typing import Optional, List, Dict, Union, Iterable

sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))
from mesh_utils import get_cell_geometries_from_file


#---------------------------------------------------------------------------------------------------------------------------------------------------------
def find_nonempty_files(
    path_to_simu_output: str,
    return_all: Optional[bool] = True
) -> Union[List[str], str]:
    """
    It may happen that the simulation outputs contain corrupted or empty file due to 
    issues in the saving process.
    To avoid comparing the initial geometry with such files, this function analyze the
    content of each simulation run's folder (e.g. simulation_outputs/simulation_1/cell_data)
    and returns the file name of the last non-empty file.

    Parameters:
    -----------

    path_to_simu_output: (str)
       The path to the vtk files saved for a given simulation id 
       (e.g., simulation_outputs/simulation_1/cell_data). 
    
    return_all: (Optional[bool] = True)
        If `True`, return all the valid files, otherwise return only the last valid one.

    Returns:
    --------

    (Union[List[str], str]):
        A list of paths to all the non-empty/non-corrupted files if `return_all==True`,
        or the path to the last valid file otherwise.
    """

    # Get all files and sort them by simu idx
    simulation_file_lst = os.listdir(path_to_simu_output)
    simulation_file_ordered_lst = sorted(simulation_file_lst, key= lambda x: int(x.split("_")[-1][:-4]))
    simulation_path_ordered_lst = [
        os.path.join(path_to_simu_output, simulation_file)
        for simulation_file in simulation_file_ordered_lst
    ]

    # Get the size of initial geometry as term of comparison
    init_size = os.path.getsize(simulation_path_ordered_lst[0])

    # Starting from the back, check whether they are empty. Stop once found the first non-empty file
    idx = len(simulation_path_ordered_lst) - 1
    while idx > 0:
        curr_file = simulation_path_ordered_lst[idx]
        file_size = os.path.getsize(curr_file)
        # to be considered non-empty, a file must be at least half of the size of the init geometry
        if file_size > init_size / 2:
            break
        
        idx -= 1

    if return_all:
        return simulation_path_ordered_lst[:(idx + 1)]
    else:
        return simulation_path_ordered_lst[idx]
#---------------------------------------------------------------------------------------------------------------------------------------------------------



#---------------------------------------------------------------------------------------------------------------------------------------------------------
def compute_derivative(
    values:  Iterable[float]
) -> float:
    '''
    Given an array of values, compute its numerical derivative using central difference method.

    Parameters:
    -----------
        
    values: (Iterable[float])
        Array of values to compute the numerical derivative of.

    Returns:
    --------

    derivative: (np.ndarray[float])
        Array containing the numerical derivative values.
    '''

    n = len(values)
    derivative = np.zeros(n)

    # n.b. timestep is the sim_iteration --> h_forward == h_backward == 1
    for i in range(1, n - 1):
        derivative[i] = (values[i + 1] - values[i - 1]) / 2

    # Compute the derivative at the edges using forward/backward differences
    derivative[0] = values[1] - values[0]
    derivative[n - 1] = values[n - 1] - values[n - 2]

    return derivative
#---------------------------------------------------------------------------------------------------------------------------------------------------------



#---------------------------------------------------------------------------------------------------------------------------------------------------------
def compute_iou(
        mesh_1: tm.Trimesh,
        mesh_2: tm.Trimesh, 
        voxel_size: float = 1e-6, 
        # cell_id = None
    ) -> float:
    """
    Voxelize the two meshes and compute the intersection over union between them

    NOTE: about voxel size
    Do not get confused with the voxel size expressed here and the original voxel size.
    In this case we don't want to set the value to the original one, as the meshes are 
    already in the reference system obtained from the original (even anisotropic) voxel 
    size. Therefore, setting a common voxel size on all the axes is not an error, since 
    allows us to stay in the same reference system of the meshes.

    Parameters:
    -----------

    mesh_1: (trimesh.Trimesh)
        The first mesh

    mesh_2: (trimesh.Trimesh)
        The second mesh

    voxel_size: (float)
        The size of the voxel to use for the voxelization

    Returns:
    --------

    iou: (float)
        The intersection over union between the two meshes
    """

    # start1 = time()
    #Get the bounds of the cell initial and final shape
    intial_cell_mesh_bounds = mesh_1.bounds
    final_cell_mesh_bounds =  mesh_2.bounds
    # print(f"Time elapsed for getting bounds: {time() - start1}")

    #Get the minimum and maximum bounds of the cell
    min_x = min(intial_cell_mesh_bounds[0,0], final_cell_mesh_bounds[0, 0])
    min_y = min(intial_cell_mesh_bounds[0,1], final_cell_mesh_bounds[0, 1])
    min_z = min(intial_cell_mesh_bounds[0,2], final_cell_mesh_bounds[0, 2])

    max_x = max(intial_cell_mesh_bounds[1,0], final_cell_mesh_bounds[1,0])
    max_y = max(intial_cell_mesh_bounds[1,1], final_cell_mesh_bounds[1,1])
    max_z = max(intial_cell_mesh_bounds[1,2], final_cell_mesh_bounds[1,2])

    #Get the global bounds of the two meshes
    global_bounds = np.array([[min_x, min_y, min_z],[max_x, max_y, max_z]])

    #Calculate the centroid of the global bounds
    global_centroid = np.mean(global_bounds, axis=0)

    #Get the longest axis of the global bounds
    global_longest_axis_length = np.max(global_bounds[1,:] - global_bounds[0,:])

    radius = int(global_longest_axis_length / voxel_size)

    # start2 = time()
    voxel_grid_1 = creation.local_voxelize(
        mesh_1, 
        point = global_centroid, 
        pitch = voxel_size, 
        radius = radius, 
        fill = True
    )
    # print(f"Time elapsed for local voxelization: {time() - start2}")
    
    voxel_grid_2 = creation.local_voxelize(
        mesh_2, 
        point = global_centroid, 
        pitch = voxel_size, 
        radius = radius, 
        fill = True
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
def shape_difference_calc(
        path_mesh_1: str, 
        path_mesh_2: str,
        cached: Dict[int, List[tm.Trimesh]]
    ) -> float:

    mesh_1_cell_lst = get_cell_geometries_from_file(path_mesh_1, cached)
    mesh_2_cell_lst = get_cell_geometries_from_file(path_mesh_2, cached)

    # assert len(mesh_1_cell_lst) == len(mesh_2_cell_lst), 'One or more cells were lost in these simulations, hence we are going to discard it!'
    loss = []
    try: 
        for cell_id in range(len(mesh_1_cell_lst)):
            iou = compute_iou(mesh_1_cell_lst[cell_id], mesh_2_cell_lst[cell_id])
            loss.append(1 - iou)
    except IndexError:
        print(
            f'Skipping comparison between cell {os.path.basename(path_mesh_1)} and {os.path.basename(path_mesh_2)}'+\
            'since one of the two files may be either empty or corrupted.'
        )

    if loss != []:
        return np.mean(loss), np.std(loss)/np.sqrt(len(mesh_1_cell_lst))
    else:
        return np.nan, np.nan
#---------------------------------------------------------------------------------------------------------------------------------------------------------


#---------------------------------------------------------------------------------------------------------------------------------------------------------
def get_shape_change_wrt_beginning(
        ordered_paths: List[str], 
        cached_data: Dict[int, List[tm.Trimesh]],
        every: Optional[int] = 3 
    ) -> List[float]:

    path_mesh_1 = ordered_paths[0]

    shape_difference_lst = []
    for i in tqdm(range(1, len(ordered_paths), every), desc='Computing shape change wrt beginning'):
        path_mesh_2 = ordered_paths[i]
        shape_diff = shape_difference_calc(path_mesh_1, path_mesh_2, cached_data)
        if shape_diff:
            shape_difference_lst.append(shape_diff)
    
    return shape_difference_lst
#---------------------------------------------------------------------------------------------------------------------------------------------------------


#---------------------------------------------------------------------------------------------------------------------------------------------------------
def get_shape_change_between_iterations(
        ordered_paths: List[str],
        cached_data: Dict[int, List[tm.Trimesh]],
        every: Optional[int] = 3 
    ) -> List[float]:

    path_mesh_1 = ordered_paths[0]

    shape_difference_lst = []

    for i in tqdm(range(1, len(ordered_paths), every), desc='Computing shape change btw iters'):
        path_mesh_2 = ordered_paths[i]
        shape_diff = shape_difference_calc(path_mesh_1, path_mesh_2, cached_data)
        if shape_diff:
            shape_difference_lst.append(shape_diff)
        path_mesh_1 = path_mesh_2
    
    return shape_difference_lst
#---------------------------------------------------------------------------------------------------------------------------------------------------------



#---------------------------------------------------------------------------------------------------------------------------------------------------------
def main(
    simulation_folder_path: str,
    sim_id: int,
    save_dir: str,
    compute_every: Optional[int] = 4
) -> None:
    '''
    Parameters:
    -----------
        simulation_folder_path: (str)
            The absolute path where the meshes of the simulation are saved (cell_data).

        sim_id: (int)
            An identifier of the simulation under analysis.

        save_dir: (str)
            The path where to save the result.

        compute_every: (Optional[int], default=4)
            Specify every how many iterations IoU is computed.
    '''

    assert not simulation_folder_path.startswith('.'), 'The path to simulation data should be absolute!'

    print(f'Analyzing simulation {simulation_folder_path}')

    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    #List all the non-empty files that have been saved during the simulation
    paths_to_simu_files = find_nonempty_files(simulation_folder_path, return_all=True)

    cached_meshes = {}

    #Compute how much the shapes of the cells have changed compare to the beginning of the simulation
    shape_difference_wrt_beginning_lst = get_shape_change_wrt_beginning(
        ordered_paths=paths_to_simu_files,
        cached_data=cached_meshes,
        every=compute_every
    )

    print(f'Shape wrt beginning: {shape_difference_wrt_beginning_lst}')
    print(f'CACHED: {[(k, len(v)) for k, v in cached_meshes.items()]}')

    # Compute the average of the last 10 derivative values (if there are enough, else return NaN)
    mean_values = [mean for mean, std in shape_difference_wrt_beginning_lst]
    if len(mean_values) >= 10:
        avg_IoU_derivative_wrt_beginning = np.mean(compute_derivative(mean_values)[-10:])
    elif 1 < len(mean_values) < 10:
        avg_IoU_derivative_wrt_beginning = np.mean(compute_derivative(mean_values))
    else:
        avg_IoU_derivative_wrt_beginning = None

    # Compute average of last 10 IoU values
    if len(mean_values) >= 10:
        avg_IoU_wrt_beginning = np.mean(mean_values[-10:])
    elif 0 < len(mean_values) < 10:
        avg_IoU_wrt_beginning = np.mean(mean_values)
    else:
        avg_IoU_wrt_beginning = None

    # Write average derivatives to file
    output_file_path = os.path.join(save_dir, 'dashboard_files/IoU_derivative_output.txt')
    with FileLock(output_file_path + ".lock"):
        with open(output_file_path, "a") as file:
            file.write(f'{sim_id} {avg_IoU_derivative_wrt_beginning} \n')

    output_file_path = os.path.join(save_dir, 'dashboard_files/IoU_output.txt')
    with FileLock(output_file_path + ".lock"):
        with open(output_file_path, "a") as file:
            file.write(f'{sim_id} {avg_IoU_wrt_beginning} \n')

    #Compute how much the shapes of the cells have changed compare to the previous iteration
    shape_difference_between_iterations_lst = get_shape_change_between_iterations(
        ordered_paths=paths_to_simu_files,
        cached_data=cached_meshes,
        every=compute_every
    )

    print(f'Shape wrt prev iters: {shape_difference_between_iterations_lst}')
    print(f'CACHED: {[(k, len(v)) for k, v in cached_meshes.items()]}')

    #Plot the results
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 6))
    fig.suptitle(f"{os.path.basename(simulation_folder_path.replace('/cell_data', '')).title()}")
    ax1, ax2 = axes

    #get means and stds
    y = [elem[0] for elem in shape_difference_wrt_beginning_lst]
    y_err = [elem[1] for elem in shape_difference_wrt_beginning_lst]
    x = [(i * compute_every) + 2 for i in range(len(shape_difference_wrt_beginning_lst))]

    ax1.errorbar(x, y, yerr=y_err, fmt='o--', capsize=5, color='blue', ecolor='grey')
    ax1.set_xlabel("File number")
    ax1.set_xticks(x[::2], x[::2], rotation=45)
    ax1.set_ylabel("Shape difference [mean(1-IoU)]")
    ax1.set_title("Shape difference wrt initial geometry")

    #get means and stds
    y = [elem[0] for elem in shape_difference_between_iterations_lst]
    y_err = [elem[1] for elem in shape_difference_between_iterations_lst]
    x = [(i * compute_every) + 1 for i in range(len(shape_difference_between_iterations_lst))]

    ax2.errorbar(x, y, yerr=y_err, fmt='o--', capsize=5, color='orange', ecolor='grey')
    ax2.set_yscale("log")
    ax2.set_xlabel("File number")
    ax2.set_xticks(x[::3], x[::3], rotation=45)
    ax2.set_ylabel("Shape difference [mean(1-IoU)]")
    ax2.set_title("Shape difference between consecutive iterations")

    plt.tight_layout()

    save_plots_dir = os.path.join(save_dir, 'shape_change_plots')
    if not os.path.exists(save_plots_dir):
        os.makedirs(save_plots_dir)

    plt.savefig(os.path.join(save_plots_dir, f"shape_change_{sim_id}.png"))
#---------------------------------------------------------------------------------------------------------------------------------------------------------



if __name__ == "__main__":

    #Check that the name of the screen was given in command line input
    assert argv.__len__() == 4 , "Wrong command line arguments: \n expected python3 shape_change.py [simulation_folder_path] [sim_id] [save_dir]"

    #Get the arguments
    simu_folder_path = argv[1].strip()
    assert os.path.exists(simu_folder_path), "The given screen folder does not exist: " + simu_folder_path
    simu_id = argv[2].strip()
    saving_dir = argv[3].strip()

    main(simu_folder_path, simu_id, saving_dir, 5)