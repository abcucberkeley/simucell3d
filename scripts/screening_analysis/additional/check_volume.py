import os
import trimesh as tm
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm
from typing import Optional, Tuple
import sys

sys.path.append("/home/fcarrara/Documents/Simulations/SimuCell3D_v2/scripts/screening_analysis/")
print(os.getcwd())
from utils import get_cell_geometries_from_file
    
#----------------------------------------------------------------------------------------
def compute_volumes_in_geometry(
        path_to_geometry_file: str,
        remove_shell: Optional[bool] = True
) -> Tuple[int, int]:
    '''
    Given the path to a .vtk file containing a list of cell meshes, return a list
    of indexes of colliding cells.

    Parameters:
    -----------
        path_to_geometry_file: (str)
            The path to a .vtk file containing a list of cell meshes.
        
        remove_shell: (Optional[bool], default=True)
            If `True` remove the cell mesh corresponding to the shell
            (i.e., the largest mesh).
    
    Returns:
    --------
        (Tuple[int, int])
            A tuple of indexes of the first 2 interpenetrating cells
            (`None` if no interpenetrating cells wer found).
    '''
    cell_mesh_lst = get_cell_geometries_from_file(path_to_geometry_file, None)

    if remove_shell:
        # remove largest mesh (shell)
        num_vertices = [len(mesh.vertices) for mesh in cell_mesh_lst]
        try:
            largest_mesh_idx, _ = max(enumerate(num_vertices), key=lambda x: x[1])
        except:
            print(f'Number of cells: {len(cell_mesh_lst)}')
            print(f'Number of vertices per mesh: {num_vertices}')
        del cell_mesh_lst[largest_mesh_idx]

    return [x.volume for x in cell_mesh_lst]
#----------------------------------------------------------------------------------------



#----------------------------------------------------------------------------------------
def get_simulation_run_volumes(
        path_to_mesh_dir: str,
) -> Tuple[np.ndarray[float], np.ndarray[float]]:
    '''
    Given a path to the directory that stores .vtk files for the different
    simulation iterations (e.g., ./simulation_outputs/simulation_{id}/cell_data),
    returns `True` if at any iteration, two meshes interpenetrate each other. 

    NOTE:
    The function stops checking once reached an iteration showing interpenetration.
    
    Parameters:
    -----------
        path_to_mesh_dir: (str)
            The path to the directory that stores .vtk files for the different
            simulation iterations (e.g., ./simulation_outputs/simulation_{id}/cell_data)

    Returns:
    --------
        The means, std errors, and the arrays of cell volumes across iterations 
    '''
    # Note: it is more likely to  have interpenetration at the last iterations.
    # Therefore, we start our search from there
    num_files = len(os.listdir(path_to_mesh_dir))

    all_volumes, vol_means, vol_errs = [], [], []
    for i in tqdm(range(1, num_files+1), desc="Computing volume"):
        mesh_file = f'result_{i}.vtk'

        curr_volumes = compute_volumes_in_geometry(os.path.join(path_to_mesh_dir, mesh_file))
        if np.any(np.asarray(curr_volumes) < 0):
            continue

        all_volumes.append(curr_volumes)
        vol_means.append(np.mean(curr_volumes))
        vol_errs.append(np.std(curr_volumes) / np.sqrt(len(curr_volumes)))
    
    return np.asarray(vol_means), np.asarray(vol_errs), np.asarray(all_volumes)
#---------------------------------------------------------------------------------------------------------


if __name__ == "__main__":

    volume_means, volume_errs, volumes = get_simulation_run_volumes(
        "/home/fcarrara/Documents/Simulations/screening_results/bronchiole_screen/simu_length_screening_v9/simulation_outputs/sim_50/cell_data/"
    )

    assert np.all(volumes >= 0), f"Volume cannot be negative!!, {volumes < 0}"

    # AVERAGE VOLUME
    x = np.arange(len(volume_means))

    fig, ax = plt.subplots()
    sns.lineplot(x=x, y=volume_means, ax=ax)
    ax.fill_between(x, volume_means - volume_errs, volume_means + volume_errs, alpha=0.2)
    ax.set_xlabel('Simulation id')
    ax.set_ylabel('Mean Volume')
    ax.set_title('Evolution of Average Volume')

    fig.savefig('./results/volume_evolution.png')

    # SINGLE CELLS VOLUMES
    num_cells = volumes.shape[1]

    num_rows = num_cells // 3 + 1
    num_cols = 3
    fig, axes = plt.subplots(nrows=num_rows, ncols=num_cols, figsize=(5*num_rows, 5*num_cols))
    fig.suptitle('Evolution of single cell volumes', fontsize=22)
    x = np.arange(1, volumes.shape[0] + 1)

    for i in range(num_cells):
        row = i // 3
        col = i % 3
        ax = axes[row, col]
        ax.plot(x, volumes[:, i], linewidth=2)
        ax.set_ylabel('Volume', fontsize=14)
        ax.set_title(f'Evolution of cell {i} volume', fontsize=18)

    plt.subplots_adjust(top=0.9, hspace=0.4, wspace=0.4)

    fig.savefig('./results/single_cells_volume_evolution.png')
