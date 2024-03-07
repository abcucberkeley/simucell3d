import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, ScalarFormatter
from sys import argv
from tqdm import tqdm
from typing import Tuple, Dict, Union
from utils.mesh_utils import get_cell_geometries_from_file
from utils.compute_IoU import find_nonempty_files, compute_iou

"""
Run this script to compute IoU between initial and final geometry for different 
configurations of adherence strength and surface tension.
This allows to find the optimal range for these parameters.

NOTE: This script collects similar data to the one needed for the screening dashboard
but it is more optimized for this task.
"""

def get_screening_parameters(
    screening_table_path: str
) -> Dict[str, np.ndarray]:
    """
    Extract arrays of screened parameters ("surface_tension" and "adherence_strength" in this case)

    Parameters:
    -----------

    screening_table_path: (str)
        Path to the .csv file used for setting up the parameter screening.

    Returns:
    --------
    
    (Dict[str, np.ndarray]):
        A dictionary that associates to each parameter name the list of screened values.
    """

    # Load csv file
    screen_table = pd.read_csv(screening_table_path)

    out_dict = {}
    out_dict["surface_tension"] = screen_table["epi_api_surface_tension"].values
    out_dict["adherence_strength"] = screen_table["epi_api_adherence_strength"].values

    return out_dict


def compute_start_vs_end_IoU_loss(
    path_to_simu_outputs: str,
) -> Tuple[np.ndarray[float], np.ndarray[int]]:
    """
    Parameters:
    -----------

    path_to_simu_outputs: (str)
        The path to "simulation_outputs" directory

    Returns:
    --------
    
    iou_loss_vals: (np.ndarray[float])
        An array that store the (1 - IoU) values for all the simulation runs

    num_iters: (np.ndarray[int])
        An array storing the number of valid iterations for each simulation run
    """

    # Get all the subdirectories related to single simulation runs
    simu_subdirs = os.listdir(path_to_simu_outputs)
    simu_subdirs = sorted(simu_subdirs, key=lambda x: int(x.split("_")[1]))

    # Load the initial geometry
    path_to_init_geom = os.path.join(path_to_simu_outputs, simu_subdirs[0], "cell_data/result_1.vtk")
    init_meshes_lst = get_cell_geometries_from_file(path_to_init_geom, None)

    # Iterate over simulation subdirs to compute IoU btw initial and final geometries
    iou_loss_vals = np.empty(len(simu_subdirs))
    num_iters = np.ones(len(simu_subdirs))
    for i in range(len(simu_subdirs)):
        print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        print(f"Computing (1 - IoU) for {simu_subdirs[i]}")
        # Get path to last valid result vtk file
        path_to_cell_data = os.path.join(path_to_simu_outputs, simu_subdirs[i], "cell_data")
        path_to_last_file = find_nonempty_files(
            path_to_simu_output=path_to_cell_data,
            return_all=False
        )
        last_sim_id = int(os.path.basename(path_to_last_file).split(".")[0].split("_")[-1])
        num_iters[i] = last_sim_id

        if last_sim_id == 1:
            iou_loss_vals[i] = np.nan
        else:
            # Get meshes in last valid geometry
            final_meshes_lst = get_cell_geometries_from_file(path_to_last_file, None)

            # Check that the number of cells is the same as in the initial geometry
            if len(init_meshes_lst) != len(final_meshes_lst):
                print((
                    f"In {simu_subdirs[i]}, the final geometry ({len(final_meshes_lst)}) has a different number "
                    f"of cells wrt the initial geometry ({len(init_meshes_lst)}), hence IoU cannot be computed."))
                iou_loss_vals[i] = np.nan
                continue
            
            # Compute the IoU average across all the cells
            curr_iou_losses = np.zeros(len(final_meshes_lst))
            for j, [init_mesh, final_mesh] in tqdm(enumerate(zip(init_meshes_lst, final_meshes_lst)), total=len(init_meshes_lst)):
                curr_iou_losses[j] = 1 - compute_iou(mesh_1=init_mesh, mesh_2=final_mesh)

            iou_loss_vals[i] = np.mean(curr_iou_losses)

    return iou_loss_vals, num_iters


def plot_IoU_loss(
    iou_losses: np.ndarray[float],
    screening_param_dict: Dict[str, np.ndarray],
    num_iterations: np.ndarray[int],
    out_dir: Union[str, None]
) -> None:
    """
    Make a scatter plot of "surface_tension" vs "adherence_strength" with points' color 
    determined by the value of IoU loss.

    Parameters:
    -----------

    iou_losses: (np.ndarray[float])
        The (1 -IoU) average values between initial and final geometry for each simulation
        run in the screening.
    
    num_iterations: (np.ndarray[int])
        An array storing the number of valid iterations for each simulation run

    screening_param_dict: (Dict[str, np.ndarray])
        A dictionary that associates to each parameter name the list of screened values

    out_dir: (str)
        The directory in which the plot is saved. If `None` the plot is not saved.
    """

    # Get screening parameters
    x = screening_param_dict["surface_tension"]
    y = screening_param_dict["adherence_strength"]

    # Get indexes of failed simulations (stopped at first iter)
    failed_mask = np.isnan(iou_losses)
    x_failed, y_failed = x[failed_mask], y[failed_mask]
    x_valid, y_valid, iou_losses_valid =  x[~failed_mask], y[~failed_mask], iou_losses[~failed_mask]

    # Make plot
    fig = plt.figure(figsize=(27, 10), layout="constrained")
    fig.suptitle("Analysis of mechanical properties screening", fontsize=24)

    # Add (1 - IoU) plot
    ax1 = fig.add_subplot(1, 3, 1)
    sc = ax1.scatter(x_valid, y_valid, c=iou_losses_valid, cmap='viridis', s=100)
    #plot failed param config as crosses
    ax1.scatter(x_failed, y_failed, marker='x', color='red', label='Failed', s=100)

    # Color bar
    cbar = plt.colorbar(sc)
    cbar.set_label("mean(1 - IoU)", fontsize=18)

    # Set axes and labels
    ax1.set_xlabel("Surface Tension (N/m)", fontsize=20)
    ax1.set_xscale("log")
    ax1.set_ylabel("Adherence Strength (Pa/m)",fontsize=20)
    ax1.set_yscale("log")
    ax1.legend()

    # Set axes ticks
    # Add minor ticks to the axes
    x_minor_locator = LogLocator(base=10.0, subs=x)
    ax1.xaxis.set_minor_locator(x_minor_locator)
    ax1.xaxis.set_minor_formatter(ScalarFormatter(useMathText=True))
    y_minor_locator = LogLocator(base=10.0, subs=y)
    ax1.yaxis.set_minor_locator(y_minor_locator)
    ax1.yaxis.set_minor_formatter(ScalarFormatter(useMathText=True))

    # Customize major ticks
    ax1.xaxis.set_major_locator(LogLocator(base=10.0))
    ax1.xaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax1.yaxis.set_major_locator(LogLocator(base=10.0))
    ax1.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))

    # Customize the appearance of ticks and labels
    ax1.tick_params(axis='x', which='minor', size=4, rotation=45)
    ax1.tick_params(axis='x', which='major', size=8, rotation=45)
    ax1.tick_params(axis='y', which='minor', size=4)
    ax1.tick_params(axis='y', which='major', size=8)


    # Add num_iterations plot
    ax2 = fig.add_subplot(1, 3, 2)
    sc = ax2.scatter(x, y, c=num_iterations, cmap='viridis', s=100)

    # Color bar
    cbar = plt.colorbar(sc)
    cbar.set_label("Number of iterations", fontsize=18)

    # Set axes and labels
    ax2.set_xlabel("Surface Tension (N/m)", fontsize=20)
    ax2.set_xscale("log")
    ax2.set_ylabel("Adherence Strength (Pa/m)",fontsize=20)
    ax2.set_yscale("log")

    # Set axes ticks
    # Add minor ticks to the axes
    x_minor_locator = LogLocator(base=10.0, subs=x)
    ax2.xaxis.set_minor_locator(x_minor_locator)
    ax2.xaxis.set_minor_formatter(ScalarFormatter(useMathText=True))
    y_minor_locator = LogLocator(base=10.0, subs=y)
    ax2.yaxis.set_minor_locator(y_minor_locator)
    ax2.yaxis.set_minor_formatter(ScalarFormatter(useMathText=True))

    # Customize major ticks
    ax2.xaxis.set_major_locator(LogLocator(base=10.0))
    ax2.xaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax2.yaxis.set_major_locator(LogLocator(base=10.0))
    ax2.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))

    # Customize the appearance of ticks and labels
    ax2.tick_params(axis='x', which='minor', size=4, rotation=45)
    ax2.tick_params(axis='x', which='major', size=8, rotation=45)
    ax2.tick_params(axis='y', which='minor', size=4)
    ax2.tick_params(axis='y', which='major', size=8)


    # Add simulation ids plot
    ax3 = fig.add_subplot(1, 3, 3)
    sim_ids = np.arange(1, len(iou_losses)+1)
    sim_ids = [str(i) for i in sim_ids]
    ax3.scatter(x, y, color='white')
    for i in range(len(x)):
        ax3.annotate(sim_ids[i], (x[i], y[i]))
    dummy_point = ax3.scatter([], [], marker='o', color='white', label="Simu ids")

    # Set axes and labels
    ax3.set_xlabel("Surface Tension (N/m)", fontsize=20)
    ax3.set_xscale("log")
    ax3.set_ylabel("Adherence Strength (Pa/m)",fontsize=20)
    ax3.set_yscale("log")
    ax3.legend(handles=[dummy_point])

    # Set axes ticks
    # Add minor ticks to the axes
    x_minor_locator = LogLocator(base=10.0, subs=x)
    ax3.xaxis.set_minor_locator(x_minor_locator)
    ax3.xaxis.set_minor_formatter(ScalarFormatter(useMathText=True))
    y_minor_locator = LogLocator(base=10.0, subs=y)
    ax3.yaxis.set_minor_locator(y_minor_locator)
    ax3.yaxis.set_minor_formatter(ScalarFormatter(useMathText=True))

    # Customize major ticks
    ax3.xaxis.set_major_locator(LogLocator(base=10.0))
    ax3.xaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax3.yaxis.set_major_locator(LogLocator(base=10.0))
    ax3.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))

    # Customize the appearance of ticks and labels
    ax3.tick_params(axis='x', which='minor', size=4, rotation=45)
    ax3.tick_params(axis='x', which='major', size=8, rotation=45)
    ax3.tick_params(axis='y', which='minor', size=4)
    ax3.tick_params(axis='y', which='major', size=8)

    # Save plot
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
        plt.savefig(os.path.join(out_dir, "mech_props_plot.png"))
    
    plt.show()

if __name__ == "__main__":
    
    assert argv.__len__() == 4, (
        "Wrong command line arguments: \n expected "
        "python3 mech_props_analysis.py [simulation_outputs_path] [screening_table_path] [output_path]"
    )
    simulation_outputs_path = argv[1]
    screening_table_path = argv[2]
    output_path = argv[3] 

    screen_params_dict = get_screening_parameters(screening_table_path)
    iou_loss_vals, num_iters = compute_start_vs_end_IoU_loss(simulation_outputs_path)
    plot_IoU_loss(
        iou_loss_vals,
        screen_params_dict,
        num_iters,
        output_path
    )
