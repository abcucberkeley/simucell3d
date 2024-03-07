import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import json
from sys import argv
from tqdm import tqdm
from typing import Iterable, Tuple, Dict

'''
Run this script to get plots of the sphericity of the cells at the end of a simulation run.

N.B. Sphericity is defined here: https://en.wikipedia.org/wiki/Sphericity
'''

#----------------------------------------------------------------------------------------
def compute_sphericity(
        areas: Iterable[float],
        volumes: Iterable[float],
    ) -> np.ndarray[float]:
    '''
    Given lists of surface areas and volumes for different cells, compute a list of
    sphericity values.

    Parameters:
    -----------
        areas: (Iterable[float])
            List of surface areas, each value associated to a different cell.

        volumes: (Iterable[float])
            List of volumes, each value associated to a different cell.

    Returns:
    --------
        sphericities: (np.ndarray[float])
            List of sphericity values, each value associated to a different cell.
    '''

    areas = np.asarray(areas)
    volumes = np.asarray(volumes)

    sphericities = (np.pi**(1/3)) * (6 * volumes)**(2/3) / areas

    return sphericities
#----------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------
def _extract_stats_from_csv(
        path_to_csv: str
    ) -> Tuple[np.ndarray[float], np.ndarray[float]]:
    '''
    Load simulation_statistics.csv dataframe and compute the value of the
    cell sphericities for the initial and final iterations of the simulation run.

    Parameters:
    -----------
        path_to_csv: (str)
            The path to simulation_statistics.csv dataframe

    Returns:
    --------
        (Tuple[np.ndarray[float], np.ndarray[float]])
            A tuple of arrays of sphericity values for the first and last iterations 
            of the simulation run.
        
    '''

    stats_df = pd.read_csv(path_to_csv)

    iter_ids = np.unique(stats_df['iteration'])

    init_iter_df = stats_df[stats_df['iteration'] == iter_ids[0]]
    init_areas = init_iter_df['area']
    init_volumes = init_iter_df['volume']

    final_iter_df = stats_df[stats_df['iteration'] == iter_ids[-1]]
    final_areas = final_iter_df['area']
    final_volumes = final_iter_df['volume']

    return (
        compute_sphericity(init_areas, init_volumes),
        compute_sphericity(final_areas, final_volumes),
    )
#----------------------------------------------------------------------------------------



#----------------------------------------------------------------------------------------
def get_screening_sphericity(
        path_to_sim_folder: str,
) -> Tuple[Dict[int, np.ndarray[float]], Dict[int, np.ndarray[float]]]:
    '''
    Parse directories of simulation_outputs in the screening folder.
    Retrieve simulation_statistics.csv dataframes and for each of them compute the value of the
    cell sphericities for the final iteration of the simulation run.

    Parameters:
    -----------
        path_to_sim_folder: (str)
            The path to `simulation_outputs` directory in the screening results folder.

    Returns:
    --------
        delta_sphericity_dict: (Dict[int, np.ndarray[float]])
            An dictionary that associates to each simulation id an array of delta sphericity values 
            between the first and last iterations of the simulation itself. 

        final_sphericity_dict: (Dict[int, np.ndarray[float]])
            An dictionary that associates to each simulation id an array of sphericity values 
            for the last iteration of the simulation itself. 
    '''

    assert os.path.exists(path_to_sim_folder), "You should provide the path to the `simulation_output` directory"+\
        "in the screening folder." 
    
    sim_subdirs = os.listdir(path_to_sim_folder)
    delta_sphericity_dict = {}
    final_sphericity_dict = {}
    for sim_subdir in tqdm(sim_subdirs, desc='Parsing simulation runs'):
        sim_id = sim_subdir.split('_')[1]
        res = _extract_stats_from_csv(
            os.path.join(path_to_sim_folder, sim_subdir, 'simulation_statistics.csv')
        )
        final_sphericity_dict[int(sim_id)] = res[1]
        delta_sphericity_dict[int(sim_id)] = res[1] - res[0]

    delta_sphericity_dict = dict(sorted(delta_sphericity_dict.items(), key=lambda x: x[0]))
    final_sphericity_dict = dict(sorted(final_sphericity_dict.items(), key=lambda x: x[0]))

    return delta_sphericity_dict, final_sphericity_dict
#----------------------------------------------------------------------------------------



#----------------------------------------------------------------------------------------
def collect_sphericity_from_screening(
       simulation_dir: str, 
       output_dir: str
) -> str:
    '''
    Wrapper for `get_screening_sphericity`. The latter is called and the results are saved
    in a dictionary, to be then used in the screening dashboard. 

    Parameters:
    -----------
    simulation_dir: (str) 
        The path to the `simulation_outputs` folder in the screening output directory.

    output_dir: (str)
        The path to the directory that will store thhe file in which the output is written.
    
    Returns:
    --------
    output_file_path: (str)
        The path to the file on which the output is written.
    '''

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    output_file_path = os.path.join(output_dir, 'sphericity_output.json')

    delta_dict, final_dict = get_screening_sphericity(simulation_dir)

    out_dict = {}
    for key in delta_dict.keys():
        out_dict[key] = (
            np.mean(delta_dict[key]), 
            np.std(delta_dict[key]),
            np.mean(final_dict[key]),
            np.std(final_dict[key])
        )

    with open(output_file_path, 'w') as f_out:
        json.dump(out_dict, f_out, indent=4)

    return output_file_path
#----------------------------------------------------------------------------------------



#----------------------------------------------------------------------------------------
def plot_sphericity(
    delta_sphericity_dict: Dict[int, np.ndarray[float]],
    final_sphericity_dict: Dict[int, np.ndarray[float]],
    save_dir: str = './plots'
) -> None:
    '''
    Plot cell sphericities for the final iteration of every simulation run in a screening.

    Parameters:
    -----------
        delta_sphericity_dict: (Dict[int, np.ndarray[float]])
            An dictionary that associates to each simulation id an array of delta sphericity values 
            between the first and last iterations of the simulation itself. 

        final_sphericity_dict: (Dict[int, np.ndarray[float]])
            An dictionary that associates to each simulation id an array of sphericity values 
            for the last iteration of the simulation itself.    

        save_dir: (str, default='./plots')     
            Where to save the produced plot.
    '''
    
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(24, 6))
    ax1, ax2 = axes
    fig.suptitle('Sphericity of cells in simulated samples', fontsize=30)

    #DELTA PLOT
    x = list(delta_sphericity_dict.keys())
    y = list(delta_sphericity_dict.values())
    y_mean = [np.mean(val) for val in y]
    y_ci = [1.96 * np.std(val)/np.sqrt(len(val)) for val in y]

    ax1.errorbar(
        x, y_mean, yerr=y_ci, fmt='o--',
        capsize=5, color='orange', ecolor='grey',
        label=r'Mean $\Delta$Sphericity w\ 95% CI'
    )
    ax1.set_xlabel("Simulation ID", fontsize=18)
    ax1.set_xticks(range(1, len(x)+1, 5), x[::5])
    ax1.set_ylabel(r"Mean $\Delta$Sphericity", fontsize=18)
    ax1.set_title(r"Mean $\Delta$Sphericity btw initial and final geometry", fontsize=24)
    ax1.legend(fontsize=14, loc='lower right')

    #FINAL PLOT
    x = list(final_sphericity_dict.keys())
    y = list(final_sphericity_dict.values())
    y_mean = [np.mean(val) for val in y]
    y_ci = [1.96 * np.std(val)/np.sqrt(len(val)) for val in y]

    ax2.errorbar(
        x, y_mean, yerr=y_ci, fmt='o--', 
        capsize=5, color='blue', ecolor='grey',
        label='Mean final sphericity w\ 95% CI'
    )
    ax2.plot(x, [1]*len(x), linestyle=':', color='red', label='Maximum sphericity achievable')
    ax2.set_xlabel("Simulation ID", fontsize=18)
    ax2.set_xticks(range(1, len(x)+1, 5), x[::5])
    ax2.set_ylabel("Mean Sphericity", fontsize=18)
    ax2.set_title("Mean Sphericity at final geometry", fontsize=24)
    ax2.legend(fontsize=14, loc='lower right')

    plt.tight_layout()

    fig.savefig(os.path.join(save_dir, 'sphericity_plot.png'))
#----------------------------------------------------------------------------------------



if __name__ == '__main__':

    print('COMPUTING CELL SPHERICITY \n')

    assert argv.__len__() == 3 , "Wrong command line arguments: \n expected python check_sphericity.py [path_to_simulation_outputs] [output_dir]"+\
        "Example: python check_sphericity.py /cluster/scratch/fcarrara/bronchiole_screening/simu_length_screening_v8/simulation_outputs/"+\
        "/cluster/scratch/fcarrara/bronchiole_screening/simu_length_screening_v8/"
    sim_out_dir = argv[1]
    output_dir = argv[2]

    # delta_sph_dict, final_sph_dict = get_screening_sphericity(sim_out_dir)

    # plot_sphericity(delta_sph_dict, final_sph_dict, save_dir=os.path.join(output_dir, 'screening_plots'))

    collect_sphericity_from_screening(
        simulation_dir=sim_out_dir, 
        output_dir=os.path.join(output_dir, 'dashboard_files')
    )