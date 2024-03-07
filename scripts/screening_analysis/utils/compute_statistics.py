import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import json
from sys import argv
from tqdm import tqdm
from typing import Iterable, Tuple, Dict, List, Literal


#----------------------------------------------------------------------------------------
def _extract_stats_from_csv(
        path_to_csv: str,
        features: Iterable[str]
    ) -> Dict[str, np.ndarray[float]]:
    '''
    Load simulation_statistics.csv dataframe and compute the value of the
    cell sphericities for the initial and final iterations of the simulation run.

    Parameters:
    -----------
        path_to_csv: (str)
            The path to simulation_statistics.csv dataframe

        features: (Iterable[str])
            A collection of feature names to extract from the statistics dataframe.

    Returns:
    --------
        out_dict: (Dict[str, np.ndarray[float])
            A dictionary containing arrays of selected statistics for the initial and final
            geometry of the given simulation run.
    '''

    stats_df = pd.read_csv(path_to_csv)

    iter_ids = np.unique(stats_df['iteration'])

    out_dict = {}

    init_iter_df = stats_df[stats_df['iteration'] == iter_ids[0]]
    final_iter_df = stats_df[stats_df['iteration'] == iter_ids[-1]]
    for feature in features:
        out_dict['init_' + feature] = init_iter_df[feature].to_numpy()
        out_dict['final_' + feature] = final_iter_df[feature].to_numpy()
    
    return out_dict 
#----------------------------------------------------------------------------------------



#----------------------------------------------------------------------------------------
def get_screening_stats(
        path_to_sim_folder: str,
        feature_names: Iterable[str]
) -> Dict[str, Dict[str, np.ndarray[float]]]:
    '''
    Parse directories of simulation_outputs in the screening folder.
    Retrieve simulation_statistics.csv dataframes and for each of them extract values of 
    selected statistics for the initial and final iteration of the simulation run.

    Parameters:
    -----------
        path_to_sim_folder: (str)
            The path to `simulation_outputs` directory in the screening results folder.

        feature_names: (Iterable[str])
            A collection of feature names to extract from the statistics dataframes.

    Returns:
    --------
        stats_dict: (Dict[str, Dict[str, np.ndarray[float]]])
            A dictionary whose keys are simulation ids as strings, and values are dictionary
            of feature names and arrays of data.   
    '''

    assert os.path.exists(path_to_sim_folder), """\
        You should provide the path to the `simulation_output` directory in the screening folder."""
    
    sim_subdirs = os.listdir(path_to_sim_folder)
    stats_dict = {}
    for sim_subdir in tqdm(sim_subdirs, desc='Parsing simulation runs'):
        sim_id = sim_subdir.split('_')[1]
        curr_stats_dict = _extract_stats_from_csv(
            os.path.join(path_to_sim_folder, sim_subdir, f'simulation_statistics.csv'),
            feature_names
        )
        stats_dict[sim_id] = curr_stats_dict

    return stats_dict
#----------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------
def _compute_sphericity(
        simulation_stats_dict: Dict[str, np.ndarray[float]]
    ) -> np.ndarray[float]:
    '''
    Given the statistics for a given simulation run, compute a list of delta sphericity values.

    Parameters:
    -----------
        simulation_stats_dict: (Dict[str, np.ndarray[float]])
            A dictionary storing the selected statistics for a given simulation run.

    Returns:
    --------
        (np.ndarray[float])
            An array storing the difference between final and initial sphericity values.
    '''

    init_areas = simulation_stats_dict['init_area']
    final_areas = simulation_stats_dict['final_area']
    init_volumes = simulation_stats_dict['init_volume']
    final_volumes = simulation_stats_dict['final_volume']
    
    def sphericity_formula(area, volume):
        return (np.pi**(1/3)) * (6 * volume)**(2/3) / area
    
    init_sph = sphericity_formula(init_areas, init_volumes)
    final_sph = sphericity_formula(final_areas, final_volumes)

    if len(init_sph) == len(final_sph):
        return (final_sph - init_sph) / init_sph * 100
    else:
        return np.nan
#----------------------------------------------------------------------------------------



#----------------------------------------------------------------------------------------
def _compute_volume_loss(
        simulation_stats_dict: Dict[str, np.ndarray[float]]
    ) -> np.ndarray[float]:
    '''
    Given the statistics for a given simulation run, compute the volume loss.

    Parameters:
    -----------
        simulation_stats_dict: (Dict[str, np.ndarray[float]])
            A dictionary storing the selected statistics for a given simulation run.

    Returns:
    --------
        volume_losses: (np.ndarray[float])
            List of sphericity values, each value associated to a different cell.
    '''

    init_volumes = simulation_stats_dict['init_volume']
    final_volumes = simulation_stats_dict['final_volume']
    if len(init_volumes) == len(final_volumes):
        volume_losses = (np.asarray(final_volumes) - np.asanyarray(init_volumes)) / np.asanyarray(init_volumes) * 100
        return volume_losses
    else:
        return np.nan
#----------------------------------------------------------------------------------------



#----------------------------------------------------------------------------------------
def collect_statistics(
    screening_stats_dict: Dict[str, Dict[str, np.ndarray[float]]],
    to_compute: List[Literal['delta_sphericity', 'volume_loss']],
    path_to_output_dir: str
) -> None:
    '''
    Given the statistics dictionary computed before, compute the selected statistics and write
    the results on file.

    Parameters:
    -----------
    screening_stats_dict: (Dict[str, Dict[str, np.ndarray[float]]])
        A dictionary whose keys are simulation ids as strings, and values are dictionary
        of feature names and arrays of data. 

    to_compute: (List[Literal['delta_sphericity', 'volume_loss']])
        A collection of statistics names to calculate from the statistics dataframe.

    path_to_output_dir: (str)
        The path to the directory in which results are saved.
    '''
    
    out_dict = {}

    avail = ['delta_sphericity', 'volume_loss']
    for name in to_compute:
        assert name in avail, f"Cannot compute {name}! Pick one from {avail}."
        out_dict[name] = {}

    # Compute statistics
    for sim_id in tqdm(screening_stats_dict.keys()):
        if 'delta_sphericity' in to_compute:
            res = _compute_sphericity(screening_stats_dict[sim_id])
            out_dict['delta_sphericity'][sim_id] = (np.mean(res), np.std(res))
    
        if 'volume_loss' in to_compute:
            res = _compute_volume_loss(screening_stats_dict[sim_id])
            out_dict['volume_loss'][sim_id] = (np.mean(res), np.std(res))

    # Write statistics to file
    for stat in to_compute:
        output_file_path = os.path.join(path_to_output_dir, f'{stat}_output.json')
        with open(output_file_path, 'w') as f_out:
            json.dump(out_dict[stat], f_out, indent=4)

    return out_dict
#----------------------------------------------------------------------------------------



#----------------------------------------------------------------------------------------
def main(
    simulation_dir: str, 
    output_dir: str,
    statistics: Iterable[str]
) -> None:
    '''
    Calls all the previous functions to collect the statistics.

    Parameters:
    -----------
    simulation_dir: (str) 
        The path to the `simulation_outputs` folder in the screening output directory.

    output_dir: (str)
        The path to the directory that will store thhe file in which the output is written.

    statistics: (Iterable[str])
        A collection of statistics names to calculate from the statistics dataframe.
    '''

    assert os.path.exists(output_dir), f"No such directory find at {output_dir}."

    stats_dict = get_screening_stats(
        path_to_sim_folder=simulation_dir, 
        feature_names=['area', 'volume']
    )

    collect_statistics(
        screening_stats_dict=stats_dict, 
        to_compute=statistics, 
        path_to_output_dir=output_dir
    )
#----------------------------------------------------------------------------------------


if __name__ == '__main__':

    print('COMPUTING CELL STATISTICS')

    assert argv.__len__() == 3 , """Wrong command line arguments: \n expected python compute_statistics.py [path_to_simulation_stats] [output_dir]
        Example: python compute_statistics.py /cluster/scratch/fcarrara/bronchiole_screening/simu_length_screening_v8/simulation_outputs
        /cluster/scratch/fcarrara/bronchiole_screening/simu_length_screening_v8/"""
    sim_out_dir = argv[1]
    output_dir = argv[2]

    main(
        simulation_dir=sim_out_dir, 
        output_dir=os.path.join(output_dir, 'dashboard_files'),
        statistics=['delta_sphericity', 'volume_loss']
    )