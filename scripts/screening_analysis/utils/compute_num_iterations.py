import os
import numpy as np
import pandas as pd
import json
from sys import argv
from tqdm import tqdm
from typing import List, Dict


#-------------------------------------------------------------------------------------------------------
def read_num_iterations_from_csv(
    path_to_summary_table: str,
    path_to_output_dir: str
) -> None:
    '''
    Read the parameter table in `summary_folder` and extract the file counts and the
    information about the crashed simulation runs.

    Parameters:
    -----------
    path_to_summary_table: (str)
        The path to the `simulation_name.csv` table from which data are extracted.

    path_to_output_dir: (str)
        The path to the directory where to save the result.
    '''

    assert os.path.exists(path_to_output_dir), f"No such directory find at {path_to_output_dir}."

    summary_df = pd.read_csv(path_to_summary_table)
    sim_ids = summary_df['sim_id'].tolist()
    file_counts = summary_df['file_count'].tolist()
    has_crashed = summary_df['has_crashed'].tolist()

    out_dict = {}
    for i, sim_id in enumerate(sim_ids):
        out_dict[sim_id] = (file_counts[i], has_crashed[i])

    output_file_path = os.path.join(path_to_output_dir, f'dashboard_files/num_iterations_output.json')
    with open(output_file_path, 'w') as f_out:
        json.dump(out_dict, f_out, indent=4)
#-------------------------------------------------------------------------------------------------------
    

if __name__ == '__main__':

    print('COMPUTING NUM ITERATIONS')

    assert argv.__len__() == 3 , """Wrong command line arguments: \n expected python compute_num_iterations.py [path_to_summary_table] [output_dir]
        Example: 
        python compute_num_iterations.py /cluster/scratch/fcarrara/bronchiole_screening/simu_length_screening_v8/summary_folder/simu_length_screening_v8.csv
        /cluster/scratch/fcarrara/bronchiole_screening/simu_length_screening_v8/"""
    sim_summary_dir = argv[1]
    output_dir = argv[2]

    read_num_iterations_from_csv(sim_summary_dir, output_dir)
