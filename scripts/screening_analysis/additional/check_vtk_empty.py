import os
from typing import List
import json
from tqdm import tqdm

'''
When launching screening on Euler, it may happens that some of the vtk files tracking the simulation
progress are not saved correctly. 
Therefore, this script serves for checking if there are any empty vtks in the outputs of a 
specific screening.
'''

def get_empty_files(
        folder_path, 
        size_threshold
    ) -> List[str]:
    """
    Analyzes files in a folder and returns the names of files whose size is under a given threshold.

    Parameters:
        folder_path (str): The path to the folder to be analyzed.
        size_threshold (int): The size threshold in bytes.

    Returns:
        list: A list of file names whose size is under the specified threshold.
    """

    if not os.path.exists(folder_path):
        print(f'Input path {folder_path} does not exist')
        return None
    
    under_threshold_files = []

    file_names = sorted(os.listdir(folder_path), key=lambda x: int(x.split('.')[0].split('_')[1]))
    for file_name in file_names:
        file_path = os.path.join(folder_path, file_name)
        if os.path.isfile(file_path) and os.path.getsize(file_path) < size_threshold:
            under_threshold_files.append(file_name)

    return under_threshold_files


if __name__ == "__main__":
    root_path = "/cluster/scratch/fcarrara/bronchiole_screening/simu_length_screening_v9/simulation_outputs"
    size_threshold = 1024  # Specify the threshold in bytes (e.g., 1024 bytes = 1KB)

    empty_files = {}
    subdirs = sorted(os.listdir(root_path), key=lambda x: int(x.split('_')[1]))
    for subdir in tqdm(subdirs, desc='Analyzing folders'):
        path_to_dir = os.path.join(root_path, subdir, 'cell_data')
        empty_files[os.path.basename(subdir)] = get_empty_files(path_to_dir, size_threshold)

    with open('empty_vtk_files.json', 'w') as file:
        json.dump(empty_files, file, indent=4)