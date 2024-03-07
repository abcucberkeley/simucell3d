import os 
import subprocess
import shutil
from sys import argv
from typing import List, Union

"""
Automatically launches different instances of `shape_change.py` in parallel on Euler.

The user is required to specify:
    - `idxs`: the simulation ids to be evaluated.
    - `root_dir`: the path to the `simulation_outputs` folder in the screening output directory.
    - `save_dir`: the path to the folder in which to save the resulting plots.
"""

#-----------------------------------------------------------------------------------------------------
def collect_shape_change_from_screening(
    simulation_dir: str,
    output_dir: str,
    indexes: Union[List[int], 'all'],    
) -> None:
    '''
    Automatically launch different instances of `compute_IoU.py` in parallel on Euler.

    Parameters:
    -----------
    simulation_dir: (str) 
        The path to the `simulation_outputs` folder in the screening output directory.

    output_dir: (str)
        The path to the directory that will store the file in which the output is written.

    indexes: (List[int]) 
        The simulation ids to evaluate. If 'all' all the indexes are taken.        
    '''

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    output_file_path_1 = os.path.join(output_dir, 'dashboard_files/IoU_derivative_output.txt')
    output_file_path_2 = os.path.join(output_dir, 'dashboard_files/IoU_output.txt')

    # Create an empty output file
    with open(output_file_path_1, "w") as file:
        pass

    with open(output_file_path_2, "w") as file:
        pass

    # Check that idxs is consistent with number of simulation for this screening
    if indexes == 'all' or len(idxs) > len(os.listdir(simulation_dir)):
        idxs = list(range(1, len(os.listdir(simulation_dir))+1))

    if not os.path.exists('./bash_scripts_shape'):
        os.makedirs('./bash_scripts_shape')

    for idx in idxs:

        simu_path = os.path.join(simulation_dir, f'simulation_{idx}/cell_data')

        slurm_script = f"""\
#!/bin/bash

#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=8196

python ./utils/compute_IoU.py {simu_path} {idx} {output_dir}
        """
        path_to_job_script = f'./bash_scripts_shape/script_{idx}.sh'

        # Write the Bash script to a file
        with open(path_to_job_script, 'w') as f:
            f.write(slurm_script)

        # Execute the slurm script
        subprocess.run(['sbatch', path_to_job_script], stdin=subprocess.PIPE)

    # Delete the jobs folder
    shutil.rmtree('./bash_scripts_shape')
#-----------------------------------------------------------------------------------------------------



if __name__ == '__main__':

    print('COMPUTING SHAPE CHANGE \n')
    
    assert argv.__len__() == 3 , "Wrong command line arguments: \n expected python compute_IoU_euler.py [simulation_dir] [output_dir]"+\
        "Example: compute_IoU_euler.py /cluster/scratch/fcarrara/bronchiole_screening/simu_length_screening_v8/simulation_outputs/"+\
        "/cluster/scratch/fcarrara/bronchiole_screening/simu_length_screening_v8/"
    root_dir = argv[1]
    save_dir = argv[2]
    idxs = 'all'

    collect_shape_change_from_screening(
        simulation_dir=root_dir,
        output_dir=save_dir,
        indexes=idxs
    )