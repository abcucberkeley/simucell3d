import os 
import subprocess
import shutil
from sys import argv
from typing import List, Union

"""
Automatically launches different instances of `check_mesh_collision.py` in parallel on Euler.

The user is required to specify:
    - `indexes`: the simulation ids to be evaluated.
    - `simulation_dir`: the path to the `simulation_outputs` folder in the screening output directory.
    - `output_file_path`: the path to the file in which to write the output.
"""

#-----------------------------------------------------------------------------------------------------
def collect_collisions_from_screening(
    simulation_dir: str,
    output_dir: str,
    indexes: Union[List[int], 'all'],    
) -> None:
    '''
    Automatically launch different instances of `check_mesh_collision.py` in parallel on Euler.

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
    
    output_file_path = os.path.join(output_dir, 'dashboard_files/interpenetration_output.txt')

    # Create an empty output file
    with open(output_file_path, "w") as file:
        pass

    # Check that indexes is consistent with number of simulation for this screening
    if indexes == 'all' or len(indexes) > len(os.listdir(simulation_dir)):
        indexes = list(range(1, len(os.listdir(simulation_dir))+1))


    # Generate job scripts for different simulation ids
    if not os.path.exists('./bash_scripts_collision'):
        os.makedirs('./bash_scripts_collision')

    for idx in indexes:
        simu_path = os.path.join(simulation_dir, f'simulation_{idx}/cell_data')

        slurm_script = f"""\
#!/bin/bash

#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=8196

python ./utils/compute_interpenetration.py {simu_path} {output_file_path}
        """
        slurm_script = slurm_script.strip()
        path_to_job_script = f'./bash_scripts_collision/script_{idx}.sh'

        # Write the Bash script to a file
        with open(path_to_job_script, 'w') as f:
            f.write(slurm_script)

        # Execute the slurm script
        subprocess.run(['sbatch', path_to_job_script], stdin=subprocess.PIPE)

    # Delete the jobs folder
    shutil.rmtree('./bash_scripts_collision')
#-----------------------------------------------------------------------------------------



if __name__ == '__main__':

    print('CHECKING CELL INTERPENERATION \n')

    assert argv.__len__() == 3, "Wrong command line arguments: \n expected python compute_interpenetration_euler.py [simulation_dir] [output_dir]"+\
        "Example: python compute_interpenetration_euler.py /cluster/scratch/fcarrara/bronchiole_screening/simu_length_screening_v8/simulation_outputs/"+\
        "/cluster/scratch/fcarrara/bronchiole_screening/simu_length_screening_v8/"
    root_dir = argv[1]
    save_dir = argv[2]
    idxs = 'all'

    collect_collisions_from_screening(
        simulation_dir=root_dir,
        output_dir=save_dir,
        indexes=idxs
    )