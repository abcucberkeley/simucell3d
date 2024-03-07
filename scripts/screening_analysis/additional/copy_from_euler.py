import subprocess
import os
from sys import argv

assert argv.__len__() == 3, "Wrong command line arguments: \n python3 copy_from_euler.py [path_to_screening_outputs_on_euler] [path_to_dest_dir_local]"
remote_outputs_path = argv[1].strip()
local_output_path = argv[2].strip()
remote_machine = "fcarrara@euler.ethz.ch"

# Get a list of subdirectories in the remote "outputs" directory
remote_subdirs = subprocess.check_output(f"ssh {remote_machine} 'ls {remote_outputs_path}'", shell=True)
remote_subdirs = remote_subdirs.decode("utf-8").split()

idxs_to_copy = list(range(1, len(remote_subdirs) + 1, 1))
# print(idxs_to_copy)

# Iterate through each remote subdirectory
for remote_subdir in remote_subdirs:
    sim_id = int(os.path.basename(remote_subdir).split("_")[1])
    if sim_id not in idxs_to_copy:
        continue 

    # Construct remote paths
    remote_cell_data_path = os.path.join(remote_outputs_path, remote_subdir, "cell_data")
    
    # Construct local paths
    local_output_subdir = os.path.join(local_output_path, remote_subdir)
    local_cell_data_path = os.path.join(local_output_subdir, "cell_data")
    
    # Create local output subdirectory if it doesn't exist
    os.makedirs(local_output_subdir, exist_ok=True)
    
    # Run scp command to copy "cell_data" from remote to local
    scp_command = f"scp -r {remote_machine}:{remote_cell_data_path} {local_cell_data_path}"
    subprocess.run(scp_command, shell=True)
