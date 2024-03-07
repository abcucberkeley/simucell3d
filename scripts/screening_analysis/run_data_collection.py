import os 
import subprocess
from typing import List, Tuple

'''
Run this script to launch all the jobs that collect data from a parameter screening run.
'''

#------------------------------------------------------------------------------------------------------
def write_sbatch_command(
        script_name: str,
        args_lst: List[str],
        hours: int = 4,
        n_cpus: int = 4,
        cpu_mem: int = 4096
    ) -> str:
    """
    Write a sbatch command to run a python script on the cluster.

    Parameters:
    -----------
    script_name: (str)
        The name of the python script.
    args_lst: (List[str])
        A list of arguments for running the script, all in str format. 
    hours: (int = 4)
        Number of hours reserved for the job.
    n_cpus: (int = 4)
        Number of cpus reserved for the job.
    cpu_mem: (int = 4096)
        Memory per cpu reserved for the job.
    """

    args_string = ""
    for arg in args_lst:
        args_string += str(arg)
        args_string += " "
    args_string = args_string[:-1]

    sbatch_cmd = f"""\
    sbatch -n 1 --cpus-per-task={n_cpus} --time={hours}:00:00 --mem-per-cpu={cpu_mem} --wrap="python {script_name} {args_string}"
    """

    return sbatch_cmd.strip()
#------------------------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------------------------
def run_command(
        cmd: str
    )-> Tuple[int, str, str]:
    """
    Run a given command.

    Parameters:
    -----------
    cmd: (str)
        The command to run.

    Returns:
    --------
    (int): 
        The return code
    """
    try:
        result = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)

        # Return the standard output and standard error as a tuple
        return result.returncode

    except subprocess.CalledProcessError as e:
        # If the command returns a non-zero exit code, handle the exception
        return e.returncode
#------------------------------------------------------------------------------------------------------



#------------------------------------------------------------------------------------------------------
def print_job_output(
    returncode: int,
) -> None:

    if returncode == 0:
        print("Command executed successfully!")
    else:
        print(f"Command failed with return code {returncode}.")
#------------------------------------------------------------------------------------------------------


if __name__ == "__main__":
    
    ##### USER SECTION #####
    ### SPECIFY PATH TO SCREENING ROOT DIRECTORY ###
    root_dir = "path/to/screening/root/directory/on/cluster"
    screening_csv_file_name = "screening_table.csv"
    ################################################

    # Set paths and create output folder
    simulation_outputs_dir = os.path.join(root_dir, "simulation_outputs")
    simulation_summary_dir = os.path.join(root_dir, "summary_folder", screening_csv_file_name)
    output_dir = root_dir
    output_files_dir = os.path.join(output_dir, 'dashboard_files')
    if not os.path.exists(output_files_dir):
        os.makedirs(output_files_dir)

    assert os.path.exists(simulation_outputs_dir), f"""\
        The provided path to 'simulation_outputs' directory does not exist {simulation_outputs_dir}
        """
    
    assert os.path.exists(simulation_summary_dir), f"""\
        The provided path to 'summary_folder/screening_results.csv' does not exist {simulation_summary_dir}
        """
    
    # SPHERICITY
    print('------------------------------------------------------------')
    print('Computing SIMULATION STATISTICS:')
    script = './utils/compute_statistics.py'
    args = [
        simulation_outputs_dir,
        output_dir
    ]
    command = write_sbatch_command(script, args)
    print(command)
    ret_code = run_command(command)
    print_job_output(ret_code)

    print('\n------------------------------------------------------------')
    print('Computing IOU:')
    script = './utils/compute_IoU_euler.py'
    args = [
        simulation_outputs_dir,
        output_dir
    ]
    command = write_sbatch_command(script, args)
    print(command)
    ret_code = run_command(command)
    print_job_output(ret_code)


    print('\n------------------------------------------------------------')
    print('Computing INTERPENETRATION:')
    script = './utils/compute_interpenetration_euler.py'
    args = [
        simulation_outputs_dir,
        output_dir
    ]
    command = write_sbatch_command(script, args)
    print(command)
    ret_code = run_command(command)
    print_job_output(ret_code)


    print('\n------------------------------------------------------------')
    print('Computing NUM ITERATIONS:')
    script = './utils/compute_num_iterations.py'
    args = [
        simulation_summary_dir,
        output_dir
    ]
    command = write_sbatch_command(script, args)
    print(command)
    ret_code = run_command(command)
    print_job_output(ret_code)
    