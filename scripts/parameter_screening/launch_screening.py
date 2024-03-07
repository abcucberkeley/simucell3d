from os import _exit, path, getcwd, chdir, system, sep, getlogin, mkdir
from sys import argv
import numpy as np
import pandas as pd
from shutil import copyfile
from parameter_file_model import *


"""
---------------------------------------------------------------------------------------------------------------------------

    This scripts is designed to launch parameter screenings on the Euler cluster. 

    Step 1: read different inputs from the command line:

        python3 launch_screening [path_to_parameters_table] [name_of_the_screen]  

    Step 2: Create a folder on the scratch with the name of the screen

    Step 3: Generate all the parameters_[sim_id].xml files in the folder on the scratch

    Step 4: Create a folder in the project root with the name of the screen
        ~/SimuCell3D/screen_name_build

    Step 5:
        Build in this folder the program with the desired configuration

    Step 6: Create a job array with the code created in this build folder

---------------------------------------------------------------------------------------------------------------------------
"""



def generate_parameter_files(simulation_folder_path, parameter_folder_path, parameter_table_df):
    """
    Generate all the parameter XML files used to parametrize the simulations
    

    Attributes
    ----------
    parameterTable (pd.DataFrame) : Each row of the table will be used to create a parameter XML file 

    Returns
    -------
    None

    """

    #First create all the parameter files
    for i, parameter_serie in parameter_table_df.iterrows():

        parameter_serie_dict = parameter_serie.to_dict()
        parameter_serie_dict["output_mesh_folder_path"] = path.join(simulation_folder_path, "simulation_{}".format(i+1))

        #Update the xml Parameter file
        file_template = get_xml_file_model(parameter_serie_dict)

        #Write the parameter file
        with open(path.join(parameter_folder_path, "parameters_{}.xml".format(i+1)), "w") as mfile: 
            mfile.write(file_template)







if __name__ == "__main__":



    #Some stuff to manage the cwd and file transfer
    project_root_dir = path.join(getcwd(), "..",  "..")
    
    #Check that the good number of arguments has beenn passed in the cmd
    assert argv.__len__() == 4, "Wrong command line arguments: \n python3 launchscreening "+\
                                "[absolute_path_to_parameters_table.csv] [name_of_the_screen] [output_folder_path]\n"+\
                                "Example: python3 launchscreening /home/srunser/parameters_table.csv screen1 /cluster/scratch/srunser"

    parameters_table_path =     argv[1].strip()
    screen_name =               argv[2].strip()
    output_folder_path =        argv[3].strip()
    output_folder_path =       path.join(output_folder_path, screen_name)
    print(output_folder_path)

    #Check that the parameters_table_path points to a correct file
    assert path.exists(parameters_table_path), "The given parameters table path does not exist: " + parameters_table_path

    if not path.exists(output_folder_path): mkdir(output_folder_path)

    #Check that the parameter table is a csv file
    assert parameters_table_path[-4:] == ".csv", "The paramater table should be a csv file and not a {} file"

    #Open the parameter table
    parameter_table_df = pd.read_csv(parameters_table_path, delimiter=",")

    #Check that there aare multiple columns in parameter_table_df, it's a quick proxy to check that the correct delimiter was used
    assert parameter_table_df.shape[1] > 1, """A wrong delimiter was used in the csv file {}, the correct delimiter is: "," """.format(parameters_table_path)

    #Check that the name of the screen does not contain whitespaces or other forbiddent characters
    forbidden_char_set = set([sep, " ", ",", ";", "\"", "\\"])
    forbidden_char_intersection = set([*screen_name]).intersection(forbidden_char_set)
    assert len(forbidden_char_intersection) == 0, "Forbidden characters were found in the screen name :" + " ".join(list(forbidden_char_intersection))

    #Store in this folder the parameter files
    parameter_folder = path.join(output_folder_path, "parameters") 
    mkdir(parameter_folder)

    #Store the simulation outputs in this directory
    simulation_folder = path.join(output_folder_path, "simulation_outputs") 
    mkdir(simulation_folder)
    
    #Store the stdout and stderr in this folder
    stdout_folder = path.join(output_folder_path, "std_out_err") 
    mkdir(stdout_folder)

    #Generate the XML parameters files based on the parameter table given in input
    generate_parameter_files(simulation_folder, parameter_folder, parameter_table_df)

    #Copy the parameters table in this directory on the scratch
    system("cp {} {}".format(parameters_table_path, output_folder_path))

    #Create a build folder with the name of the screen in the project root directory
    build_folder_path = path.join(project_root_dir, "build_{}".format(screen_name))
    assert not path.exists(build_folder_path), "The build folder {} already exists.".format(build_folder_path)
    mkdir(build_folder_path)
    chdir(build_folder_path)

    ##Compile the program in the build folder in Release mode
    system("""cmake -DCMAKE_BUILD_TYPE=Release .. ; make""")

    #Get the username of the person runningt the script
    user_name = getlogin().split("@")[0]

    #Create the command to launch the job array
    command_str   = """sbatch -n 1 --cpus-per-task=4 --time=48:00:00 --job-name="{screen_name}" --array=1-{nb_jobs} --mem-per-cpu=2048 --mail-type=END --mail-user={user_name}@bsse.ethz.ch --wrap="export OMP_NUM_THREADS=4; ./simucell3d {parameter_path}_\$SLURM_ARRAY_TASK_ID.xml &> {stdout_path}_\$SLURM_ARRAY_TASK_ID.txt" """.format(
        screen_name = screen_name,
        user_name = user_name,
        nb_jobs = parameter_table_df.shape[0],
        stdout_path = path.join(stdout_folder,"stdout"),
        parameter_path = path.join(parameter_folder,"parameters")
    )

    print(command_str)

    #Start the job array
    system(command_str)
    print("{} jobs have been launched on the cluster".format(parameter_table_df.shape[0]))



    



    





    

