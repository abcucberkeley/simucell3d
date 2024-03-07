from os import path, _exit
import pandas as pd 

path_to_script_folder = path.dirname(path.realpath(__file__))

#Get the path to the simucell3d python library
path_to_build = path.join(path_to_script_folder, "..", "..", "..", "..", "build")

#Get the path where the simulation results are stored
path_to_simulation_results = path.join(path_to_build, "test_cell_internalization_output", "simulation_0", "cell_statistics_data.csv")
sim_stats_df = pd.read_csv(path_to_simulation_results)

print(sim_stats_df[["cell_id", "volume", "area", "cell_contact_area_fraction", "pressure"]])





