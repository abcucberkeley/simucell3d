from os import getcwd, path, listdir, mkdir, _exit, popen




path_to_output_folder = r"/cluster/scratch/srunser/cell_triplet_test"


sim_folder_lst = listdir(path_to_output_folder)
summary_folder_path = path.join(path_to_output_folder, "summary_cell_triplet") 
if not path.exists(summary_folder_path): mkdir(summary_folder_path)

for sim_folder in sim_folder_lst:
    if not sim_folder.startswith("simulation_"): continue
    sim_id = sim_folder.split("_")[1]
    sim_folder_path = path.join(path_to_output_folder, sim_folder, sim_folder, "face_data")
    sim_file_lst = listdir(sim_folder_path)

    if "result_100.vtk" in sim_file_lst:
        popen("cp " + path.join(sim_folder_path, "result_100.vtk") + " " + path.join(summary_folder_path, "simulation_" + sim_id + ".vtk"))



