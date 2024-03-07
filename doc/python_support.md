# Run SimuCell3D from python

Several python bindings have been implemented such that simulations can be calibrated, launched and analyzed from python. (Tested with python 3.10.6).

### Installation
Compile the code with the python bindinds enabled:
```
cd path/to/SimuCell3D/build
cmake -DENABLE_PYTHON_BINDINGS=TRUE -DPYTHON_EXECUTABLE=$(which python3) -DCMAKE_BUILD_TYPE=Release .. && make -j6
```

### How to start a simulation?
An example of how to start a simulation from python can be found in the file `scripts/python_bindings/launch_simulation_example.py`. 
Here is a quick summaryon how it works:

First you need to make sure that the path to the build folder is correct:
```
path_to_build = "/path/to/SimuCell3D/build"
```
Then you need to define the global simulation parameters (time step, etc):
```
global_simulation_parameters = simucell3d.global_simulation_parameters()
global_simulation_parameters.output_folder_path_ =      path.join(path_to_build, "simulation_results")
global_simulation_parameters.input_mesh_path_ =         "data/input_meshes/cube.vtk"
global_simulation_parameters.damping_coefficient_ =     5e-10
global_simulation_parameters.simulation_duration_ =     3e-06
global_simulation_parameters.sampling_period_ =         1e-06
global_simulation_parameters.time_step_ =               1e-07
global_simulation_parameters.min_edge_len_ =            4e-07
global_simulation_parameters.contact_cutoff_adhesion_ = 2e-07
global_simulation_parameters.contact_cutoff_repulsion_ = 2e-07

```

Then you can define the parameters of your face types:
```
face_type_1 = simucell3d.face_type_parameters()
face_type_1.name_                = "apical"     
face_type_1.face_type_global_id_ = 0
face_type_1.adherence_strength_  = 2.2e9
face_type_1.repulsion_strength_  = 1e9
face_type_1.surface_tension_     = 7e-4 
face_type_1.bending_modulus_     = 0. 

face_type_2 = simucell3d.face_type_parameters()
face_type_2.name_                = "lateral"     
face_type_2.face_type_global_id_ = 1
face_type_2.adherence_strength_  = 2.2e9
face_type_2.repulsion_strength_  = 1e9
face_type_2.surface_tension_     = 7e-4 
face_type_2.bending_modulus_     = 0. 
```

And finally you can define the parameters of your cell types:
```
cell_type_1 = simucell3d.cell_type_parameters()
cell_type_1.name_                       = "epithelial"     
cell_type_1.global_type_id_             = 0 
cell_type_1.mass_density_               = 1e3
cell_type_1.bulk_modulus_               = 2.5e3    
cell_type_1.max_pressure_               = 5e3  
cell_type_1.avg_growth_rate_            = 5e-14
cell_type_1.std_growth_rate_            = 0  
cell_type_1.target_isoperimetric_ratio_ = 150
cell_type_1.area_elasticity_modulus_    = 0   
cell_type_1.cortex_thickness_           = 1.9e-7   
cell_type_1.cortex_dynamic_viscosity_   = 1e-2   
cell_type_1.avg_division_vol_           = 3e-16    
cell_type_1.std_division_vol_           = 0   
cell_type_1.min_vol_                    = 5e-17   
cell_type_1.add_face_type(face_type_1) #Add the apical face type
cell_type_1.add_face_type(face_type_2) #Add the lateral face type
```

Now, the simulation can be launched by calling the python wrapper:
```
simulation_wrapper = simucell3d.simucell3d_wrapper(
    global_simulation_parameters, 
    [cell_type_1], 
    False #=VERBOSE
)
```

Once the simulation is finished, you can check if the simulation was successful:
```
simulation_outputs.RETURN_CODE_ == 0 #Means that it was successful
```

And you can convert the statistics collected on the cell to a pandas dataframe:
```
cell_stat_df = get_cell_statistics_data(simulation_outputs)
```

Note that you can add new statistics by adding new functions to the vector named `file_data_mapper_lst` of the file `include/io/mesh_data.hpp`.

