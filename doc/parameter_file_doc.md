# How to configure a simulation?

<br>

## Compilation time parameters

The contact model, time integration method, and the polarization technique can be chosen by modifying the preprocessor macros of the file `inlude/global_configuration.hpp`. To make the changes effective, the code has to be recompiled.

<br>
<br>


## Runtime parameters

The numerical parameters of the simulation (`time_step` etc) as well as the cell 
biomechanical parameters (`cell_bulk_modulus` etc) are stored in an XML file `parameters.xml`. This file has to be passed as an argument to the executable `simucell3d` : 

```
cd path/to/SimuCell3D/build
./simucell3d path/to/parameters.xml
```

The modifications of the parameters in the XML file are effective without recompiling the code. This file is structured in the following way:

<br>

The first section of the file contains the numerical parameters: 

```
<numerical_parameters>
    <!-- The mesh file containing the initial configuration <! -->
    <input_mesh_file_path>data/input_meshes/2_cubes.vtk</input_mesh_file_path> 

    <!-- The path where the mesh files will be written <! -->
    <output_mesh_folder_path>./simulation_results</output_mesh_folder_path> 

    <!-- The damping coefficients [M / T] <! -->
    <damping_coefficient>1e6</damping_coefficient> 

    ...

    <!-- The interaction cutoff distance above which cell-cell interactions are not considered anymore[L]<! -->
    <contact_cutoff_adhesion>2.5e-07</contact_cutoff_adhesion> 
    <contact_cutoff_repulsion>2.5e-07</contact_cutoff_repulsion> 

</numerical_parameters>
```

The path to the input mesh file is given in the XML tag <input_mesh_file_path>. This input mesh file has to follow a specific format which is is explained in more detail [here](./input_mesh_format.md). 
The `sampling_period` corresponds to the time interval between the saving of the mesh files. 
The `min_edge_length` tag corresponds to the minimum edge length allowed for the cell meshes. A small value gives small triangles and therefore a higher geometrical resolution, 
The 'contact_cutoff_adhesion' and 'contact_cutoff_repulsion' tags correspond to the distances between the faces of adjacent cells above which the contact forces are zero. 

<br> 

The second part of the XML file contains the mechanical parameters of each cell type:

```
<cell_types>

    <cell_type>
        [data of cell type 1]
    </cell_type


    <cell_type>
        [data of cell type 2]
    </cell_type

    ...

</cell_types>
```

The data of each cell type is stored in the following way: 

```
<!-- The name of the cell type <! -->
<cell_type_name>epithelial</cell_type_name> 

<!-- The ID of the cell type <! -->
<global_cell_id>0</global_cell_id>  

<!-- The cell density [M /(L^3)]<! -->
<cell_mass_density>1.0e3</cell_mass_density> 

<!-- The bulk modulus of the cells [M /(L T^2)] <! -->
<cell_bulk_modulus>6e3</cell_bulk_modulus> 

<!-- The maximal inner pressure of a cell. Set to INF to have no pressure limit [M / (T^2 L)] <! -->
<max_inner_pressure>INF</max_inner_pressure> 

<!-- Growth rate = cellDivisionTargetVolume / cellCycleDuration [L^3 / T] Default: 3e-12 <! -->
<avg_growth_rate>0</avg_growth_rate> 
<std_growth_rate>0</std_growth_rate> 

<!-- Isoperimetric ratio = cellArea^3 / cellVolume^2 <! -->
<target_isoperimetric_ratio>150</target_isoperimetric_ratio> 

<!-- The elastic modulus of the cell membranes [M / (T^2 L^2)] <! -->
<area_elasticity_modulus>2e-17</area_elasticity_modulus> 

<!-- The volume that the cell has to reach before division [L^3] (Drawn from a normal distribution)<! -->
<avg_division_volume>2e-16</avg_division_volume> 
<std_division_volume>0</std_division_volume> 

<!-- Growth rate = cellDivisionTargetVolume / cellCycleDuration [L^3 / T]<! -->
<!-- If the cell reach this volume we delete it[L^3]<! -->
<min_vol>1e-18</min_vol> 

<!-- The cell might have different face types with different mechanical properties <! -->
<face_types>     
    <face_type>           
        <global_face_id>0</global_face_id>

        <face_type_name>apical</face_type_name>

        <!-- Modulate the adherence between cells [(M / (T^2 L^2)] <! -->
        <adherence_strength>1e9</adherence_strength>

        <!-- Modulate the repulsion between cells [(M / (T^2 L^2)] <! -->
        <repulsion_strength>1e9</repulsion_strength> 

        <!-- The surface tension of the cell membranes [M / T^2] <! -->
        <surface_tension>5e-5</surface_tension> 

        <!-- The bending elasticity of the cell surface [M L^2 / T^2] <! -->
        <bending_modulus>5e-19</bending_modulus> 
    </face_type>

    <face_type>           
        <global_face_id>1</global_face_id>

        <face_type_name>lateral</face_type_name>

        <!-- Modulate the adherence between cells [M / (T^2 L^2)] <! -->
        <adherence_strength>2.2e9</adherence_strength>

        <!-- Modulate the repulsion between cells [(M / (T^2 L^2)] <! -->
        <repulsion_strength>1e9</repulsion_strength> 

        <!-- The surface tension of the cell membranes [M / T^2] <! -->
        <surface_tension>1e-3</surface_tension> 

        <!-- The bending elasticity of the cell surface [M L^2 / T^2] <! -->
        <bending_modulus>2e-18</bending_modulus> 
    </face_type>
</cell_type>
```
