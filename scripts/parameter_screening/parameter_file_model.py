"""
    This file contains a model of the xml file used to parametrize the simulations.  This model
    is used to generate the various xml files used to run the simulations of the screening.
"""

def get_xml_file_model(par):
    """
        This function generates the xml file model used to parametrize the simulations.

        Parameters:
        -----------
        
        par (dict):
            The parameters of the simulation.

        Returns:
        --------
        
        xml_file_model (str):
            The xml file model.
    """

    xml_file_model = """
<?xml version="1.0"?>

<numerical_parameters>
    <!-- The mesh file containing the initial configuration <! -->
    <input_mesh_file_path>{input_mesh_file_path}</input_mesh_file_path> 

    <!-- The path where the mesh files will be written <! -->
    <output_mesh_folder_path>{output_mesh_folder_path}</output_mesh_folder_path>

    <!-- If equals 0 the cells will not be triangulated at the beginning of the simulation \<! -->
    <perform_initial_triangulation>{perform_initial_triangulation}</perform_initial_triangulation> 
 
    <!-- The damping coefficients [M / T] <! -->
    <damping_coefficient>{damping_coefficient}</damping_coefficient> 

    <!-- The total simulation time [T] <! -->
    <simulation_duration>{simulation_duration}</simulation_duration> 

    <!-- The time between each mesh file saving [T]<! -->
    <sampling_period>{sampling_period}</sampling_period> 

    <!-- The simulation time step [T]<! -->
    <time_step>{time_step}</time_step> 

    <!-- The minimum edge length allowed [L] <! -->
    <min_edge_length>{min_edge_length}</min_edge_length> 

    <!-- The interaction cutoff distance above which cell-cell interactions are not considered anymore[L]<! -->
    <contact_cutoff_adhesion>{contact_cutoff_adhesion}</contact_cutoff_adhesion> 
    <contact_cutoff_repulsion>{contact_cutoff_repulsion}</contact_cutoff_repulsion> 

    <!-- The maximum force in Newtons that can be applied to a node [M L / T^2]<! -->
    <!-- It can be set to "inf" if you do not wish to cap the forces of the nodes <! -->
    <enable_edge_swap_operation>{enable_edge_swap_operation}</enable_edge_swap_operation> 
</numerical_parameters>



<cell_types>
    <cell_type>
        <!-- The name of the cell type <! -->
        <cell_type_name>epithelial</cell_type_name> 

        <!-- The ID of the cell type <! -->
        <global_cell_id>0</global_cell_id>  

        <!-- The cell density [M /(L^3)]<! -->
        <cell_mass_density>{epi_mass_density}</cell_mass_density> 

        <!-- The bulk modulus of the cells [M /(L T^2)] <! -->
        <cell_bulk_modulus>{epi_bulk_modulus}</cell_bulk_modulus> 

        <!-- The maximal inner pressure of a cell. Set to INF to have no pressure limit [M / (T^2 L)] <! -->
        <max_inner_pressure>{epi_max_inner_pressure}</max_inner_pressure> 

        <!-- Growth rate = cellDivisionTargetVolume / cellCycleDuration [L^3 / T] Default: 3e-12 <! -->
        <avg_growth_rate>{epi_avg_growth_rate}</avg_growth_rate> 
        <std_growth_rate>{epi_std_growth_rate}</std_growth_rate> 

        <!-- Isoperimetric ratio = cellArea^3 / cellVolume^2 <! -->
        <target_isoperimetric_ratio>{epi_target_isoperimetric_ratio}</target_isoperimetric_ratio> 

        <!-- The elastic modulus of the cell membranes [M / (T^2 L^2)] <! -->
        <area_elasticity_modulus>{epi_area_elasticity_modulus}</area_elasticity_modulus> 

        <!-- Used to make sure that the angles of the triangular faces are all approximately 60 deg [M L^2 / T^2] <! -->
        <angle_regularization_factor>{epi_angle_regularization_factor}</angle_regularization_factor> 


        <!-- The volume that the cell has to reach before division [L^3] (Drawn from a normal distribution)<! -->
        <avg_division_volume>{epi_avg_division_volume}</avg_division_volume> 
        <std_division_volume>{epi_std_division_volume}</std_division_volume> 

        <!-- The forces at which the links between the epithelial cells are broken  [M L / T^2] <! -->
        <surface_coupling_max_curvature>{epi_surface_coupling_max_curvature}</surface_coupling_max_curvature>


        <!-- If the cell reach this volume we delete it[L^3]<! -->
        <min_vol>{epi_min_vol}</min_vol> 

        <!-- The cell might have different face types with different mechanical properties <! -->
        <face_types>     
            <face_type>           
                <global_face_id>0</global_face_id>

                <face_type_name>apical</face_type_name>

                <!-- Modulate the adherence between cells [(M / (T^2 L^2)] <! -->
                <adherence_strength>{epi_api_adherence_strength}</adherence_strength>

                <!-- Modulate the repulsion between cells [(M / (T^2 L^2)] <! -->
                <repulsion_strength>{epi_api_repulsion_strength}</repulsion_strength> 

                <!-- The surface tension of the cell membranes [M / T^2] <! -->
                <surface_tension>{epi_api_surface_tension}</surface_tension> 

                <!-- The bending elasticity of the cell surface [M L^2 / T^2] <! -->
                <bending_modulus>{epi_api_bending_modulus}</bending_modulus> 
            </face_type>

            <face_type>           
                <global_face_id>1</global_face_id>

                <face_type_name>lateral</face_type_name>

                <!-- Modulate the adherence between cells [(M / (T^2 L^2)] <! -->
                <adherence_strength>{epi_lat_adherence_strength}</adherence_strength>

                <!-- Modulate the repulsion between cells [(M / (T^2 L^2)] <! -->
                <repulsion_strength>{epi_lat_repulsion_strength}</repulsion_strength> 

                <!-- The surface tension of the cell membranes [M / T^2] <! -->
                <surface_tension>{epi_lat_surface_tension}</surface_tension> 

                <!-- The bending elasticity of the cell surface [M L^2 / T^2] <! -->
                <bending_modulus>{epi_lat_bending_modulus}</bending_modulus> 
            </face_type>

           <face_type>           
                <global_face_id>2</global_face_id>

                <face_type_name>basal</face_type_name>

                <!-- Modulate the adherence between cells [(M / (T^2 L^2)] <! -->
                <adherence_strength>{epi_bas_adherence_strength}</adherence_strength>

                <!-- Modulate the repulsion between cells [(M / (T^2 L^2)] <! -->
                <repulsion_strength>{epi_bas_repulsion_strength}</repulsion_strength> 

                <!-- The surface tension of the cell membranes [M / T^2] <! -->
                <surface_tension>{epi_bas_surface_tension}</surface_tension> 

                <!-- The bending elasticity of the cell surface [M L^2 / T^2] <! -->
                <bending_modulus>{epi_bas_bending_modulus}</bending_modulus> 
            </face_type>
        </face_types>  
    </cell_type>

    <cell_type>
        <!-- The ID of the cell type <! -->
        <global_cell_id>1</global_cell_id>  

        <!-- The name of the cell type <! -->
        <cell_type_name>ecm</cell_type_name> 

                <!-- The cell density [M /(L^3)]<! -->
        <cell_mass_density>{ecm_mass_density}</cell_mass_density> 

        <!-- The bulk modulus of the cells [M /(L T^2)] <! -->
        <cell_bulk_modulus>{ecm_bulk_modulus}</cell_bulk_modulus> 

        <!-- The maximal inner pressure of a cell. Set to INF to have no pressure limit [M / (T^2 L)] <! -->
        <max_inner_pressure>{ecm_max_inner_pressure}</max_inner_pressure> 

        <!-- Growth rate = cellDivisionTargetVolume / cellCycleDuration [L^3 / T] Default: 3e-12 <! -->
        <avg_growth_rate>{ecm_avg_growth_rate}</avg_growth_rate> 
        <std_growth_rate>{ecm_std_growth_rate}</std_growth_rate> 

        <!-- Isoperimetric ratio = cellArea^3 / cellVolume^2 <! -->
        <target_isoperimetric_ratio>{ecm_target_isoperimetric_ratio}</target_isoperimetric_ratio> 

        <!-- The elastic modulus of the cell membranes [M / (T^2 L^2)] <! -->
        <area_elasticity_modulus>{ecm_area_elasticity_modulus}</area_elasticity_modulus> 

        <!-- Used to make sure that the angles of the triangular faces are all approximately 60 deg [M L^2 / T^2] <! -->
        <angle_regularization_factor>{ecm_angle_regularization_factor}</angle_regularization_factor> 

        <!-- The volume that the cell has to reach before division [L^3] (Drawn from a normal distribution)<! -->
        <avg_division_volume>{ecm_avg_division_volume}</avg_division_volume> 
        <std_division_volume>{ecm_std_division_volume}</std_division_volume> 

        <!-- The forces at which the links between the epithelial cells are broken  [M L / T^2] <! -->
        <surface_coupling_max_curvature>{ecm_surface_coupling_max_curvature}</surface_coupling_max_curvature>

        <!-- If the cell reach this volume we delete it[L^3]<! -->
        <min_vol>{ecm_min_vol}</min_vol> 

        <!-- The cell might have different face types with different mechanical properties <! -->
        <face_types>     
            <face_type>           
                <global_face_id>3</global_face_id>

                <face_type_name>ecm_face</face_type_name>

                <!-- Modulate the adherence between cells [(M / (T^2 L^2)] <! -->
                <adherence_strength>{ecm_adherence_strength}</adherence_strength>

                <!-- Modulate the repulsion between cells [(M / (T^2 L^2)] <! -->
                <repulsion_strength>{ecm_repulsion_strength}</repulsion_strength> 

                <!-- The surface tension of the cell membranes [M / T^2] <! -->
                <surface_tension>{ecm_surface_tension}</surface_tension> 

                <!-- The bending elasticity of the cell surface [M L^2 / T^2] <! -->
                <bending_modulus>{ecm_bending_modulus}</bending_modulus> 
            </face_type>
        </face_types>  
    </cell_type>


    <cell_type>
        <!-- The ID of the cell type <! -->
        <global_cell_id>2</global_cell_id>  

        <!-- The name of the cell type <! -->
        <cell_type_name>lumen</cell_type_name> 

         <!-- The cell density [M /(L^3)]<! -->
        <cell_mass_density>{lum_mass_density}</cell_mass_density> 

        <!-- The bulk modulus of the cells [M /(L T^2)] <! -->
        <cell_bulk_modulus>{lum_bulk_modulus}</cell_bulk_modulus> 

        <!-- The maximal inner pressure of a cell. Set to INF to have no pressure limit [M / (T^2 L)] <! -->
        <max_inner_pressure>{lum_max_inner_pressure}</max_inner_pressure> 

        <!-- Growth rate = cellDivisionTargetVolume / cellCycleDuration [L^3 / T] Default: 3e-12 <! -->
        <avg_growth_rate>{lum_avg_growth_rate}</avg_growth_rate> 
        <std_growth_rate>{lum_std_growth_rate}</std_growth_rate> 

        <!-- Isoperimetric ratio = cellArea^3 / cellVolume^2 <! -->
        <target_isoperimetric_ratio>{lum_target_isoperimetric_ratio}</target_isoperimetric_ratio> 

        <!-- The elastic modulus of the cell membranes [M / (T^2 L^2)] <! -->
        <area_elasticity_modulus>{lum_area_elasticity_modulus}</area_elasticity_modulus> 

        <!-- Used to make sure that the angles of the triangular faces are all approximately 60 deg [M L^2 / T^2] <! -->
        <angle_regularization_factor>{lum_angle_regularization_factor}</angle_regularization_factor> 

        <!-- The volume that the cell has to reach before division [L^3] (Drawn from a normal distribution)<! -->
        <avg_division_volume>{lum_avg_division_volume}</avg_division_volume> 
        <std_division_volume>{lum_std_division_volume}</std_division_volume>

        <!-- The forces at which the links between the epithelial cells are broken  [M L / T^2] <! -->
        <surface_coupling_max_curvature>{lum_surface_coupling_max_curvature}</surface_coupling_max_curvature>
 
        <!-- If the cell reach this volume we delete it[L^3]<! -->
        <min_vol>{lum_min_vol}</min_vol> 

        <!-- The cell might have different face types with different mechanical properties <! -->
        <face_types>     
            <face_type>           
                <global_face_id>4</global_face_id>

                <face_type_name>lumen_face</face_type_name>

                <!-- Modulate the adherence between cells [(M / (T^2 L^2)] <! -->
                <adherence_strength>{lum_adherence_strength}</adherence_strength>

                <!-- Modulate the repulsion between cells [(M / (T^2 L^2)] <! -->
                <repulsion_strength>{lum_repulsion_strength}</repulsion_strength> 

                <!-- The surface tension of the cell membranes [M / T^2] <! -->
                <surface_tension>{lum_surface_tension}</surface_tension> 

                <!-- The bending elasticity of the cell surface [M L^2 / T^2] <! -->
                <bending_modulus>{lum_bending_modulus}</bending_modulus> 
            </face_type>
        </face_types>  
    </cell_type>

    <cell_type>
        <!-- The ID of the cell type <! -->
        <global_cell_id>3</global_cell_id>  

        <!-- The name of the cell type <! -->
        <cell_type_name>nucleus</cell_type_name> 

        <!-- The cell density [M /(L^3)]<! -->
        <cell_mass_density>{nuc_mass_density}</cell_mass_density> 

        <!-- The bulk modulus of the cells [M /(L T^2)] <! -->
        <cell_bulk_modulus>{nuc_bulk_modulus}</cell_bulk_modulus> 

        <!-- The maximal inner pressure of a cell. Set to INF to have no pressure limit [M / (T^2 L)] <! -->
        <max_inner_pressure>{nuc_max_inner_pressure}</max_inner_pressure> 

        <!-- Growth rate = cellDivisionTargetVolume / cellCycleDuration [L^3 / T] Default: 3e-12 <! -->
        <avg_growth_rate>{nuc_avg_growth_rate}</avg_growth_rate> 
        <std_growth_rate>{nuc_std_growth_rate}</std_growth_rate> 

        <!-- Isoperimetric ratio = cellArea^3 / cellVolume^2 <! -->
        <target_isoperimetric_ratio>{nuc_target_isoperimetric_ratio}</target_isoperimetric_ratio> 

        <!-- The elastic modulus of the cell membranes [M / (T^2 L^2)] <! -->
        <area_elasticity_modulus>{nuc_area_elasticity_modulus}</area_elasticity_modulus> 

        <!-- Used to make sure that the angles of the triangular faces are all approximately 60 deg [M L^2 / T^2] <! -->
        <angle_regularization_factor>{nuc_angle_regularization_factor}</angle_regularization_factor> 

        <!-- The volume that the cell has to reach before division [L^3] (Drawn from a normal distribution)<! -->
        <avg_division_volume>{nuc_avg_division_volume}</avg_division_volume> 
        <std_division_volume>{nuc_std_division_volume}</std_division_volume> 

        <!-- The forces at which the links between the epithelial cells are broken  [M L / T^2] <! -->
        <surface_coupling_max_curvature>{nuc_surface_coupling_max_curvature}</surface_coupling_max_curvature>

        <!-- If the cell reach this volume we delete it[L^3]<! -->
        <min_vol>{nuc_min_vol}</min_vol> 

        <!-- The force applied to the nucleus to simulmate the INM [M L / T^2]<! -->
        <INMForce>0</INMForce> 

        <!-- The cell might have different face types with different mechanical properties <! -->
        <face_types>     
            <face_type>           
                <global_face_id>5</global_face_id>

                <face_type_name>nucleus_face</face_type_name>

                <!-- Modulate the adherence between cells [(M / (T^2 L^2)] <! -->
                <adherence_strength>{nuc_adherence_strength}</adherence_strength>

                <!-- Modulate the repulsion between cells [(M / (T^2 L^2)] <! -->
                <repulsion_strength>{nuc_repulsion_strength}</repulsion_strength> 

                <!-- The surface tension of the cell membranes [M / T^2] <! -->
                <surface_tension>{nuc_surface_tension}</surface_tension> 

                <!-- The bending elasticity of the cell surface [M L^2 / T^2] <! -->
                <bending_modulus>{nuc_bending_modulus}</bending_modulus> 
            </face_type>
        </face_types>  
    </cell_type>

    <cell_type>
        <!-- The ID of the cell type <! -->
        <global_cell_id>4</global_cell_id>  

        <!-- The name of the cell type <! -->
        <cell_type_name>static_cell</cell_type_name> 

               <!-- The cell density [M /(L^3)]<! -->
        <cell_mass_density>{sta_mass_density}</cell_mass_density> 

        <!-- The bulk modulus of the cells [M /(L T^2)] <! -->
        <cell_bulk_modulus>{sta_bulk_modulus}</cell_bulk_modulus> 

        <!-- The maximal inner pressure of a cell. Set to INF to have no pressure limit [M / (T^2 L)] <! -->
        <max_inner_pressure>{sta_max_inner_pressure}</max_inner_pressure> 

        <!-- Growth rate = cellDivisionTargetVolume / cellCycleDuration [L^3 / T] Default: 3e-12 <! -->
        <avg_growth_rate>{sta_avg_growth_rate}</avg_growth_rate> 
        <std_growth_rate>{sta_std_growth_rate}</std_growth_rate> 

        <!-- Isoperimetric ratio = cellArea^3 / cellVolume^2 <! -->
        <target_isoperimetric_ratio>{sta_target_isoperimetric_ratio}</target_isoperimetric_ratio> 

        <!-- The elastic modulus of the cell membranes [M / (T^2 L^2)] <! -->
        <area_elasticity_modulus>{sta_area_elasticity_modulus}</area_elasticity_modulus> 

        <!-- Used to make sure that the angles of the triangular faces are all approximately 60 deg [M L^2 / T^2] <! -->
        <angle_regularization_factor>{sta_angle_regularization_factor}</angle_regularization_factor> 

        <!-- The volume that the cell has to reach before division [L^3] (Drawn from a normal distribution)<! -->
        <avg_division_volume>{sta_avg_division_volume}</avg_division_volume> 
        <std_division_volume>{sta_std_division_volume}</std_division_volume> 

        <!-- The forces at which the links between the epithelial cells are broken  [M L / T^2] <! -->
        <surface_coupling_max_curvature>{sta_surface_coupling_max_curvature}</surface_coupling_max_curvature>

        <!-- If the cell reach this volume we delete it[L^3]<! -->
        <min_vol>{sta_min_vol}</min_vol> 

        <!-- The cell might have different face types with different mechanical properties <! -->
        <face_types>     
            <face_type>           
                <global_face_id>6</global_face_id>

                <face_type_name>static_face</face_type_name>

                <!-- Modulate the adherence between cells [(M / (T^2 L^2)] <! -->
                <adherence_strength>{sta_adherence_strength}</adherence_strength>

                <!-- Modulate the repulsion between cells [(M / (T^2 L^2)] <! -->
                <repulsion_strength>{sta_repulsion_strength}</repulsion_strength> 

                <!-- The surface tension of the cell membranes [M / T^2] <! -->
                <surface_tension>{sta_surface_tension}</surface_tension> 

                <!-- The bending elasticity of the cell surface [M L^2 / T^2] <! -->
                <bending_modulus>{sta_bending_modulus}</bending_modulus> 
            </face_type>
        </face_types>  
    </cell_type>
</cell_types>
    """.format(
        input_mesh_file_path =              par["input_mesh_file_path"],
        output_mesh_folder_path =           par["output_mesh_folder_path"],
        perform_initial_triangulation =                 par["perform_initial_triangulation"],
        damping_coefficient =               par["damping_coefficient"],
        simulation_duration =               par["simulation_duration"],
        sampling_period =                   par["sampling_period"],
        time_step =                         par["time_step"],
        min_edge_length =                   par["min_edge_length"],
        contact_cutoff_adhesion =           par["contact_cutoff_adhesion"],
        contact_cutoff_repulsion =          par["contact_cutoff_repulsion"],
        enable_edge_swap_operation =        par["enable_edge_swap_operation"],   
        #--------------------------------------------------------------------------------------------------

        #--------------------------------------------------------------------------------------------------
        epi_mass_density =                      par["epi_mass_density"],
        epi_bulk_modulus =                      par["epi_bulk_modulus"],
        epi_max_inner_pressure =                par["epi_max_inner_pressure"],
        epi_avg_growth_rate =                   par["epi_avg_growth_rate"],
        epi_std_growth_rate =                   par["epi_std_growth_rate"],
        epi_target_isoperimetric_ratio =        par["epi_target_isoperimetric_ratio"],
        epi_angle_regularization_factor =       par["epi_angle_regularization_factor"],
        epi_area_elasticity_modulus =           par["epi_area_elasticity_modulus"],
        epi_avg_division_volume =               par["epi_avg_division_volume"],
        epi_std_division_volume =               par["epi_std_division_volume"],
        epi_surface_coupling_max_curvature =    par["epi_surface_coupling_max_curvature"],
        
        epi_min_vol =                           par["epi_min_vol"],
        epi_api_adherence_strength =            par["epi_api_adherence_strength"],
        epi_api_repulsion_strength =            par["epi_api_repulsion_strength"],
        epi_api_surface_tension =               par["epi_api_surface_tension"],
        epi_api_bending_modulus =               par["epi_api_bending_modulus"],

        epi_lat_adherence_strength =            par["epi_lat_adherence_strength"],
        epi_lat_repulsion_strength =            par["epi_lat_repulsion_strength"],
        epi_lat_surface_tension =               par["epi_lat_surface_tension"],
        epi_lat_bending_modulus =               par["epi_lat_bending_modulus"],

        epi_bas_adherence_strength =            par["epi_bas_adherence_strength"],
        epi_bas_repulsion_strength =            par["epi_bas_repulsion_strength"],
        epi_bas_surface_tension =               par["epi_bas_surface_tension"],
        epi_bas_bending_modulus =               par["epi_bas_bending_modulus"],
        #--------------------------------------------------------------------------------------------------

        #--------------------------------------------------------------------------------------------------
        ecm_mass_density =                      par["ecm_mass_density"],
        ecm_bulk_modulus =                      par["ecm_bulk_modulus"],
        ecm_max_inner_pressure =                par["ecm_max_inner_pressure"],
        ecm_avg_growth_rate =                   par["ecm_avg_growth_rate"],
        ecm_std_growth_rate =                   par["ecm_std_growth_rate"],
        ecm_target_isoperimetric_ratio =        par["ecm_target_isoperimetric_ratio"],
        ecm_angle_regularization_factor =       par["ecm_angle_regularization_factor"],
        ecm_area_elasticity_modulus =           par["ecm_area_elasticity_modulus"],
        ecm_avg_division_volume =               par["ecm_avg_division_volume"],
        ecm_std_division_volume =               par["ecm_std_division_volume"],
        ecm_surface_coupling_max_curvature =    par["ecm_surface_coupling_max_curvature"],
        ecm_min_vol =                           par["ecm_min_vol"],

        ecm_adherence_strength =                par["ecm_adherence_strength"],
        ecm_repulsion_strength =                par["ecm_repulsion_strength"],
        ecm_surface_tension =                   par["ecm_surface_tension"],
        ecm_bending_modulus =                   par["ecm_bending_modulus"],
        #--------------------------------------------------------------------------------------------------

        #--------------------------------------------------------------------------------------------------
        lum_mass_density =                      par["lum_mass_density"],
        lum_bulk_modulus =                      par["lum_bulk_modulus"],
        lum_max_inner_pressure =                par["lum_max_inner_pressure"],
        lum_avg_growth_rate =                   par["lum_avg_growth_rate"],
        lum_std_growth_rate =                   par["lum_std_growth_rate"],
        lum_target_isoperimetric_ratio =        par["lum_target_isoperimetric_ratio"],
        lum_angle_regularization_factor =       par["lum_angle_regularization_factor"],
        lum_area_elasticity_modulus =           par["lum_area_elasticity_modulus"],
        lum_avg_division_volume =               par["lum_avg_division_volume"],
        lum_std_division_volume =               par["lum_std_division_volume"],
        lum_surface_coupling_max_curvature =    par["lum_surface_coupling_max_curvature"],
        lum_min_vol =                           par["lum_min_vol"],

        lum_adherence_strength =                par["lum_adherence_strength"],
        lum_repulsion_strength =                par["lum_repulsion_strength"],
        lum_surface_tension =                   par["lum_surface_tension"],
        lum_bending_modulus =                   par["lum_bending_modulus"],
        #--------------------------------------------------------------------------------------------------

        #--------------------------------------------------------------------------------------------------
        nuc_mass_density =                      par["nuc_mass_density"],
        nuc_bulk_modulus =                      par["nuc_bulk_modulus"],
        nuc_max_inner_pressure =                par["nuc_max_inner_pressure"],
        nuc_avg_growth_rate =                   par["nuc_avg_growth_rate"],
        nuc_std_growth_rate =                   par["nuc_std_growth_rate"],
        nuc_target_isoperimetric_ratio =        par["nuc_target_isoperimetric_ratio"],
        nuc_angle_regularization_factor =       par["nuc_angle_regularization_factor"],
        nuc_area_elasticity_modulus =           par["nuc_area_elasticity_modulus"],
        nuc_avg_division_volume =               par["nuc_avg_division_volume"],
        nuc_std_division_volume =               par["nuc_std_division_volume"],
        nuc_surface_coupling_max_curvature =    par["nuc_surface_coupling_max_curvature"],
        nuc_min_vol =                           par["nuc_min_vol"],

        nuc_adherence_strength =                par["nuc_adherence_strength"],
        nuc_repulsion_strength =                par["nuc_repulsion_strength"],
        nuc_surface_tension =                   par["nuc_surface_tension"],
        nuc_bending_modulus =                   par["nuc_bending_modulus"],
        #--------------------------------------------------------------------------------------------------

        #--------------------------------------------------------------------------------------------------
        sta_mass_density =                      par["sta_mass_density"],
        sta_bulk_modulus =                      par["sta_bulk_modulus"],
        sta_max_inner_pressure =                par["sta_max_inner_pressure"],
        sta_avg_growth_rate =                   par["sta_avg_growth_rate"],
        sta_std_growth_rate =                   par["sta_std_growth_rate"],
        sta_target_isoperimetric_ratio =        par["sta_target_isoperimetric_ratio"],
        sta_angle_regularization_factor =       par["sta_angle_regularization_factor"],
        sta_area_elasticity_modulus =           par["sta_area_elasticity_modulus"],
        sta_avg_division_volume =               par["sta_avg_division_volume"],
        sta_std_division_volume =               par["sta_std_division_volume"],
        sta_surface_coupling_max_curvature =    par["sta_surface_coupling_max_curvature"],
        sta_min_vol =                           par["sta_min_vol"],

        sta_adherence_strength =                par["sta_adherence_strength"],
        sta_repulsion_strength =                par["sta_repulsion_strength"],
        sta_surface_tension =                   par["sta_surface_tension"],
        sta_bending_modulus =                   par["sta_bending_modulus"],
        #--------------------------------------------------------------------------------------------------

        #--------------------------------------------------------------------------------------------------
    )
    return xml_file_model