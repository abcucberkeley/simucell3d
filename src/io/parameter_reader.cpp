#include "parameter_reader.hpp"


//----------------------------------------------------------------------------------------------------------------------
//Constructor   
parameter_reader::parameter_reader(const std::string& parameter_file_path){

    //The path to the parameter file must be an absolute path
    tinyxml2::XMLError read_success = xml_doc.LoadFile(parameter_file_path.c_str());

    //If the file is not found
    if(read_success != tinyxml2::XML_SUCCESS){
        throw parameter_reader_exception("The parameter file could not be found at the path: " + parameter_file_path +
        "\nThe path to the parameter file must be an absolute path.");
    }
}
//----------------------------------------------------------------------------------------------------------------------


//----------------------------------------------------------------------------------------------------------------------
//Given a particular position in the XML file (tinyxml2::XMLElement), and a given XML markup, this function
//returns the value of the XML markup as a string. If the XML markup does not exist, it returns an empty
std::optional<std::string> parameter_reader::get_string_value(const tinyxml2::XMLElement* e,  const std::string XML_markup, const bool to_lower_case /*=false*/) const noexcept{

    //Get the XML markup
    auto sting_value = e->FirstChildElement(XML_markup.c_str());

    //If the XML marrkup was not found
    if(sting_value == nullptr) return std::nullopt;

    //Convert the XML markup to a string
    std::string str = sting_value->GetText();
    return (to_lower_case) ? lower_string(str) : str;
}
//----------------------------------------------------------------------------------------------------------------------


//----------------------------------------------------------------------------------------------------------------------
//Tries to select a section of the xl file, if the section is not found throw an exception
tinyxml2::XMLElement* parameter_reader::select_section(const std::string& section_name) noexcept(false){

    tinyxml2::XMLElement* element = xml_doc.FirstChildElement(section_name.c_str());

    if(element == nullptr)  throw parameter_reader_exception("The section \""+section_name+"\" was not found in the parameter file.");

    return element;
}
//----------------------------------------------------------------------------------------------------------------------






//----------------------------------------------------------------------------------------------------------------------
//Read the simulation parameters from the parameter file
global_simulation_parameters parameter_reader::read_numerical_parameters() noexcept(false){

    //Store the collected parameters in this structure
    global_simulation_parameters sim_parameters;

    //Get the section where the input output variables of the simulation are stored
    auto io_section = select_section("numerical_parameters");

    //Get the input mesh file path
    auto input_mesh_file_path_opt = get_string_value(io_section, "input_mesh_file_path");
    if(!input_mesh_file_path_opt.has_value()) throw parameter_reader_exception("The xml markup \"input_mesh_file_path\" was not found in the parameter file.");
    sim_parameters.input_mesh_path_ = input_mesh_file_path_opt.value();

    //Get the output mesh folder path
    auto output_mesh_folder_path_opt = get_string_value(io_section, "output_mesh_folder_path");
    if(!output_mesh_folder_path_opt.has_value()) throw parameter_reader_exception("The xml markup \"output_mesh_folder_path\" was not found in the parameter file.");
    sim_parameters.output_folder_path_ = output_mesh_folder_path_opt.value();

    //Get the damping coefficient
    auto damping_coefficient_opt = get_string_value(io_section, "damping_coefficient");
    if(!damping_coefficient_opt.has_value()) throw parameter_reader_exception("The xml markup \"damping_coefficient\" was not found in the parameter file.");
    sim_parameters.damping_coefficient_ = std::stod(damping_coefficient_opt.value());
    if(sim_parameters.damping_coefficient_ < 0.0){throw parameter_reader_exception("The damping_coefficient_ must be strictly positive");}

    //Get the boolean indicating if the cells should be triangulated before the start of the simulation
    auto perform_initial_triangulation_opt = get_string_value(io_section, "perform_initial_triangulation");
    if(!perform_initial_triangulation_opt.has_value()) throw parameter_reader_exception("The xml markup \"perform_initial_triangulation\" was not found in the parameter file.");
    sim_parameters.perform_initial_triangulation_ = (std::stoi(perform_initial_triangulation_opt.value()) == 0) ? false : true;

    //Get the simulation duration
    auto simulation_duration_opt = get_string_value(io_section, "simulation_duration");
    if(!simulation_duration_opt.has_value()) throw parameter_reader_exception("The xml markup \"simulation_duration\" was not found in the parameter file.");
    sim_parameters.simulation_duration_ = std::stod(simulation_duration_opt.value());
    if(sim_parameters.simulation_duration_ <= 0.0){throw parameter_reader_exception("The simulation_duration_ must be strictly positive");}

    //Get the simulation time step
    auto time_step_opt = get_string_value(io_section, "time_step");
    if(!time_step_opt.has_value()) throw parameter_reader_exception("The xml markup \"time_step\" was not found in the parameter file.");
    sim_parameters.time_step_ = std::stod(time_step_opt.value());
    if(sim_parameters.time_step_ <= 0.0){throw parameter_reader_exception("The time_step_ must be strictly positive");}

    //Get the simulation sampling period
    auto sampling_period_opt = get_string_value(io_section, "sampling_period");
    if(!sampling_period_opt.has_value()) throw parameter_reader_exception("The xml markup \"sampling_period\" was not found in the parameter file.");
    sim_parameters.sampling_period_ = std::stod(sampling_period_opt.value());
    if(sim_parameters.sampling_period_ <= 0.0){throw parameter_reader_exception("The sampling_period_ must be strictly positive");}
    if(sim_parameters.sampling_period_ < sim_parameters.time_step_ ){throw parameter_reader_exception("The sampling_period_ must be greater than the time step");}


    //Get the min_edge_length
    auto min_edge_length_opt = get_string_value(io_section, "min_edge_length");
    if(!min_edge_length_opt.has_value()) throw parameter_reader_exception("The xml markup \"min_edge_length\" was not found in the parameter file.");
    sim_parameters.min_edge_len_ = std::stod(min_edge_length_opt.value());
    if(sim_parameters.min_edge_len_ <= 0.0){throw parameter_reader_exception("The min_edge_len_ must be strictly positive");}

    //Get the adhesion distance cutoff
    auto contact_cutoff_adhesion_opt = get_string_value(io_section, "contact_cutoff_adhesion");
    if(!contact_cutoff_adhesion_opt.has_value()) throw parameter_reader_exception("The xml markup \"contact_cutoff_adhesion\" was not found in the parameter file.");
    sim_parameters.contact_cutoff_adhesion_ = std::stod(contact_cutoff_adhesion_opt.value());
    if(sim_parameters.contact_cutoff_adhesion_ <= 0.0){throw parameter_reader_exception("The contact_cutoff_adhesion_ must be strictly positive");}

    //Get the repulsion distance cutoff
    auto contact_cutoff_repulsion_opt = get_string_value(io_section, "contact_cutoff_repulsion");
    if(!contact_cutoff_repulsion_opt.has_value()) throw parameter_reader_exception("The xml markup \"contact_cutoff_repulsion\" was not found in the parameter file.");
    sim_parameters.contact_cutoff_repulsion_ = std::stod(contact_cutoff_repulsion_opt.value());
    if(sim_parameters.contact_cutoff_repulsion_ <= 0.0){throw parameter_reader_exception("The contact_cutoff_repulsion_ must be strictly positive");}

    //Get the boolean indicating if the edge swap operation should be enabled
    auto enable_edge_swap_operation_opt = get_string_value(io_section, "enable_edge_swap_operation");
    if(!enable_edge_swap_operation_opt.has_value()) throw parameter_reader_exception("The xml markup \"enable_edge_swap_operation\" was not found in the parameter file.");
    sim_parameters.enable_edge_swap_operation_ = (std::stoi(enable_edge_swap_operation_opt.value()) == 0) ? false : true;

    return sim_parameters;
}   
//----------------------------------------------------------------------------------------------------------------------


//----------------------------------------------------------------------------------------------------------------------
//Read the biomechanical parameters of each cell type from the parameter file
std::vector<std::shared_ptr<cell_type_parameters>> parameter_reader::read_biomechanical_parameters() noexcept(false){

    //Get the cell_types root section
    auto cell_type_root_section = select_section("cell_types");

    //Make sure that at least one cell type  has been defined
    if(cell_type_root_section->FirstChildElement("cell_type") == nullptr){
        throw parameter_reader_exception("No cell type has been defined in the parameter file.");
    }

    //Store all the cell types in this vector
    std::vector<std::shared_ptr<cell_type_parameters>> cell_type_lst;

    //Loop over the cell types that are defined in the file
    unsigned cell_type_id = 0;
    for(tinyxml2::XMLElement* cell_type_section = cell_type_root_section->FirstChildElement("cell_type"); 
        cell_type_section != NULL; 
        cell_type_section = cell_type_section->NextSiblingElement("cell_type")){

        //Load the parameters of the cell type
        std::shared_ptr<cell_type_parameters> cell_parameters = read_cell_type_parameters(cell_type_id, cell_type_section);

        //Get the face type section
        tinyxml2::XMLElement* face_type_root_section = cell_type_section->FirstChildElement("face_types");
        if(face_type_root_section == nullptr){
            throw parameter_reader_exception("The xml markup \"face_types\" was not found in the cell type + "+cell_parameters->name_+" in the parameter file.");
        }

        //Make sure that at least one face type has been defined for the current cell type
        if(face_type_root_section->FirstChildElement("face_type") == nullptr){
            throw parameter_reader_exception("No face type has been defined for the cell type "+cell_parameters->name_+" in the parameter file.");
        }


        //Loop over the face types defined for the current cell type
        unsigned face_type_id = 0;
        for(tinyxml2::XMLElement* face_type_section = face_type_root_section->FirstChildElement("face_type"); face_type_section != NULL; face_type_section = face_type_section->NextSiblingElement("face_type")){

            //Load the parameters of the face type
            face_type_parameters face_parameters = read_face_type_parameters(cell_parameters->name_, face_type_id, face_type_section);

            //Add the face type parameters to the cell type parameters
            cell_parameters->face_types_.push_back(face_parameters);
            face_type_id++; 
        }

        cell_type_lst.push_back(cell_parameters);
        cell_type_id++;
    }
    return cell_type_lst;
}
//----------------------------------------------------------------------------------------------------------------------





//----------------------------------------------------------------------------------------------------------------------
//Load the parameters of the given face type
face_type_parameters parameter_reader::read_face_type_parameters(
    const std::string& cell_type_name,
    const unsigned face_type_id,     
    tinyxml2::XMLElement* face_type_section
) noexcept(false){

    assert(face_type_section != nullptr);

    //Store the collected parameters in this structure
    face_type_parameters face_parameters; 

    //First try to get the name of the face type
    auto face_type_name_opt = get_string_value(face_type_section, "face_type_name");
    if(!face_type_name_opt.has_value()) throw parameter_reader_exception("The xml markup \"face_type_name\" of the face type nb "+
    std::to_string(face_type_id)+" of the cell "+cell_type_name+" was not found in the parameter file.");
    face_parameters.name_ = face_type_name_opt.value();

    //Get the global face id (this is the face type id printed in the output mesh)
    auto global_face_id_opt = get_string_value(face_type_section, "global_face_id");
    if(!global_face_id_opt.has_value()) throw parameter_reader_exception("The xml markup \"global_face_id\" of the face "+face_parameters.name_+
    " of the cell "+cell_type_name+" was not found in the parameter file.");
    face_parameters.face_type_global_id_ = std::stoi(global_face_id_opt.value());
    if(face_parameters.face_type_global_id_ < 0) throw parameter_reader_exception("The face_type_global_id_ of the face "+face_parameters.name_+" of the cell "+cell_type_name+" is negative.");


    //Get the surface tension of the face type
    auto surface_tension_opt = get_string_value(face_type_section, "surface_tension");
    if(!surface_tension_opt.has_value()) throw parameter_reader_exception("The xml markup \"surface_tension\" of the face "+face_parameters.name_+
    " of the cell "+cell_type_name+" was not found in the parameter file.");
    face_parameters.surface_tension_ = std::stod(surface_tension_opt.value());
    if(face_parameters.surface_tension_ < 0.0) throw parameter_reader_exception("The surface_tension of the face "+face_parameters.name_+" of the cell "+cell_type_name+" is negative.");


    //Get the adherence strength of the face type
    auto adherence_strength_opt = get_string_value(face_type_section, "adherence_strength");
    if(!adherence_strength_opt.has_value()) throw parameter_reader_exception("The xml markup \"adherence_strength\" of the face "+face_parameters.name_+
    " of the cell "+cell_type_name+" was not found in the parameter file.");
    face_parameters.adherence_strength_ = std::stod(adherence_strength_opt.value());
    if(face_parameters.adherence_strength_ < 0.0) throw parameter_reader_exception("The adherence strength of the face "+face_parameters.name_+" of the cell "+cell_type_name+" is negative.");

    //Get the repulsion strength of the face type
    auto repulsion_strength_opt = get_string_value(face_type_section, "repulsion_strength");
    if(!repulsion_strength_opt.has_value()) throw parameter_reader_exception("The xml markup \"repulsion_strength\" of the face "+face_parameters.name_+
    " of the cell "+cell_type_name+" was not found in the parameter file.");
    face_parameters.repulsion_strength_ = std::stod(repulsion_strength_opt.value());
    if(face_parameters.repulsion_strength_ < 0.0) throw parameter_reader_exception("The repulsion strength of the face "+face_parameters.name_+" of the cell "+cell_type_name+" is negative.");


    //Get the bending modulus of the face type
    auto bending_modulus_opt = get_string_value(face_type_section, "bending_modulus");
    if(!bending_modulus_opt.has_value()) throw parameter_reader_exception("The xml markup \"bending_modulus\" of the face "+face_parameters.name_+
    " of the cell "+cell_type_name+" was not found in the parameter file.");
    face_parameters.bending_modulus_ = std::stod(bending_modulus_opt.value());
    if(face_parameters.bending_modulus_ < 0.0) throw parameter_reader_exception("The bending_modulus of the face "+face_parameters.name_+" of the cell "+cell_type_name+" is negative.");

    return face_parameters;
}
//----------------------------------------------------------------------------------------------------------------------




//----------------------------------------------------------------------------------------------------------------------
//Load the parameters of the given cell type but not its face type parameters
std::shared_ptr<cell_type_parameters>  parameter_reader::read_cell_type_parameters(
    const unsigned cell_type_id,     
    tinyxml2::XMLElement* cell_type_section
) noexcept(false){

    assert(cell_type_section != nullptr);

    //Store the collected parameters in this structure
    std::shared_ptr<cell_type_parameters>  cell_parameters = std::make_shared<cell_type_parameters>();

    //First try to get the name of the cell type
    auto cell_type_name_opt = get_string_value(cell_type_section, "cell_type_name");
    if(!cell_type_name_opt.has_value()) throw parameter_reader_exception("The xml markup \"cell_type_name\" of the cell type nb "+
    std::to_string(cell_type_id)+" was not found in the parameter file.");
    cell_parameters->name_ = cell_type_name_opt.value();

    //Get the global cell id (this is the cell type id printed in the output mesh)
    auto global_cell_id_opt = get_string_value(cell_type_section, "global_cell_id");
    if(!global_cell_id_opt.has_value()) throw parameter_reader_exception("The xml markup \"global_cell_id\" of the cell "+cell_parameters->name_+
    " was not found in the parameter file.");
    cell_parameters->global_type_id_ = std::stoi(global_cell_id_opt.value());

    //Get the cell mass density
    auto cell_mass_density_opt = get_string_value(cell_type_section, "cell_mass_density");
    if(!cell_mass_density_opt.has_value()) throw parameter_reader_exception("The xml markup \"cell_mass_density\" of the cell "+cell_parameters->name_+
    " was not found in the parameter file.");
    cell_parameters->mass_density_ = std::stod(cell_mass_density_opt.value());

    //Get the cell bulk modulus
    auto cell_bulk_modulus_opt = get_string_value(cell_type_section, "cell_bulk_modulus");
    if(!cell_bulk_modulus_opt.has_value()) throw parameter_reader_exception("The xml markup \"cell_bulk_modulus\" of the cell "+cell_parameters->name_+
    " was not found in the parameter file.");
    cell_parameters->bulk_modulus_ = std::stod(cell_bulk_modulus_opt.value());

    //Get the max inner pressure of the cell
    auto max_inner_pressure_opt = get_string_value(cell_type_section, "max_inner_pressure", true);
    if(!max_inner_pressure_opt.has_value()) throw parameter_reader_exception("The xml markup \"max_inner_pressure\" of the cell "+cell_parameters->name_+
    " was not found in the parameter file.");
    std::string max_inner_pressure_str = max_inner_pressure_opt.value();
    cell_parameters->max_pressure_ = (max_inner_pressure_str == "inf") ? std::numeric_limits<double>::infinity(): std::stod(max_inner_pressure_opt.value());


    //Get the cell area elasticity modulus
    auto area_elasticity_modulus_opt = get_string_value(cell_type_section, "area_elasticity_modulus");
    if(!area_elasticity_modulus_opt.has_value()) throw parameter_reader_exception("The xml markup \"area_elasticity_modulus\" of the cell "+cell_parameters->name_+
    " was not found in the parameter file.");
    cell_parameters->area_elasticity_modulus_ = std::stod(area_elasticity_modulus_opt.value());

    //Get the cell avg division volume
    auto avg_division_volume_opt = get_string_value(cell_type_section, "avg_division_volume", true);
    if(!avg_division_volume_opt.has_value()) throw parameter_reader_exception("The xml markup \"avg_division_volume\" of the cell "+cell_parameters->name_+
    " was not found in the parameter file.");
    cell_parameters->avg_division_vol_ = (avg_division_volume_opt.value() == "inf") ? std::numeric_limits<double>::infinity(): std::stod(avg_division_volume_opt.value());

    //Get the cell division volume standard deviation
    auto division_volume_std_dev_opt = get_string_value(cell_type_section, "std_division_volume");
    if(!division_volume_std_dev_opt.has_value()) throw parameter_reader_exception("The xml markup \"std_division_volume\" of the cell "+cell_parameters->name_+
    " was not found in the parameter file.");
    cell_parameters->std_division_vol_ = std::stod(division_volume_std_dev_opt.value());

    //Get the cell average growth rate
    auto avg_growth_rate_opt = get_string_value(cell_type_section, "avg_growth_rate");
    if(!avg_growth_rate_opt.has_value()) throw parameter_reader_exception("The xml markup \"avg_growth_rate\" of the cell "+cell_parameters->name_+
    " was not found in the parameter file.");
    cell_parameters->avg_growth_rate_ = std::stod(avg_growth_rate_opt.value());

    //Get the cell growth rate standard deviation
    auto growth_rate_std_dev_opt = get_string_value(cell_type_section, "std_growth_rate");
    if(!growth_rate_std_dev_opt.has_value()) throw parameter_reader_exception("The xml markup \"avg_growth_rate\" of the cell "+cell_parameters->name_+
    " was not found in the parameter file.");
    cell_parameters->std_growth_rate_ = std::stod(growth_rate_std_dev_opt.value());

    //Get the cell target isoperimetric ratio
    auto target_isoperimetric_ratio_opt = get_string_value(cell_type_section, "target_isoperimetric_ratio");
    if(!target_isoperimetric_ratio_opt.has_value()) throw parameter_reader_exception("The xml markup \"target_isoperimetric_ratio\" of the cell "+cell_parameters->name_+
    " was not found in the parameter file.");
    cell_parameters->target_isoperimetric_ratio_ = std::stod(target_isoperimetric_ratio_opt.value());
    if(cell_parameters->target_isoperimetric_ratio_ <= 0.) throw parameter_reader_exception("The xml markup \"target_isoperimetric_ratio\" of the cell "+cell_parameters->name_+
    " must be strictly positive.");

    //Get the cell angle regularization factor
    auto angle_regularization_factor_opt = get_string_value(cell_type_section, "angle_regularization_factor");
    if(!angle_regularization_factor_opt.has_value()) throw parameter_reader_exception("The xml markup \"angle_regularization_factor\" of the cell "+cell_parameters->name_+
    " was not found in the parameter file.");
    cell_parameters->angle_regularization_factor_ = std::stod(angle_regularization_factor_opt.value());

    //Get the cell minimum volume
    auto min_vol_opt = get_string_value(cell_type_section, "min_vol");
    if(!min_vol_opt.has_value()) throw parameter_reader_exception("The xml markup \"min_vol\" of the cell "+cell_parameters->name_+
    " was not found in the parameter file.");
    cell_parameters->min_vol_ = std::stod(min_vol_opt.value());


    //Get the surface coupling breaking force
    auto surface_coupling_max_curvature_opt = get_string_value(cell_type_section, "surface_coupling_max_curvature");
    if(!surface_coupling_max_curvature_opt.has_value()) throw parameter_reader_exception("The xml markup \"surface_coupling_max_curvature\" of the cell "+cell_parameters->name_+
    " was not found in the parameter file.");
    cell_parameters->surface_coupling_max_curvature_ = std::stod(surface_coupling_max_curvature_opt.value());
    if(cell_parameters->target_isoperimetric_ratio_ < 0.) throw parameter_reader_exception("The xml markup \"surface_coupling_max_curvature\" of the cell "+cell_parameters->name_+
    " must be strictly positive.");

    return cell_parameters;
}
//----------------------------------------------------------------------------------------------------------------------
