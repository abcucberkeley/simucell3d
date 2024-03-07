#include "test_parameter_reader.hpp"

//---------------------------------------------------------------------------------------------------------
int read_numerical_parameters_test(){
    
    //Open the xml parameter file
    parameter_reader xml_reader(std::string(PROJECT_SOURCE_DIR) + "/test/test_io/test_parameter_reader/test_file.xml");

    //Read the numerical parameters of the simulation
    global_simulation_parameters sim_parameters = xml_reader.read_numerical_parameters();

    //Make sure that all the values are correct
    bool t1 = sim_parameters.output_folder_path_ == "./simulation_results";
    bool t2 = sim_parameters.input_mesh_path_ == "../data/Geometries/Cube.vtk";
    bool t3 = sim_parameters.damping_coefficient_ == 1000000.0;
    bool t4 = sim_parameters.simulation_duration_ == 0.001;
    bool t5 = sim_parameters.sampling_period_ == 1e-06;
    bool t6 = sim_parameters.time_step_ == 1e-07;
    bool t7 = sim_parameters.min_edge_len_ == 2.5e-07;
    bool t8 = sim_parameters.contact_cutoff_adhesion_ == 2.5e-07;
    bool t9 = sim_parameters.contact_cutoff_repulsion_ == 2.5e-07;
    bool t10 = sim_parameters.perform_initial_triangulation_ == true;

    std::cout << "t1 = " << t1 << std::endl;
    std::cout << "t2 = " << t2 << std::endl;
    std::cout << "t3 = " << t3 << std::endl;
    std::cout << "t4 = " << t4 << std::endl;
    std::cout << "t5 = " << t5 << std::endl;
    std::cout << "t6 = " << t6 << std::endl;
    std::cout << "t7 = " << t7 << std::endl;
    std::cout << "t8 = " << t8 << std::endl;
    std::cout << "t9 = " << t9 << std::endl;
    std::cout << "t10 = " << t10 << std::endl;

    return !(t1 && t2 && t3 && t4 && t5 && t6 && t7 && t8 && t9 && t10);
}
//---------------------------------------------------------------------------------------------------------




//---------------------------------------------------------------------------------------------------------
int test_parameter_reader::read_cell_type_parameters_test() const{

    //Open the xml parameter file
    parameter_reader xml_reader(std::string(PROJECT_SOURCE_DIR) + "/test/test_io/test_parameter_reader/test_file.xml");


    //Get the cell_types root section
    tinyxml2::XMLElement* cell_type_root_section;
    try{
        cell_type_root_section = xml_reader.select_section("cell_types");
    }
    catch(const parameter_reader_exception& e){
        std::cout << "The section \"cell_types\" could not be found" << std::endl;
        return 1; 
    }
  

    //Try to get the first cell_type
    auto cell_type_section = cell_type_root_section->FirstChildElement("cell_type");

    if(cell_type_section == nullptr){
        std::cout << "No section \"cell_type\" could be found" << std::endl;
        return 1;
    }

    //Load the parameters of the given cell type but not its face type parameters
    auto cell_type_param  = xml_reader.read_cell_type_parameters(0,  cell_type_section);

    //Make sure the parameters are correct
    bool t1 = cell_type_param->name_ == "epithelial";
    bool t2 = cell_type_param->global_type_id_ == 0;
    bool t3 = cell_type_param->mass_density_ == 1.0e3;
    bool t4 = cell_type_param->bulk_modulus_ == 2500;
    bool t5 = cell_type_param->max_pressure_ == std::numeric_limits<double>::infinity();
    bool t6 = cell_type_param->area_elasticity_modulus_ == 0;
    bool t7 = cell_type_param->avg_division_vol_ == 2e-16;
    bool t8 = cell_type_param->std_division_vol_ == 0;
    bool t9 = cell_type_param->avg_growth_rate_ == 1e-12;
    bool t10 = cell_type_param->std_growth_rate_ == 0;
    bool t11 = cell_type_param->min_vol_ == 1e-18;
    bool t12 = cell_type_param->target_isoperimetric_ratio_ == 130;
    bool t15 = cell_type_param->angle_regularization_factor_ == 5e-17;


    std::cout << "t1 = " << t1 << std::endl;
    std::cout << "t2 = " << t2 << std::endl;
    std::cout << "t3 = " << t3 << std::endl;
    std::cout << "t4 = " << t4 << std::endl;
    std::cout << "t5 = " << t5 << std::endl;
    std::cout << "t6 = " << t6 << std::endl;   
    std::cout << "t7 = " << t7 << std::endl;
    std::cout << "t8 = " << t8 << std::endl;
    std::cout << "t9 = " << t9 << std::endl;
    std::cout << "t10 = " << t10 << std::endl;
    std::cout << "t11 = " << t11 << std::endl;
    std::cout << "t12 = " << t12 << std::endl;


    return !(t1 && t2 && t3 && t4 && t5 && t6 && t7 && t8 && t9 && t10 && t11 && t12 && t15);

}
//---------------------------------------------------------------------------------------------------------






//---------------------------------------------------------------------------------------------------------
int test_parameter_reader::read_face_type_parameters_test() const{

    //Open the xml parameter file
    parameter_reader xml_reader(std::string(PROJECT_SOURCE_DIR) + "/test/test_io/test_parameter_reader/test_file.xml");


    //Get the cell_types root section
    tinyxml2::XMLElement* cell_type_root_section;
    try{
        cell_type_root_section = xml_reader.select_section("cell_types");
    }
    catch(const parameter_reader_exception& e){
        std::cout << "The section \"cell_types\" could not be found" << std::endl;
        return 1; 
    }
  

    //Try to get the first cell_type
    auto cell_type_section = cell_type_root_section->FirstChildElement("cell_type");

    if(cell_type_section == nullptr){
        std::cout << "No section \"cell_type\" could be found" << std::endl;
        return 1;
    }

    //Get the face_types section
    tinyxml2::XMLElement* face_type_root_section = cell_type_section->FirstChildElement("face_types");
    if(face_type_root_section== nullptr){
        std::cout << "The face_types section has not been defined for the first cell type in the parameter file." << std::endl;
        return 1;
    }


    //Get the first face type
    tinyxml2::XMLElement* face_type_section = face_type_root_section->FirstChildElement("face_type");
    if(face_type_section == nullptr){
        std::cout << "No face type has been defined for the first cell type in the parameter file." << std::endl;
        return 1; 
    }

    //Load the parameters of the given face type
    auto face_type_param  = xml_reader.read_face_type_parameters("epithelial", 0,  face_type_section);

    //Make sure the parameters are correct
    bool t1 = face_type_param.name_ == "apical";
    bool t2 = face_type_param.face_type_global_id_ == 0;
    bool t3 = face_type_param.surface_tension_ == 1e-3;
    bool t4 = face_type_param.adherence_strength_ == 2.2e9;
    bool t5 = face_type_param.repulsion_strength_ == 1e9;
    bool t6 = face_type_param.bending_modulus_ == 2e-18;

    std::cout << "t1 = " << t1 << std::endl;
    std::cout << "t2 = " << t2 << std::endl;
    std::cout << "t3 = " << t3 << std::endl;
    std::cout << "t4 = " << t4 << std::endl;
    std::cout << "t5 = " << t5 << std::endl;
    std::cout << "t6 = " << t6 << std::endl;

    return !(t1 && t2 && t3 && t4 && t5 && t6);

}
//---------------------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------------------

//Read the whole file and make random tests
int test_parameter_reader::read_biomechanical_parameters_test() const{

    //Open the xml parameter file
    parameter_reader xml_reader(std::string(PROJECT_SOURCE_DIR) + "/test/test_io/test_parameter_reader/test_file.xml");

    //Load all the parameters of the cells
    std::vector<std::shared_ptr<cell_type_parameters>> cell_type_param_lst = xml_reader.read_biomechanical_parameters();

    if(cell_type_param_lst.size() != 5) return 1;


    //Make sure some of the parameters are correct, I don't have the strength to check all of them
    bool t1 = cell_type_param_lst[1]->name_ == "ecm";
    bool t2 = cell_type_param_lst[1]->face_types_.size() == 1;
    bool t3 = t2 ? cell_type_param_lst[1]->face_types_[0].name_ == "ecm_face" : false;
    bool t4 = t2 ? cell_type_param_lst[1]->face_types_[0].adherence_strength_ == 2.2e9 : false;
    
    bool t5 = cell_type_param_lst[4]->target_isoperimetric_ratio_ == 216.;
    bool t6 = cell_type_param_lst[4]->name_ == "static_cell";
    bool t7 = cell_type_param_lst[4]->face_types_.size() == 1;
    bool t8 = t7 ? cell_type_param_lst[4]->face_types_[0].name_ == "static_face" : false;
    bool t9 = t7 ? cell_type_param_lst[4]->face_types_[0].adherence_strength_ == 0 : false;
    bool t10 = t7 ? cell_type_param_lst[4]->face_types_[0].repulsion_strength_ == 1e9 : false;

    
    std::cout << "t1 = " << t1 << std::endl;
    std::cout << "t2 = " << t2 << std::endl;
    std::cout << "t3 = " << t3 << std::endl;
    std::cout << "t4 = " << t4 << std::endl;
    std::cout << "t5 = " << t5 << std::endl;
    std::cout << "t6 = " << t6 << std::endl;
    std::cout << "t7 = " << t7 << std::endl;
    std::cout << "t8 = " << t8 << std::endl;
    std::cout << "t9 = " << t9 << std::endl;
    std::cout << "t10 = " << t10 << std::endl;

    return !(t1 && t2 && t3 && t4 && t5 && t6 && t7 && t8 && t9 && t10);
}
//---------------------------------------------------------------------------------------------------------






//---------------------------------------------------------------------------------------------------------
// The main function
int main (int argc, char** argv){

    //Check that the command line input is correctly formatted
    assert(argc == 2); 

    //Get the name of the test to run
    std::string test_name = argv[1];

    test_parameter_reader tester;
    
    //Run the selected test
    if (test_name == "read_numerical_parameters_test")           return read_numerical_parameters_test();
    if (test_name == "read_cell_type_parameters_test")           return tester.read_cell_type_parameters_test();
    if (test_name == "read_face_type_parameters_test")       return tester.read_face_type_parameters_test();
    if (test_name == "read_biomechanical_parameters_test")       return tester.read_biomechanical_parameters_test();

    std::cout << "TEST NAME :" << test_name << " DOES NOT EXIST" << std::endl;
    return 1;
}
//---------------------------------------------------------------------------------------------------------



