#ifndef DEF_PARAMETER_READER
#define DEF_PARAMETER_READER

#include <string>
#include <fstream>
#include <iostream>
#include <streambuf>
#include <sstream>
#include <cassert>
#include <vector>
#include <map>
#include <optional>
#include <memory>
#include "tinyxml2.h"
#include <algorithm>



#include "utils.hpp"
#include "custom_structures.hpp"
#include "custom_exception.hpp"


/*
    This class contains all the methods that read the input XML paramater file. 
    All the data contain in this file are stored in specific structures. These 
    structures are then passed to the Solver class and Cell classes to set the 
    different parameters of the simulation.
*/

class parameter_reader{
   private:
        //The tinyxml2 object that will be used to parse the parameter file
        tinyxml2::XMLDocument xml_doc;

        //Class used for testing
        friend class test_parameter_reader;

   public:
        parameter_reader() = delete;                                    //default constructor
        parameter_reader(const parameter_reader& v) = delete;           //copy constructor
        parameter_reader(parameter_reader&& v) = delete;                //move constructor
        parameter_reader& operator=(const parameter_reader& v) = delete;//copy assignment operator
        parameter_reader& operator=(parameter_reader&& v) = delete;     //move assignment operator 
        
        //Instantiate the reader by providing the inpt mesh file
        parameter_reader(const std::string& paramter_file_path) noexcept(false);

        //Given a particular position in the XML file (tinyxml2::XMLElement), rand a given XML markup, this function
        //returns the value of the XML markup as a string. If the XML markup does not exist, it returns an empty
        std::optional<std::string> get_string_value(const tinyxml2::XMLElement* e,  const std::string XML_markup, const bool to_lower_case = false) const noexcept; 

        //Tries to select a section of the xl file, if the section is not found throw an exception
        tinyxml2::XMLElement* select_section(const std::string& section_name) noexcept(false);

        //Read the simulation parameters from the parameter file
        global_simulation_parameters read_numerical_parameters() noexcept(false); 

        //Load all the cell types and their face types
        std::vector<std::shared_ptr<cell_type_parameters>> read_biomechanical_parameters() noexcept(false);

        //Load the parameters of the given face type
        face_type_parameters read_face_type_parameters(
            const std::string& cell_type_name,
            const unsigned face_type_id,     
            tinyxml2::XMLElement* face_type_section
        ) noexcept(false);


        //Load the parameters of the given cell type but not its face type parameters
        std::shared_ptr<cell_type_parameters> read_cell_type_parameters(
            const unsigned cell_type_id,     
            tinyxml2::XMLElement* cell_type_section
        ) noexcept(false);


};

#endif