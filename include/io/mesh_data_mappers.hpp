#ifndef DEF_MESH_DATA_MAPPERS
#define DEF_MESH_DATA_MAPPERS

#include <vector>
#include <string>
#include <functional>


#include "utils.hpp"
#include "cell.hpp"
#include "face.hpp"

//---------------------------------------------------------------------------------------
//The following 2 structures are used to map the data of the cells and faces on the output VTK mesh
struct cell_data_mapper{

    //The name of the data that should be mapped, ex: "cell area"
    std::string value_name_;

    //The type of the data that should be mapped, ex: "double"
    std::string value_type_;

    //The function that extracts the data from the cell
    std::function<std::string(cell_ptr)>  value_extractor_;

    //default constructor
    cell_data_mapper() = default; 

    //trivial constructor
    cell_data_mapper(const std::string value_name, const std::string value_type, std::function<std::string(cell_ptr)> value_extractor):
        value_name_(value_name), value_type_(value_type), value_extractor_(value_extractor){};

};

struct face_data_mapper{

    //The name of the data that should be mapped, ex: "cell area"
    std::string value_name_;

    //The type of the data that should be mapped, ex: "double"
    std::string value_type_;

    //The function that extracts the data from the cell
    std::function<std::string(const face& f)>  value_extractor_;

    //default constructor
    face_data_mapper() = default; 

    //trivial constructor
    face_data_mapper(const std::string value_name, const std::string value_type, std::function<std::string(const face& f)> value_extractor):
        value_name_(value_name), value_type_(value_type), value_extractor_(value_extractor){};

};
//---------------------------------------------------------------------------------------




struct node_data_mapper{

    //The name of the data that should be mapped, ex: "cell area"
    std::string value_name_;

    //The type of the data that should be mapped, ex: "double"
    std::string value_type_;

    //The function that extracts the data from the cell
    std::function<std::string(const node& n)>  value_extractor_;

    //default constructor
    node_data_mapper() = default; 

    //trivial constructor
    node_data_mapper(const std::string value_name, const std::string value_type, std::function<std::string(const node& n)> value_extractor):
        value_name_(value_name), value_type_(value_type), value_extractor_(value_extractor){};

};
//---------------------------------------------------------------------------------------




#endif