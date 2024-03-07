#ifndef DEF_MESH_READER
#define DEF_MESH_READER


#include <string>
#include <fstream>
#include <iostream>
#include <streambuf>
#include <sstream>
#include <regex>
#include <cassert>
#include <vector>
#include <set>
#include <map>

#include "utils.hpp"
#include "custom_structures.hpp"
#include "custom_exception.hpp"


/*
    Class that contains all the functions to read an initial geometry file. The initial mesh should 
    be contained in a VTK file (.vtk) with ASCII encoding. 
*/

class mesh_reader{
   private:
      std::string file_content_;

   public:
      mesh_reader() = delete;                               //default constructor
      mesh_reader(const mesh_reader& v) = delete;           //copy constructor
      mesh_reader(mesh_reader&& v) = delete;                //move constructor
      mesh_reader& operator=(const mesh_reader& v) = delete;//copy assignment operator
      mesh_reader& operator=(mesh_reader&& v) = delete;     //move assignment operator 
      
      //Instantiate the reader by providing the inpt mesh file
      mesh_reader(const std::string& mesh_file_path, const bool verbose = true) noexcept(false);

      //Read the mesh file and return the data in a list of mesh objects
      std::vector<mesh> read() const noexcept(false);

      //Get the position of the mesh nodes
      std::vector<double> get_node_pos() const noexcept(false);

      //Return in a structure the faces of each cell
      std::vector<std::vector<unsigned>>  read_cell_faces() const noexcept(false);

      //Transform the loaded point data and face connectivity data into a vector of mesh objects
      static std::vector<mesh>  get_cell_mesh(
         const std::vector<double>& node_pos, 
         const std::vector<std::vector<unsigned>>& face_connectivity_lst_lst
      ) noexcept(false);

      //Find the cell_type_id array in the VTK file and extract the cell type of each cell
      std::vector<short>  get_cell_types() const noexcept(false);

};

#endif