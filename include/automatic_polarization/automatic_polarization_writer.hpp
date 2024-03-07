#ifndef DEF_AUTOMATIC_POLARIZATION_WRITER
#define DEF_AUTOMATIC_POLARIZATION_WRITER


#include <iostream>
#include <vector>
#include <array>
#include <omp.h>
#include <string>


#include "global_configuration.hpp"

#include "utils.hpp"
#include "custom_structures.hpp"
#include "custom_exception.hpp"

#include "mesh_writer.hpp"

#include "cell.hpp"
#include "face.hpp"
#include "node.hpp"
#include "uspg_3d.hpp"
#include "vec3.hpp"


#include "automatic_polarizer.hpp"


/*

This class is used to write in a vtk file the type of each voxel of the discretized space created 
by the automatic polarizer. 

The different types are: 
   -1: undefined
    0: boundary voxels, 
    1: cytoplasm,
    2: exterior space
    3: lumen,
    
*/



class automatic_polarization_writer{

        protected:

        public:
                automatic_polarization_writer() = default;                                         //default constructor
                automatic_polarization_writer(const automatic_polarization_writer& v) = delete;           //copy constructor
                automatic_polarization_writer(automatic_polarization_writer&& v) = delete;                //move constructor
                automatic_polarization_writer& operator=(const automatic_polarization_writer& v) = delete;//copy assignment operator
                automatic_polarization_writer& operator=(automatic_polarization_writer&& v) = default;    //move assignment operator 

                //Generate a cubic cell where min_corner contains {min_x, min_y, min_z} and max_corner contains {max_x, max_y, max_z}
                static cell_ptr generate_cubic_cell(const vec3& min_corner, const double edge_length, const unsigned cell_id) noexcept;

                //Write the whole automatic_polarizer discretization grid in a vtk file
                static void write(const std::string& output_path, const automatic_polarizer& polarizer) noexcept(false);

                //Create a cell for each voxel of the grid
                static std::vector<cell_ptr> convert_grid_voxels_to_cells(const uspg_3d<unsigned short>& grid) noexcept;

                //Write the type of each voxel in the file 
                static void write_voxel_type_data(
                    std::ofstream& polarization_data_file,
                    const uspg_3d<unsigned short>& grid
                );


};

#endif