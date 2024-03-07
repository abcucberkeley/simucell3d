#include "automatic_polarization_writer.hpp"





//-----------------------------------------------------------------------------------------------------------------
//Generate a cubic cell where min_corner is {min_x, min_y, min_z} 
cell_ptr automatic_polarization_writer::generate_cubic_cell(const vec3& min_corner, const double edge_length, const unsigned cell_id) noexcept{

    const auto [min_x, min_y, min_z] = min_corner.to_array();
    const double max_x = min_x + edge_length;
    const double max_y = min_y + edge_length;
    const double max_z = min_z + edge_length;

    std::vector<double> cell_node_pos_lst{
        min_x,min_y,min_z,  
        max_x,min_y,min_z,  
        max_x,min_y,max_z,
        min_x,min_y,max_z,  
        min_x,max_y,min_z,  
        max_x,max_y,min_z,
        min_x,max_y,max_z,  
        max_x,max_y,max_z
    };


    std::vector<unsigned> cell_face_connectivity{
        0, 1, 3, //0
        2, 3, 1, //1
        0, 4, 1, //2
        5, 1, 4, //3
        0, 3, 4, //4
        6, 4, 3, //5
        1, 5, 2, //6
        7, 2, 5, //7
        5, 4, 7, //8
        6, 7, 4, //9
        3, 2, 6, //10
        7, 6, 2  //11
    };

    return std::make_shared<cell>(cell_node_pos_lst, cell_face_connectivity, cell_id);
}
//-----------------------------------------------------------------------------------------------------------------






//-----------------------------------------------------------------------------------------------------------------
//Write the whole automatic_polarizer discretization grid in a vtk file
void automatic_polarization_writer::write(const std::string& output_path, const automatic_polarizer& polarizer) noexcept(false){

    //Starts by getting a reference of the grid
    const auto& grid = polarizer.get_grid();

    //Create the file where the discretization grid will be storeda
    std::ofstream polarization_data_file(output_path.c_str(), std::ofstream::out);

    //Throw an error if the file was not opened correctly
    if(!polarization_data_file){throw mesh_writer_exception("Unable to write file: " + output_path);}

    //Create a cell for each voxel of the grid
    std::vector<cell_ptr> voxel_lst = automatic_polarization_writer::convert_grid_voxels_to_cells(grid);

    //Write all the voxels in a vtk file
    mesh_writer::write_cell_data_file(
        polarization_data_file,
        voxel_lst,
        /*rebase=*/ false
    );

    //Write the type of each voxel in the file 
    write_voxel_type_data(polarization_data_file, grid);
    
    polarization_data_file.close();
}
//-----------------------------------------------------------------------------------------------------------------





//-----------------------------------------------------------------------------------------------------------------
//Create a cell for each voxel of the grid
std::vector<cell_ptr> automatic_polarization_writer::convert_grid_voxels_to_cells(const uspg_3d<unsigned short>& grid) noexcept{

    //Get the grid voxel edge length
    const double voxel_size = grid.get_voxel_size();

    //Get the voxel origin
    const auto [min_x, min_y, min_z] = grid.get_min_corner();

    //Get the number of voxels in each direction
    const auto [nb_voxels_x, nb_voxels_y, nb_voxels_z] = grid.get_nb_voxels();

    //Store the created cells in this vector
    std::vector<cell_ptr> voxel_lst;

    for(unsigned voxel_x_index = 0; voxel_x_index < nb_voxels_x; voxel_x_index++){
    for(unsigned voxel_y_index = 0; voxel_y_index < nb_voxels_y; voxel_y_index++){
    for(unsigned voxel_z_index = 0; voxel_z_index < nb_voxels_z; voxel_z_index++){

        //Compute the origin of the voxel
        const vec3 voxel_min_coord(
            min_x + voxel_x_index * voxel_size,
            min_y + voxel_y_index * voxel_size,
            min_z + voxel_z_index * voxel_size
        );

        //Create the cell corresponding to the voxel
        voxel_lst.push_back(
            automatic_polarization_writer::generate_cubic_cell(voxel_min_coord, voxel_size, voxel_lst.size())
        );
    }
    }
    }
    return voxel_lst;
}
//-----------------------------------------------------------------------------------------------------------------


//-----------------------------------------------------------------------------------------------------------------
//Write the type of each voxel in the file 
void automatic_polarization_writer::write_voxel_type_data(
    std::ofstream& polarization_data_file,
    const uspg_3d<unsigned short>& grid
){

    //Get the number of voxels in each direction
    const auto [nb_voxels_x, nb_voxels_y, nb_voxels_z] = grid.get_nb_voxels();
    const size_t nb_voxels_total = nb_voxels_x * nb_voxels_y * nb_voxels_z;


    polarization_data_file << "\nCELL_DATA " + std::to_string(nb_voxels_total) + "\n";
    polarization_data_file << "FIELD FieldData 1\n";
    polarization_data_file << "voxel_type_id 1 " + std::to_string(nb_voxels_total) + " int\n";

    size_t iteration = 0;
    for(unsigned voxel_x_index = 0; voxel_x_index < nb_voxels_x; voxel_x_index++){
    for(unsigned voxel_y_index = 0; voxel_y_index < nb_voxels_y; voxel_y_index++){
    for(unsigned voxel_z_index = 0; voxel_z_index < nb_voxels_z; voxel_z_index++){

        //Get the type of the voxel
        const std::optional<unsigned short> voxel_type = grid.get_voxel_content(voxel_x_index, voxel_y_index, voxel_z_index);
        const std::string voxel_value = voxel_type.has_value() ? std::to_string(voxel_type.value()) : "-1";

        const std::string sep = ((iteration+1)%9==0) ? "\n" : " ";
        polarization_data_file << voxel_value + sep;

        iteration++;

    }
    }
    }

}
//-----------------------------------------------------------------------------------------------------------------
