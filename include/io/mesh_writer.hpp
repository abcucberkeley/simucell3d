#ifndef DEF_MESH_WRITER
#define DEF_MESH_WRITER


#include <string>
#include <fstream>
#include <iostream>
#include <cassert>
#include <vector>
#include <numeric>
#include <omp.h>


#include "cell.hpp"
#include "utils.hpp"
#include "custom_structures.hpp"
#include "custom_exception.hpp"
#include "mesh_data.hpp"

/*
    Write the cell meshes in .vtk files with ASCII form. The MeshWriter class is a static class, there is no need to instantiate it.
*/



class mesh_writer{


    public:

        mesh_writer() = delete;                               //default constructor
        mesh_writer(const mesh_writer& v) = delete;           //copy constructor
        mesh_writer(mesh_writer&& v) = delete;                //move constructor
        mesh_writer& operator=(const mesh_writer& v) = delete;//copy assignment operator
        mesh_writer& operator=(mesh_writer&& v) = delete;     //move assignment operator 

        //Write in parallel the cell data in a vtk file
        static void write( 
            const std::string& output_path_cell_data, 
            const std::string& output_path_face_data, 
            std::vector<cell_ptr>& cell_lst
        ) noexcept(false);

        //Write all the cell_data in a vtk file
        static void write_cell_data_file(   std::ofstream&  cell_data_file,
                                            const std::vector<mesh>& mesh_lst) noexcept(false);

        static void write_cell_data_file( std::ofstream&  cell_data_file,
                                        const std::vector<cell_ptr>& cell_lst,
                                        bool rebase_bool = true) noexcept(false);



        static void write_cell_data_file(const std::string& file_path,
                                        const std::vector<mesh>& mesh_lst
                                    ) noexcept(false);

        
       static void write_cell_data_file(const std::string& file_path,
                                        const std::vector<cell_ptr>& cell_lst,
                                        const bool rebase = true
                                    ) noexcept(false);

        static void write_face_data_file(const std::string& file_path,
                                const std::vector<cell_ptr>& cell_lst
                            ) noexcept(false);



        static void write_face_data_file(   std::ofstream&  face_data_file,
                                            const std::vector<cell_ptr>& cell_lst) noexcept(false);

        //This function can be called with a cell_lst or a mesh_lst
        template<typename T1>
        static void write_point_data(   std::ofstream& data_file, 
                                        const std::vector<T1>& object_lst,
                                        const std::function<std::vector<double>(const T1&)> node_pos_extractor  
                                    ) noexcept(false);




        static void write_cell_data(    std::ofstream& data_file, 
                                        const std::vector<mesh>& mesh_lst,
                                        const std::vector<size_t>& node_id_offset_lst
                                    ) noexcept(false);

        static void write_cell_data(    std::ofstream& data_file, 
                                        const std::vector<cell_ptr>& cell_lst,
                                        const std::vector<size_t>& node_id_offset_lst
                            ) noexcept(false);



        static void write_face_data(    std::ofstream& data_file, 
                                        const std::vector<cell_ptr>& cell_lst,
                                        const std::vector<size_t>& node_id_offset_lst
                            ) noexcept(false);



        static void add_cell_data_arrays_to_mesh(
            std::ofstream& data_file,                    
            const std::vector<cell_ptr>& cell_lst
        ) noexcept;


        static void add_face_data_arrays_to_mesh(
            std::ofstream& data_file,                    
            const std::vector<cell_ptr>& cell_lst
        ) noexcept;


        static void add_node_data_arrays_to_mesh(
            std::ofstream& data_file,                    
            const std::vector<cell_ptr>& cell_lst
        ) noexcept;


        static void add_cell_data_arrays_to_mesh(std::ofstream& data_file, const std::vector<mesh>& mesh_lst) noexcept{}
        static void add_face_data_arrays_to_mesh(std::ofstream& data_file, const std::vector<mesh>& mesh_lst) noexcept{}



};

#endif