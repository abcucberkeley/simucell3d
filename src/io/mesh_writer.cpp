#include "mesh_writer.hpp"



//-----------------------------------------------------------------------------------------------------------------
//Write in parallel the cell data in a vtk file
void mesh_writer::write( 
    const std::string& output_path_cell_data, 
    const std::string& output_path_face_data, 
    std::vector<cell_ptr>& cell_lst) noexcept(false){

    assert(cell_lst.size());
    

    //Wrap the call to cell::rebase() since it might throw an exception
    const std::function<void(cell_ptr)> rebase_func  = [=](cell_ptr c) -> void {c->rebase();};

    //Call the rebase function in parallel, rethrow the exception if any is thrown
    parallel_exception_handler(cell_lst, rebase_func);


    //Create an exception pointer
    std::exception_ptr e_ptr;

    //Write the cell and face mesh files simultaneously
    #pragma omp parallel sections
    {
        #pragma omp section
        { 
            try{
                //Create the file 
                std::ofstream cell_data_file(output_path_cell_data.c_str(), std::ofstream::out);

                //Throw an error if the file was not opened correctly
                if(!cell_data_file){throw mesh_writer_exception("Unable to write file: " + output_path_cell_data);}

                //Write the cell surface meshes in the vtk file
                mesh_writer::write_cell_data_file(cell_data_file, cell_lst, false); 

                //Add the cell data arrays to the vtk file (cell_volume, cell_area .. etc)
                mesh_writer::add_cell_data_arrays_to_mesh(cell_data_file, cell_lst);

                //Close the file
                cell_data_file.close();
            }

            //Capture the exception thrown by write_cell_data_file
            catch(...){
            
                #pragma omp critical
                {
                    //Store the exception
                    e_ptr = std::current_exception(); 
                }
            }      
        }
        #pragma omp section
        { 
            try{

                //Create the file 
                std::ofstream face_data_file(output_path_face_data.c_str(), std::ofstream::out);

                //Throw an error if the file was not opened correctly
                if(!face_data_file){throw mesh_writer_exception("Unable to write file: " + output_path_face_data);}

                //Write the cell surface meshes in the vtk file
                mesh_writer::write_face_data_file(face_data_file, cell_lst);

                //Add the face data arrays to the vtk file (face_area, face_surface_tension .. etc)
                mesh_writer::add_face_data_arrays_to_mesh(face_data_file, cell_lst);
                mesh_writer::add_node_data_arrays_to_mesh(face_data_file, cell_lst);
                
                //Close the file 
                face_data_file.close();
            }

            //Capture the exception thrown by write_face_data_file
            catch(...){
            
                #pragma omp critical
                {
                    //Store the exception
                    e_ptr = std::current_exception(); 
                }
            }           
        }
    }

    //If an exception was thrown during the execution of the parallel loop, rethrow it
    if (e_ptr) std::rethrow_exception(e_ptr);

}

//-----------------------------------------------------------------------------------------------------------------


//-----------------------------------------------------------------------------------------------------------------
void mesh_writer::write_cell_data_file(const std::string& file_path,
                                const std::vector<mesh>& mesh_lst
                            ) noexcept(false){

    //Create the file 
    std::ofstream cell_data_file(file_path.c_str(), std::ofstream::out);

    //Throw an error if the file was not opened correctly
    if(!cell_data_file){throw mesh_writer_exception("Unable to write file: " + file_path);}

    mesh_writer::write_cell_data_file(cell_data_file, mesh_lst);

    cell_data_file.close();

}
//-----------------------------------------------------------------------------------------------------------------


//-----------------------------------------------------------------------------------------------------------------
void mesh_writer::write_face_data_file(const std::string& file_path,
                        const std::vector<cell_ptr>& cell_lst
                    ) noexcept(false){

    std::ofstream face_data_file(file_path.c_str(), std::ofstream::out);


    //Throw an error if the file was not opened correctly
    if(!face_data_file){throw mesh_writer_exception("Unable to write file: " + file_path);}

    mesh_writer::write_cell_data_file(face_data_file, cell_lst);

    face_data_file.close();
           
}
//-----------------------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------------------
void mesh_writer::write_cell_data_file(const std::string& file_path,
                                const std::vector<cell_ptr>& cell_lst, 
                                const bool rebase /*=true*/
                            ) noexcept(false){

        //Create the file 
        std::ofstream cell_data_file(file_path.c_str(), std::ofstream::out);

        //Throw an error if the file was not opened correctly
        if(!cell_data_file){throw mesh_writer_exception("Unable to write file: " + file_path);}

        mesh_writer::write_cell_data_file(cell_data_file, cell_lst, rebase);

        cell_data_file.close();
}
//-----------------------------------------------------------------------------------------------------------------




//-----------------------------------------------------------------------------------------------------------------
//Write all the cell_data in a vtk file
void mesh_writer::write_cell_data_file( std::ofstream&  cell_data_file,
                                        const std::vector<mesh>& mesh_lst) noexcept(false){

                                            
    assert(mesh_lst.size());

    //Write the header
    const std::string header = "# vtk DataFile Version 4.2\nvtk output\nASCII\nDATASET UNSTRUCTURED_GRID\n";
    cell_data_file << header;

    //Fnctionn that returns the number of nodes in a mesh
    std::function<size_t(const mesh&)> get_mesh_nb_nodes = [](const mesh& m) -> size_t {return m.node_pos_lst.size() / 3;};

    //Get a serie of partial sum of the number of nodes in the different meshes
    std::vector<size_t> node_id_offset_lst = partial_sum_vector(mesh_lst, get_mesh_nb_nodes);
    node_id_offset_lst.insert(node_id_offset_lst.begin(), 0);

    //Get the total nb of nodes in the mesh
    const size_t nb_nodes = node_id_offset_lst.back();
    cell_data_file << "POINTS "+std::to_string(nb_nodes)+" float\n";

    //Function that returns the position of the nodes in a mesh
    std::function<std::vector<double>(const mesh&)> node_pos_extractor = [](const mesh& m){return m.node_pos_lst;};

    //Write the node coordinates in the output mesh file
    mesh_writer::write_point_data(cell_data_file, mesh_lst, node_pos_extractor);

    //Write the connectivity of the faces in the mesh file
    write_cell_data(cell_data_file, mesh_lst, node_id_offset_lst);

}
//-----------------------------------------------------------------------------------------------------------------







//-----------------------------------------------------------------------------------------------------------------
//Write all the cell_data in a vtk file
void mesh_writer::write_cell_data_file( std::ofstream& cell_data_file, 
                                        const std::vector<cell_ptr>& cell_lst, 
                                        bool rebase_bool) noexcept(false){
    assert(cell_lst.size());
    

    //Call the rebase function on all the cells to make sure they use all of their available nodes and faces
    if(rebase_bool) std::for_each(cell_lst.begin(), cell_lst.end(), [&cell_lst](const cell_ptr c)-> void {c->rebase();});
    

    //Write the header
    const std::string header = "# vtk DataFile Version 4.2\nvtk output\nASCII\nDATASET UNSTRUCTURED_GRID\n";
    cell_data_file << header;

    //Function that returns the number of nodes in a cell

    std::function<size_t(const mesh&)> get_mesh_nb_nodes = [](const mesh& m) -> size_t {return m.node_pos_lst.size() / 3;};

    std::function<size_t(const cell_ptr&)> get_cell_nb_nodes = [](const cell_ptr& c) -> size_t {return c->get_node_lst().size();};

    //Get a serie of partial sum of the number of nodes in the different meshes
    std::vector<size_t> node_id_offset_lst = partial_sum_vector(cell_lst, get_cell_nb_nodes);
    node_id_offset_lst.insert(node_id_offset_lst.begin(), 0);

    //Get the total nb of nodes in the mesh
    const size_t total_nb_nodes = node_id_offset_lst.back();
    cell_data_file << "POINTS "+std::to_string(total_nb_nodes)+" float\n";

    //Function that returns the position of the nodes in a cell
    std::function<std::vector<double>(const cell_ptr&)> node_pos_extractor =  [](const cell_ptr& c){return c->get_flat_node_coord_lst();};
    

    //Write the node coordinates in the output mesh file
    mesh_writer::write_point_data(cell_data_file, cell_lst, node_pos_extractor);

    //Write the connectivity of the faces in the mesh file
    write_cell_data(cell_data_file, cell_lst, node_id_offset_lst);


    

}
//-----------------------------------------------------------------------------------------------------------------






//-----------------------------------------------------------------------------------------------------------------
//Write all the cell_data in a vtk file
template<typename T1>
void mesh_writer::write_point_data(   std::ofstream& data_file, 
                                        const std::vector<T1>& object_lst,
                                        const std::function<std::vector<double>(const T1&)> node_pos_extractor  
                                    ) noexcept(false){


    assert(object_lst.size());
    


    //Keep track of the number of coordinates added in the output file
    size_t i = 0;

    //Loop over the meshes or the cells
    for(const T1& object: object_lst){

        //Get all the position of all the nodes in a given mesh or cell
        //Convert each coordinate to a string in scientific notation
        for(const double coord: node_pos_extractor(object)){
            i++;

            //Store the string in the output mesh file
            std::string sep = (i % 9 == 0) ? "\n" : " ";

            if(!std::isfinite(coord)){throw mesh_writer_exception("NaN found in the point coordinates. The simulation is unstable. Reducing the time step might help");}

            data_file << format_number(coord, "%.4e") + sep;
            
        }
    }
}
//-----------------------------------------------------------------------------------------------------------------


//-----------------------------------------------------------------------------------------------------------------
//Write the conectivity of each mes object in the VTK file
void mesh_writer::write_cell_data(std::ofstream& data_file, 
                                        const std::vector<mesh>& mesh_lst,
                                        const std::vector<size_t>& node_id_offset_lst
                                    ) noexcept(false){

    //First, we calculate the number of integers required to store all the connectivity data 
    //of the cell surfaces
    assert(mesh_lst.size() > 0);

    size_t nb_cells = std::count_if(mesh_lst.begin(), mesh_lst.end(), [](const mesh& m){return m.face_point_ids.size() > 0;});

    //Store the number of integers required to store the mesh of each cell
    std::vector<size_t> cell_int_size;
    cell_int_size.reserve(nb_cells);

    //Loop over the meshes structures or the cell objects
    for(const mesh& m: mesh_lst){

        //Check that there are faces left in the mesh
        if(m.face_point_ids.size() == 0) continue;

        //One integer per face will be necessary plus one integer to indicate the nb of faces
        size_t nb_int_cell = static_cast<size_t>(1) + m.face_point_ids.size();

        //Sum the number of nodes of each face
        nb_int_cell += std::accumulate(m.face_point_ids.begin(), m.face_point_ids.end(), 0, 
            [](size_t i, const auto& f) -> size_t {return i + f.size();}
        );

        //Store the nb of integer required to specify the cell mesh
        cell_int_size.push_back(nb_int_cell);
    }

    //Compute the number of integers required to represent the whole tissue
    const size_t nb_integer_tissue = (nb_cells == 0 ) ? 0 : std::accumulate(cell_int_size.begin(), cell_int_size.end(), mesh_lst.size(), std::plus<size_t>());
    
    //Write the header of the cell section
    data_file <<"\n\nCELLS "+ std::to_string(nb_cells) +" "+ std::to_string(nb_integer_tissue)+"\n";

    //Reloop over the cells
    for(size_t i = 0; i < nb_cells; i++){

        const mesh& m = mesh_lst[i];

        //Check that there are faces left in the mesh
        if(m.face_point_ids.size() == 0) continue;

        //Get the number of faces in the cell
        const size_t cell_nb_faces = m.face_point_ids.size();

        //Convert from local node ids to global node ids
        const size_t node_id_offset = node_id_offset_lst[i];

        //Starts to write the cell information
        data_file << std::to_string(cell_int_size[i]) +" "+ std::to_string(cell_nb_faces) + " ";

        //Loop over the faces and save them in the file
        for(const std::vector<unsigned>& face: m.face_point_ids){
            data_file << std::to_string(face.size()) +" ";
            
            for(auto node_id: face){
                data_file << std::to_string(node_id + node_id_offset) + " ";
            }
        }
        data_file << "\n";
    }

    //Indicate that each cell is a polyhedron
    data_file << "\nCELL_TYPES " + std::to_string(nb_cells) + "\n";
    for(int i = 0; i < nb_cells; i++){data_file << "42\n";}

    


}
//-----------------------------------------------------------------------------------------------------------------







//-----------------------------------------------------------------------------------------------------------------
//Write the conectivity of each mes object in the VTK file
void mesh_writer::write_cell_data(std::ofstream& data_file, 
                                        const std::vector<cell_ptr>& cell_lst,
                                        const std::vector<size_t>& node_id_offset_lst
                                    ) noexcept(false){
    assert(cell_lst.size());
    

    //First, we calculate the number of integers required to store all the connectivity data 
    //of the cell surfaces
    size_t nb_cells = cell_lst.size();
    assert(nb_cells > 0);

    //Store the number of integers required to store the mesh of each cell
    std::vector<size_t> cell_int_size;
    cell_int_size.reserve(nb_cells);


    //Loop over the meshes structures or the cell objects
    for(const cell_ptr c: cell_lst){


        //Check that there are faces left in the mesh
        assert(c->get_node_lst().size() > 0);
        assert(c->get_face_lst().size() > 0);

        //One integer per face will be necessary plus one integer to indicate the nb of faces
        size_t nb_int_cell = 1 + c->get_face_lst().size() * 4;
        
        //Store the nb of integer required to specify the cell mesh
        cell_int_size.push_back(nb_int_cell);
    }

    //Compute the number of integers required to represent the whole tissue
    const size_t nb_integer_tissue = std::accumulate(cell_int_size.begin(), cell_int_size.end(), cell_lst.size(), std::plus<size_t>());
    
    //Write the header of the cell section
    data_file <<"\n\nCELLS "+ std::to_string(nb_cells) +" "+ std::to_string(nb_integer_tissue)+"\n";

    //Reloop over the cells
    for(size_t i = 0; i < cell_lst.size(); i++){
        
        const cell_ptr c = cell_lst[i];

        //Get the number of faces in the cell
        const size_t cell_nb_faces = c->get_face_lst().size();

        //Convert from local node ids to global node ids
        const size_t node_id_offset = node_id_offset_lst[i];

        //Starts to write the cell information
        data_file << std::to_string(cell_int_size[i]) +" "+ std::to_string(cell_nb_faces) + " ";

        //Loop over the faces and save them in the file
        for(const face& f: c->get_face_lst()){
            data_file << std::to_string(3) +" ";
            
            for(auto node_id: f.get_node_ids()){
                data_file << std::to_string(node_id + node_id_offset) + " ";
            }
        }
        data_file << "\n";
    }

    //Indicate that each cell is a polyhedron
    data_file << "\nCELL_TYPES " + std::to_string(nb_cells) + "\n";
    for(int i = 0; i < nb_cells; i++){data_file << "42\n";}
}
//-----------------------------------------------------------------------------------------------------------------



//-----------------------------------------------------------------------------------------------------------------
//This method add all the cell data given by the extraction functions in /include/io/mesh_data.hpp cell_data_mapper_lst
//to the vtk file
void mesh_writer::add_cell_data_arrays_to_mesh(
    std::ofstream& data_file,                    
    const std::vector<cell_ptr>& cell_lst
) noexcept{

    //Loop over the mesh mapping function defined "output_mes_data.hpp"
    if(!cell_data_mapper_lst.empty()){

        size_t nb_cells = cell_lst.size();

        data_file << "\nCELL_DATA " + std::to_string(nb_cells) + "\n";
        data_file << "FIELD FieldData " + std::to_string(cell_data_mapper_lst.size());

        //Loop over the mapping functions
        for(const auto& data_mapper: cell_data_mapper_lst){

            data_file << "\n" + data_mapper.value_name_ + " 1 " + std::to_string(nb_cells) + " " + data_mapper.value_type_ + "\n";
  
            //Loop over the cells
            for(size_t i = 0; i < cell_lst.size(); i++){
                cell_ptr c = cell_lst[i];

                const std::string sep = ((i+1)%9==0) ? "\n" : " ";
                data_file << data_mapper.value_extractor_(c) + sep;
            }
        }
    }
}
//-----------------------------------------------------------------------------------------------------------------




//-----------------------------------------------------------------------------------------------------------------
//Write the data associated with each face (area, adhesion_energy, repulsion energy .. etc)
void mesh_writer::write_face_data_file(   std::ofstream& face_data_file,   
                                          const std::vector<cell_ptr>& cell_lst) noexcept(false){

    assert(cell_lst.size());
    //The beginning of the function is the same as the write_cell_data function

    //Write the header
    const std::string header = "# vtk DataFile Version 4.2\nvtk output\nASCII\nDATASET UNSTRUCTURED_GRID\n";
    face_data_file << header;

    //Function that returns the number of nodes in a cell

    std::function<size_t(const mesh&)> get_mesh_nb_nodes = [](const mesh& m) -> size_t {return m.node_pos_lst.size() / 3;};

    std::function<size_t(const cell_ptr&)> get_cell_nb_nodes = [](const cell_ptr& c) -> size_t {return c->get_node_lst().size();};

    //Get a serie of partial sum of the number of nodes in the different meshes
    std::vector<size_t> node_id_offset_lst = partial_sum_vector(cell_lst, get_cell_nb_nodes);
    node_id_offset_lst.insert(node_id_offset_lst.begin(), 0);

    //Get the total nb of nodes in the mesh
    const size_t total_nb_nodes = node_id_offset_lst.back();
    face_data_file << "POINTS "+std::to_string(total_nb_nodes)+" float\n";

    //Function that returns the position of the nodes in a cell
    std::function<std::vector<double>(const cell_ptr&)> node_pos_extractor =  [](const cell_ptr& c){return c->get_flat_node_coord_lst();};
    

    //Write the node coordinates in the output mesh file
    mesh_writer::write_point_data(face_data_file, cell_lst, node_pos_extractor);

    //Write the connectivity of the faces in the mesh file
    write_face_data(face_data_file, cell_lst, node_id_offset_lst);



}
//-----------------------------------------------------------------------------------------------------------------




//-----------------------------------------------------------------------------------------------------------------
void mesh_writer::write_face_data(    std::ofstream& data_file, 
                                const std::vector<cell_ptr>& cell_lst,
                                const std::vector<size_t>& node_id_offset_lst
                    ) noexcept(false){


    //First, we calculate the number of integers required to store all the connectivity data 
    //of the cell surfaces
    assert(cell_lst.size() > 0);

    //Calculate the number of integer the vtk file will have to store
    const size_t nb_integer_mesh = std::accumulate(cell_lst.begin(), cell_lst.end(), 0, [](size_t acc, const cell_ptr c){return acc + c->get_face_lst().size() * 4;});
    
    //Calculate the number of faces in the tissue
    const size_t nb_faces = std::accumulate(cell_lst.begin(), cell_lst.end(), 0, [](size_t acc, const cell_ptr c){return acc + c->get_face_lst().size();});
    const size_t nb_nodes = std::accumulate(cell_lst.begin(), cell_lst.end(), 0, [](size_t acc, const cell_ptr c){return acc + c->get_nb_of_nodes();});


    //Save the cell information
    data_file <<"\nCELLS "+ std::to_string(nb_faces) +" "+ std::to_string(nb_integer_mesh)+"\n";

    //Loop over the cells
    for(size_t i = 0; i < cell_lst.size(); i++){

        //Get the cell and the offset to convert the local node ids to global node ids
        cell_ptr c =  cell_lst[i];
        size_t node_id_offset = node_id_offset_lst[i];
       
        //Loop over the faces
        for(const face& f: c->get_face_lst()){
            if(!f.is_used()) continue;

            //Get the local ids of the 3 nodes composing the face
            auto [local_id_n1, local_id_n2, local_id_n3] = f.get_node_ids();

            //Convert the local ids to global ids
            size_t global_id_n1 = local_id_n1 + node_id_offset;
            size_t global_id_n2 = local_id_n2 + node_id_offset;
            size_t global_id_n3 = local_id_n3 + node_id_offset;

            //Write the face information
            data_file << "3 " << global_id_n1 << " " << global_id_n2 << " " << global_id_n3 << "\n";
        }
    }

    //Store the VTK cell types
    data_file << "\nCELL_TYPES " + std::to_string(nb_faces) + "\n";
    for(size_t i = 0; i < nb_faces; i++){data_file << "7\n";  } 
}

//-----------------------------------------------------------------------------------------------------------------



//-----------------------------------------------------------------------------------------------------------------
//This method add all the cell data given by the extraction functions in /include/io/mesh_data.hpp cell_data_mapper_lst
//to the vtk file
void mesh_writer::add_face_data_arrays_to_mesh(
    std::ofstream& data_file,                    
    const std::vector<cell_ptr>& cell_lst
) noexcept{

    //Calculate the number of faces in the tissue
    const size_t nb_faces = std::accumulate(cell_lst.begin(), cell_lst.end(), 0, [](size_t acc, const cell_ptr c){return acc + c->get_face_lst().size();});

    //Check that all the cells have a type otherwise do not write the cell data
    bool all_have_type = std::all_of(cell_lst.begin(), cell_lst.end(), [](const cell_ptr c){return c->get_cell_type() != nullptr;});

    //Loop over the mesh mapping function defined "output_mes_data.hpp"
    if(!face_data_mapper_lst.empty() && all_have_type){

        data_file << "\nCELL_DATA "+std::to_string(nb_faces) + "\nFIELD FieldData " + std::to_string(face_data_mapper_lst.size());

        //Loop over the mapping functions
        for(const auto& data_mapper: face_data_mapper_lst){

            data_file << "\n" + data_mapper.value_name_ + " 1 " + std::to_string(nb_faces) + " " + data_mapper.value_type_ + "\n";
  
            //Loop over the cells
            size_t i = 0;
            for(cell_ptr c: cell_lst){

                //Loop over the faces of the cells
                for(const face& f: c->get_face_lst()){

                    //Get the value of the data to be written
                    const std::string sep = ((i+1)%9==0) ? "\n" : " ";
                    data_file << data_mapper.value_extractor_(f) + sep;
                    i++;
                }
            }
        }

        //If the normals should be written in the mesh filee
        #if WRITE_NORMALS_IN_OUTPUT_MESH_FILE 

            //Store the face information
            data_file << "\nVECTORS face_normals double\n";

            size_t i = 0;
            for(cell_ptr c: cell_lst){
        
                //Loop over the faces
                for(const face& f: c->get_face_lst()){

                    //Get the local ids of the 3 nodes composing the face
                    auto [normal_dx, normal_dy, normal_dz] = f.get_normal().to_array();
                    data_file << format_number(normal_dx, "%.3e") +" "+ format_number(normal_dy, "%.3e") +" "+ format_number(normal_dz, "%.3e") +" ";

                    //Line return every 3 faces
                    if((i+1)%3 ==0){data_file<< "\n";}
                    i++;
                }
            }

        #endif
    }
}
//-----------------------------------------------------------------------------------------------------------------





//--------------------------------------------------------------------------------------------
//This method add to the mesh all the point data given by the extraction functions in /include/io/mesh_data.hpp node_data_mapper_lst
void mesh_writer::add_node_data_arrays_to_mesh(
    std::ofstream& data_file,                    
    const std::vector<cell_ptr>& cell_lst
) noexcept{

    //Calculate the number of faces in the tissue
    const size_t nb_nodes = std::accumulate(cell_lst.begin(), cell_lst.end(), 0, [](size_t acc, const cell_ptr c){return acc + c->get_nb_of_nodes();});

    //Loop over the mesh mapping function defined "output_mes_data.hpp"
    if(!node_data_mapper_lst.empty()){

        //Write the point data
        data_file << "\nPOINT_DATA " + std::to_string(nb_nodes) + "\n";
        data_file << "FIELD FieldData " + std::to_string(node_data_mapper_lst.size()) + "\n";

        //Loop over the mapping functions
        for(const auto& data_mapper: node_data_mapper_lst){

            data_file << "\n" + data_mapper.value_name_ + " 1 " + std::to_string(nb_nodes) + " " + data_mapper.value_type_ + "\n";

            //Loop over the cells
            size_t i = 0;
            for(cell_ptr c: cell_lst){

                //Loop over the nodes of the cells
                for(const node& n: c->get_node_lst()){

                    //Get the value of the data to be written
                    const std::string sep = ((++i)%9==0) ? "\n" : " ";
                    data_file << data_mapper.value_extractor_(n) + sep;
                }
            }
        }




    //If the normals should be written in the mesh filee
    #if WRITE_NORMALS_IN_OUTPUT_MESH_FILE && CONTACT_MODEL_INDEX == 1

        //Save the nodes normals
        data_file << "\nVECTORS node_normals double\n";
        
        size_t i = 0;
        for(cell_ptr c: cell_lst){

            //Loop over the faces of the cells
            for(const node& n: c->get_node_lst()){

                const vec3 node_normal = n.get_normal();

                const std::string sep = ((++i)%3==0) ? "\n" : " ";
                data_file << format_number(node_normal.dx(), "%.3e") +" "+ format_number(node_normal.dy(), "%.3e") +" "+ format_number(node_normal.dz(), "%.3e") + sep;
            }
        }

    #endif

    } 
}
//-----------------------------------------------------------------------------------------------------------------