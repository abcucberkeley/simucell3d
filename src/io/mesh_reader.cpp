#include "mesh_reader.hpp"
          


//-------------------------------------------------------------------------
mesh_reader::mesh_reader(const std::string& mesh_file_path, const bool verbose) noexcept(false){

    //Read the input mesh file and store it in a string
    std::ifstream stream_absolute(mesh_file_path);

    //All the file content will be passed to this stream buffer 
    std::stringstream buffer;
    

    //If the given path does not exist
    if(stream_absolute.good()) {
        buffer << stream_absolute.rdbuf();
    }
    //Else try to load the file with a relative path
    else{
        std::string relative_path = std::string(PROJECT_SOURCE_DIR) + "/" + mesh_file_path;
        std::ifstream stream_relative(relative_path);

        if(verbose) std::cerr << "WARNING: the absolute mesh file path: " << mesh_file_path << " was not found" << std::endl;
        if(verbose) std::cerr << "The program will try instead the relative path: " << relative_path << std::endl << std::endl;

        if(stream_relative.good()) {
            buffer << stream_relative.rdbuf();
        }
        else{
            throw mesh_reader_exception("ERROR: input mesh file " + mesh_file_path + " not found"); 
        }
    }

    //Load the file content in a string
    file_content_ = buffer.str();

    //Get the version of the vtk file
    std::smatch match_1;
    std::regex rgx_1(R"(# vtk DataFile Version (\d*\.?\d*))");
    if(!std::regex_search(file_content_, match_1, rgx_1)){
        throw mesh_reader_exception("ERROR: header with .vtk file version not found"); 
    }

    //Make sure its the correct version
    double file_version = std::stod(match_1[1].str());
    if(file_version != 4.2){
        if(verbose) std::cerr << "WARNING: input vtk mesh file has a version of " << file_version << " instead of 4.2. This may cause problems" << std::endl;
    }
}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
//The main function returning the mesh_lst
std::vector<mesh> mesh_reader::read() const noexcept(false){

    //Get the position of the nodes
    const std::vector<double> node_pos = get_node_pos();

    //Get the connectivity data of the faces
    const std::vector<std::vector<unsigned>>  face_connectivity  = read_cell_faces();

    //Use this 2 dataset to create the cell objects
    return mesh_reader::get_cell_mesh(node_pos, face_connectivity); 
}


//-------------------------------------------------------------------------



//-------------------------------------------------------------------------
//Get the position of the mesh nodes
std::vector<double> mesh_reader::get_node_pos() const noexcept(false){


    //Regular expression to retrieve the number of points in the mesh
    std::smatch match_1;
    std::regex rgx_1(R"(POINTS ([0-9]+) ([a-z]+))");

    //Return an error if the number of points is not found
    if(!std::regex_search(file_content_, match_1, rgx_1)){
        throw mesh_reader_exception("ERROR: input mesh file, line with number of points not found"); 
    }

    //Extract the number of points
    int nb_nodes = std::stoi(match_1[1].str());

    //Make sure they have the correct format
    const std::string coord_type = match_1[2].str();
    if(coord_type != "float" && coord_type != "double"){
        throw mesh_reader_exception("ERROR: input mesh file, nodes have incorrect format: " + coord_type); 
    }
    
    //The position in the text where the node coordinates start
    auto coord_pos_start = match_1.position(0) + match_1.length(0);


    //Determines where the node coordinates terminates in the text
    std::smatch match_2;
    std::regex  rgx_2(R"(([A-Za-z_]){2,})");
    std::string sub_string = file_content_.substr(coord_pos_start);

    //If only the points are stored in the file
    std::string point_coord_text;
    if(!std::regex_search(sub_string, match_2, rgx_2)){
        point_coord_text = sub_string;
    }
    else{
        auto coord_pos_end =  match_2.position(0);
        point_coord_text = sub_string.substr(0,coord_pos_end);
    }


    //Store all the node positions in this vector
    std::vector<double> node_pos;
    node_pos.reserve(nb_nodes * 3);

    //Very complex regular expression that matches all decimal numbers
    std::regex rgx_3(R"(([-\+]?[\d.]+(?:[e|E][-\+]?\d+)?))");
    std::smatch matches_3;
    std::string::const_iterator search_start(point_coord_text.cbegin() );

    //Read the node coordinates and save them
    while(std::regex_search(search_start, point_coord_text.cend(), matches_3, rgx_3)){
        
        //Make sure the position value is in a valid format
        double node_pos_value;

        //Throw an error if the node coordinate could not be converetd
        try{
            node_pos_value = std::stod(matches_3[0]); 
        }
        catch(const std::invalid_argument& e){
            throw mesh_reader_exception("ERROR: input mesh file, impossible node value conversion: " + std::string(matches_3[0])); 
        }

        //Make sure the number is finite and not NA
        if(!std::isfinite(node_pos_value)){
            throw mesh_reader_exception("ERROR: input mesh file, node value is not finite: " + std::string(matches_3[0])); 
        }

        //If everything is okay save the node value
        node_pos.push_back(node_pos_value);

        //Restart the number search at the next position
        search_start = matches_3.suffix().first;
    }

    //Check that the correct number of nodes have been loaded
            //Make sure the number is finite and not NA
        if((int) node_pos.size() / 3 != nb_nodes){
            throw mesh_reader_exception("ERROR mesh_reader: not all the nodes were loaded. Some node coordinates probably are NA values"); 
        }

    //Return the loaded node coordinates
    return node_pos;

}
//-------------------------------------------------------------------------



//-------------------------------------------------------------------------
//Return in a 2D vector all the connectivity information of the cell faces
std::vector<std::vector<unsigned>> mesh_reader::read_cell_faces() const noexcept(false){

    //--------------------------------------------------------------------------------------------------
    //This first section makes sure all the cellls are stored as polyhedron #format 42
    std::smatch match_0;
    std::regex rgx_0(R"(CELL_TYPES ([0-9]+))");

    if(!std::regex_search(file_content_, match_0, rgx_0)){
        throw mesh_reader_exception("ERROR: input mesh file, could not find CELL_TYPES section"); 
    }
    int nb_cells = std::stoi(match_0[1].str());

    //Find the starting position of the cell type section
    auto type_pos_start = match_0.position(0) + match_0.length(0);
    std::string sub_string = file_content_.substr(type_pos_start);



    //Try to determine the end position in the text of the CELL section 
    std::smatch match_1;
    std::regex  rgx_1(R"(([A-Z]))");

    //If only the points are stored in the file
    std::string cell_type_section_str;
    if(!std::regex_search(sub_string, match_1, rgx_1)){
        cell_type_section_str = sub_string;
    }
    else{
        auto coord_pos_end =  match_1.position(0);
        cell_type_section_str = sub_string.substr(0,coord_pos_end);
    }


    //Loop over the cell type of each cell 
    std::regex rgx_2(R"([0-9]+)");
    std::smatch matches_2;
    std::string::const_iterator search_start(cell_type_section_str.cbegin() );

    //Use this integer to check that the CELL_TYPE section is complete
    int i = 0; 


    //Read the node coordinates and save them
    while(std::regex_search(search_start, cell_type_section_str.cend(), matches_2, rgx_2)){
        
        //Make sure the position value is in a valid format
        int cell_type;

        //Throw an error if the node coordinate could not be converetd
        try{
            cell_type = std::stoi(matches_2[0]); 
        }
        catch(const std::invalid_argument& e){
            throw mesh_reader_exception("ERROR: input mesh file, impossible node value conversion: " + std::string(matches_2[0])); 
        }

        //Make sure the number is finite and not NA
        if(cell_type != 42){
            throw mesh_reader_exception("ERROR: cell is not a polyhedron (type 42) but has type: " + std::string(matches_2[0])); 
        }

        i++;

        //Restart the number search at the next position
        search_start = matches_2.suffix().first;
    }


    if(i != nb_cells){
        throw mesh_reader_exception("ERROR: not all cells have a CELL_TYPE defined. nb_cells = " + 
            std::to_string(nb_cells) + ", nb counted "+ std::to_string(i)
        ); 
    }
    //--------------------------------------------------------------------------------------------------



    //--------------------------------------------------------------------------------------------------
    //This second section get the section of the input mesh file containing the face connectivity information
    //Get the position of the line starting by CELLS
    std::smatch match_3;
    std::regex rgx_3(R"(CELLS ([0-9]+) ([0-9])+)");
    if(!std::regex_search(file_content_, match_3, rgx_3)){
        throw mesh_reader_exception("ERROR mesh_reader: CELLS line not found. Input file probably invalid"); 
    }
    auto cell_pos_start = match_3.position(0) + match_3.length(0);

    //Get the number of cells and the number of integers required to store all the faces
    int nb_ints  = std::stoi(match_3[2].str());


    //Try to determine the end position in the text of the CELL section 
    std::smatch match_4;
    std::regex  rgx_4(R"([A-Z])");
    sub_string = file_content_.substr(cell_pos_start);
    
    if(!std::regex_search(sub_string, match_4, rgx_4)){
        throw mesh_reader_exception("ERROR: input mesh file, could not find the end of the CELL section"); 
    }
    auto cell_pos_end = match_4.position(0);


    //If the position of the CELL section has been identified in the text 
    //Copy all the face data in a string
    std::string cell_text = sub_string.substr(1,cell_pos_end);
    //--------------------------------------------------------------------------------------------------


    //--------------------------------------------------------------------------------------------------  
    //This final section extract the face connectivity information and store it for each face in a vector  
    //Store temporaly all the cell faces connection in this 2D vector
    std::vector<std::vector<unsigned>> cell_face_lst_lst;
    std::istringstream iss(cell_text);

    //Each line corresponds to the faces of an individual cell
    int nb_cell_counter = 0; 
    for (std::string line; std::getline(iss, line); )
    {
        //If the line just contains a return line character or a whitespace
        if (line.size() <= 3) continue;

        //Find the first integer of the string
        std::smatch match_5;
        std::regex  rgx_5(R"(^[0-9]+)");

        if(!std::regex_search(line, match_5, rgx_5)){
            throw mesh_reader_exception("ERROR: CELL section of the input mesh file corrupted"); 
        }

        //The number of integer to store all the face connectivity of the cell 
        int nb_cell_data = std::stoi(match_5[0]); 

        std::vector<unsigned> cell_face_lst;
        cell_face_lst.reserve(nb_cell_data);
        auto line_pos_start = match_5.position(0) + match_5.length(0);
        auto line_substring = line.substr(line_pos_start);

        std::smatch match_6;
        std::regex  rgx_6(R"([0-9]+)");
        search_start = line_substring.cbegin();

        //Store the face information of the cell
        while(std::regex_search(search_start, line_substring.cend(), match_6, rgx_6)){
        
            //Throw an error if the node coordinate could not be converetd
            try{
                cell_face_lst.push_back(std::stoi(match_6[0])); 
            }
            catch(const std::invalid_argument& e){
                throw mesh_reader_exception("ERROR: could not read face data of cell " + std::to_string(nb_cell_counter)); 
            }

            //Restart the number search at the next position
            search_start = match_6.suffix().first;
        }


        if(nb_cell_data != cell_face_lst.size()){
            throw mesh_reader_exception("ERROR: could not read the data of cell " + std::to_string(nb_cell_counter)); 
        }

        //Use this counter to make sure all the cells were loaded
        cell_face_lst_lst.push_back(cell_face_lst);
        nb_cell_counter++;
    }
    //--------------------------------------------------------------------------------------------------    




    return cell_face_lst_lst;

}
//-------------------------------------------------------------------------



//-------------------------------------------------------------------------
//Transform the loaded point data and face connectivity data into a vector of mesh objects
std::vector<mesh>  mesh_reader::get_cell_mesh(
        const std::vector<double>& node_pos, 
        const std::vector<std::vector<unsigned>>& face_connectivity_lst_lst
    ) noexcept(false){

    //The mesh correspoonding to each cell will be stored in this vector
    std::vector<mesh> mesh_lst;
    mesh_lst.reserve(face_connectivity_lst_lst.size());

    //Loop over the cells
    for(const auto& cell_face_lst_1D: face_connectivity_lst_lst){

        //Create a new mesh object to store the cell
        mesh m;

        //Get the number of faces in the cell
        const int nb_faces = cell_face_lst_1D[0];
        m.face_point_ids.reserve(nb_faces);

        //Keep track of the nodes used by the face of the given cell
        std::set<unsigned> cell_node_set;

        //Use an iterator to keep track of the position of the loop in the cell face_lst
        //A lot of jumping will be necessary    
        for(auto it = std::next(cell_face_lst_1D.begin(), 1); it != cell_face_lst_1D.end(); ){

            //The first value of the iterator gives the number of nodes in the face
            int nb_nodes_in_face = *it;

            //Make sure that the data of the face is present
            if(std::distance(cell_face_lst_1D.begin(), it + 1 + nb_nodes_in_face) > cell_face_lst_1D.size()){
                throw mesh_reader_exception("ERROR: input mesh file, face connectivity data corrupted"); 
            }

            //Move the iterator one position forward
            ++it;

            //Create a vector to store the face connectivity data
            std::vector<unsigned> face_node_lst;
            face_node_lst.reserve(nb_nodes_in_face);

            //Copy the connectivity data in the vector created for that 
            std::copy(it, it+nb_nodes_in_face, std::back_inserter(face_node_lst));

            //Keep track of the nodes used by the cell
            cell_node_set.insert(face_node_lst.begin(), face_node_lst.end());

            if(face_node_lst.size() != nb_nodes_in_face){
                throw mesh_reader_exception(
                    "ERROR: nb nodes in face not correct "+ std::to_string(nb_nodes_in_face) + " / " + std::to_string(face_node_lst.size())
                    );  
                }

            //Move the iterator to the data of the next face
            it  = std::next(it, nb_nodes_in_face);

            //Store the face in the cell mesh
            m.face_point_ids.push_back(face_node_lst);
        }


        //At this point all the face information has been stored in the mesh object
        //Now replace the global node ids by local ids.

        //CHeck that the correct number of faces has been loaded
        if(m.face_point_ids.size() != nb_faces){
            throw mesh_reader_exception(
                "ERROR: nb faces in mesh not correct "+ std::to_string(m.face_point_ids.size()) + " / " + std::to_string(nb_faces)
            );  
        }


        std::map<unsigned, unsigned> global_to_local_node_id_map; 

        //Loop over the used nodes and store their local and global node ids
        unsigned local_node_id = 0;
        for(unsigned global_node_id: cell_node_set){
            global_to_local_node_id_map.insert({global_node_id, local_node_id}); 
            local_node_id++;

            //Copy the node position in the mesh.node_pos_lst
            std::copy(node_pos.begin()+global_node_id*3, node_pos.begin()+global_node_id*3+3, std::back_inserter(m.node_pos_lst));
        }

        //Go through the faces stored in the cell and modify the global node ids to the local node ids
        for(auto& f: m.face_point_ids){
            for(size_t i = 0; i< f.size(); i++){
                f[i] = global_to_local_node_id_map[f[i]];
            }
        }

        //Store the cell mesh in the list of meshes
        mesh_lst.push_back(m);
    }

    if(mesh_lst.size() != face_connectivity_lst_lst.size()){
        throw mesh_reader_exception(
            "ERROR: something went wrong while loading the meshes "+ std::to_string(mesh_lst.size()) + " / " + std::to_string(face_connectivity_lst_lst.size())
            );  
        }


    return mesh_lst;
}
//-----------------------------------------------------------------------------------------------



//-----------------------------------------------------------------------------------------------
//Find the cell_type_id array in the VTK file and extract the cell type of each cell
std::vector<short>  mesh_reader::get_cell_types() const noexcept(false){

    //Tries to locate in the file where the array is stored
    std::smatch match_1;
    std::regex rgx_1(R"([C|c]ell_type_id [0-9]+ [0-9]+ [a-zA-Z]+)");
    if(!std::regex_search(file_content_, match_1, rgx_1)){
        throw mesh_reader_exception("ERROR mesh_reader: could not locate the cell_type_id array in the input VTK file"); 
    }

    auto cell_id_start = match_1.position(0) + match_1.length(0);

    //Find the end of the cell_type_id array
    std::smatch match_2;
    std::regex  rgx_2(R"([A-Za-z])");
    std::string sub_string = file_content_.substr(cell_id_start);


    //If the cell_type_id array is the last array of the file
    std::string cell_id_array_str;
    if(!std::regex_search(sub_string, match_2, rgx_2)){
        cell_id_array_str = sub_string;
    }
    else{
        auto coord_pos_end =  match_2.position(0);
        cell_id_array_str = sub_string.substr(0,coord_pos_end);
    }


    //Store all the node positions in this vector
    std::vector<short> cell_id_lst;

    //Very complex regular expression that matches all decimal numbers
    std::regex rgx_3(R"(([0-9]+))");
    std::smatch matches_3;
    std::string::const_iterator search_start(cell_id_array_str.cbegin() );

    //Read the node coordinates and save them
    while(std::regex_search(search_start, cell_id_array_str.cend(), matches_3, rgx_3)){
        
        //Make sure the position value is in a valid format
        short  cell_id_value;

        //Throw an error if the node coordinate could not be converetd
        try{
            cell_id_value = std::stoi(matches_3[0]); 
        }
        catch(const std::invalid_argument& e){
            throw mesh_reader_exception("ERROR: input mesh file, impossible cell type id conversion: " + std::string(matches_3[0])); 
        }

        //Make sure the number is finite and not NA
        if(!std::isfinite(static_cast<double>(cell_id_value))){
            throw mesh_reader_exception("ERROR: input mesh file, the cell type id is not finite: " + std::string(matches_3[0])); 
        }

        //If everything is okay save the node value
        cell_id_lst.push_back(cell_id_value);

        //Restart the number search at the next position
        search_start = matches_3.suffix().first;
    }

    return cell_id_lst;

}