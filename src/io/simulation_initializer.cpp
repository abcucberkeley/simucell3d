#include "simulation_initializer.hpp"

/*
   This class contains all the methods that takes care of the initialization of the simulation.
   It first reads the parameters from the XML file, then it creates the cells and finally it triangulates the cell surfaces.
   All the simulation parameters are stored in structures that can then  be used by the solver class to run the simulation.
*/



//-----------------------------------------------------------------------------------------------------------
//Initialize the simulation by directly providing the structures containing the simulation and cell parameters
simulation_initializer::simulation_initializer(
    const global_simulation_parameters& sim_parameters,
    const std::vector<cell_type_param_ptr>& cell_type_param_lst,
    const bool verbose
)  noexcept(false): verbose_(verbose), sim_parameters_(sim_parameters), perform_initial_triangulation_(sim_parameters.perform_initial_triangulation_){

    //Start the initialization process
    run(cell_type_param_lst);
}
//-----------------------------------------------------------------------------------------------------------



//-----------------------------------------------------------------------------------------------------------
//Initialize the simulation with the path to an XML parameter file
simulation_initializer::simulation_initializer(
    const std::string& parameter_file_path, 
    const bool verbose
    ) noexcept(false) : verbose_(verbose){
    //Read the parameters
    parameter_reader xml_reader(parameter_file_path);

    //Get the parameters that will be passed to the solver
    sim_parameters_ = xml_reader.read_numerical_parameters();

    //If the cells should be triangulated or not at the beginning of the simulation
    perform_initial_triangulation_ = sim_parameters_.perform_initial_triangulation_;

    //Get the list of cell_type parameters
    std::vector<cell_type_param_ptr> cell_type_param_lst = xml_reader.read_biomechanical_parameters();

    //Start the initialization process
    run(cell_type_param_lst);
}
//-----------------------------------------------------------------------------------------------------------




//-----------------------------------------------------------------------------------------------------------
//Starts the initialization process
void simulation_initializer::run(const std::vector<cell_type_param_ptr>& cell_type_param_lst) noexcept(false){

    //Make sure that all the cell types have at least one face type defined
    if(!std::all_of(cell_type_param_lst.begin(), cell_type_param_lst.end(), [](cell_type_param_ptr ctp){return ctp->face_types_.size() > 0;})){
        throw intialization_exception("All cell types must have at least one face type defined");
    }

    //Read the mesh file
    mesh_reader mesh_reader(sim_parameters_.input_mesh_path_, verbose_);
    std::vector<mesh> cell_mesh_lst = mesh_reader.read();
    const size_t nb_cells = cell_mesh_lst.size();

    //Get the type ids of the cells
    std::vector<short> cell_type_id_lst = mesh_reader.get_cell_types();

    //Make sure the number of cells and the number of cell types are the same
    if(nb_cells != cell_type_id_lst.size()){
        throw intialization_exception("The number of cells and the number of cell types are not the same, in the input mesh file");
    }

    //Fill the cell_lst with nullptr
    cell_lst_.resize(nb_cells, nullptr);


    //Wraps the call to the triangulation function in this lambda function
    //This lambda can then be passed to the parrallel_exception_manager that will
    //triangulate the cells in parallel and rethrow properly any exception that might occur
    const std::function<void(size_t)> triangulate_cell_lambda_wrapper  = [&, this](size_t cell_id) -> void {

        //Get the type of the cell
        const short cell_type_id = cell_type_id_lst[cell_id];
        
        //Make sure the cell type id makes sense
        if(cell_type_id >= cell_type_param_lst.size()){
            throw intialization_exception("The cell type id of cell " + std::to_string(cell_id) + " ("+std::to_string(cell_type_id)+") is not valid");
        }
    
        //Get the cell mesh
        mesh& cell_mesh = cell_mesh_lst[cell_id];

        //Tries to triangulate the cell surface
        cell_ptr c = triangulate_surface(cell_mesh, cell_id, cell_type_param_lst[cell_type_id]);
        assert(c != nullptr);

        //Add the cell to the list of cells
        cell_lst_[cell_id] = c;

    
    };

    //Generate the a vector with the ids of the cells
    std::vector<size_t> cell_id_lst(nb_cells);
    std::iota(cell_id_lst.begin(), cell_id_lst.end(), 0);

    //Call the triangulation function in parallel, rethrow the exception if any is thrown
    parallel_exception_handler(cell_id_lst, triangulate_cell_lambda_wrapper);


}
//-----------------------------------------------------------------------------------------------------------



//-----------------------------------------------------------------------------------------------------------
//This function makes sure that the cells properly triangulated and then create an instance of the cells 
//with the correct type. The currently implemented cell types are:
//  - 0: epithelial cell
//  - 1: ECM cell
//  - 2: lumen cell
//  - 3: nucleus cell
//  - 4: static cell
cell_ptr simulation_initializer::triangulate_surface(
        const mesh& cell_mesh, 
        const size_t cell_id,
        cell_type_param_ptr cell_type
    ) noexcept(false){
    cell_ptr c0 = nullptr;


    #pragma omp critical
    {
        if(verbose_) std::cout << "Triangulating cell " << cell_id << " of type: " << cell_type->name_ << std::endl;
    }
    
       
    //Try to triangulate the cell surface 10 times before throwing an exception
    constexpr short max_nb_tries = 10;
    for(short i = 0; i < max_nb_tries; ++i){

        try{

            //Create a mesh structure to store the triangulated cell surface
            mesh cell_surface_mesh;
            if(perform_initial_triangulation_){
                cell_surface_mesh = initial_triangulation::triangulate_surface(sim_parameters_.min_edge_len_, 3. * sim_parameters_.min_edge_len_, cell_mesh, cell_id);
            }
            else{
                
                //Make sure that the input cell mesh is triangulated
                if(!std::all_of(cell_mesh.face_point_ids.begin(), cell_mesh.face_point_ids.end(), [](const auto& f) -> bool {return f.size() == 3;})){
                    throw intialization_exception("If you disable the initial triangulation, the input cell meshes must already be triangulated");
                }

                cell_surface_mesh = cell_mesh;
            }
    
            //Create a cell object of the correct type
            switch(cell_type->global_type_id_){
                case 0: {c0 = std::make_shared<epithelial_cell> (cell_surface_mesh, cell_id, cell_type);    break;}
                case 1: {c0 = std::make_shared<ecm_cell>        (cell_surface_mesh, cell_id, cell_type);    break;}
                case 2: {c0 = std::make_shared<lumen_cell>      (cell_surface_mesh, cell_id, cell_type);    break;}
                case 3: {c0 = std::make_shared<nucleus_cell>    (cell_surface_mesh, cell_id, cell_type);    break;}
                case 4: {c0 = std::make_shared<static_cell>     (cell_surface_mesh, cell_id, cell_type);    break;}

                default:{
                    throw intialization_exception("The cell type id of cell " + std::to_string(cell_id) + 
                    " ("+std::to_string(cell_type->global_type_id_)+") is not valid");
                }
            }
        
            //Initialize the cell properties
            c0->initialize_cell_properties();
            break;
        
        }
        catch(const std::exception& e){
            #pragma omp critical
            {
                 std::cerr << "Triangulation of cell "+ std::to_string(cell_id) + " will be restarted. "
                + std::string("Reason: ") + e.what() << std::endl;
            }   
        }
        if(i == max_nb_tries - 1) throw intialization_exception("The cell " + std::to_string(cell_id) + " could not be triangulated");  
    } 

    return c0;
}

//-----------------------------------------------------------------------------------------------------------








