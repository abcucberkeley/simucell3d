#include "cell.hpp"


/*
    The cell class is a mediator (see mediator pattern design). It manages the relationships between the nodes, the edges
    and the faces and make sure that the integrity of the cell surface is maintained. This limits
    the dependencies between the nodes and faces, and therefore reduces the memory usage.
*/


//---------------------------------------------------------------------------------------------------------
//Instantiate the cell by directly giving it its nodes and faces as well as its parameters
cell::cell(
            const std::vector<node>& node_lst,
            const std::vector<face>& face_lst,
            unsigned cell_id,
            cell_type_param_ptr cell_type_parameters
        ) noexcept :cell_id_(cell_id), cell_type_(cell_type_parameters)
{

    node_lst_ = node_lst;
    face_lst_ = face_lst;
}

//---------------------------------------------------------------------------------------------------------



//---------------------------------------------------------------------------------------------------------
//Construct the cell based on the node position and the node  ids of the faces
//The node_position vector is structured in the following way [x0,y0,z0,x1,y1,z1,...,xn,yn,zn]
//The face_node_ids vector is structured in the following way [f0_n0,f0_n1, f0_n2,f1_n0, f1_n1,f1_n2,...,fn_n0,fn_n1,fn_n2]
//Where f0_n0 is the id of the first node of the first face, f0_n1 is the id of the second node of the first face, etc...
cell::cell(
            const std::vector<double>&  node_position,
            const std::vector<unsigned>& face_node_ids,
            unsigned cell_id,
            cell_type_param_ptr cell_type_parameters
        ) noexcept :cell_id_(cell_id), cell_type_(cell_type_parameters)
{

    assert(node_position.size()  % 3 == 0);
    assert(face_node_ids.size()  % 3 == 0);

    //Make sure the min and max face id are in the range of available nodes
    const auto [min_node_id, max_node_id] = std::minmax_element(face_node_ids.begin(), face_node_ids.end());
    assert(*min_node_id == 0); assert(*max_node_id <=  (node_position.size() / 3)  -1);

    //Get the nb of nodes and faces
    auto nb_nodes = static_cast<size_t>(node_position.size() / 3);
    auto nb_faces = static_cast<size_t>(face_node_ids.size() / 3);

    //Check that the number of nodes or faces do not exceed the limits allowed by the identifier system
    assert(nb_nodes < std::numeric_limits<unsigned>::max());
    assert(nb_faces < std::numeric_limits<unsigned>::max());

    //Starts by reserving enough size in the 2 vectors
    node_lst_.reserve(nb_nodes);
    face_lst_.reserve(nb_faces);

    //Generate the node objects
    for(unsigned node_id = 0; node_id < nb_nodes; node_id++){
        node_lst_.emplace_back(node_position[node_id*3    ],node_position[node_id*3 + 1],node_position[node_id*3 + 2],node_id);
    }


    //Now generate the face objects
    for(unsigned face_id = 0; face_id < nb_faces; face_id++){
        face_lst_.emplace_back(face_node_ids[face_id*3    ],face_node_ids[face_id*3 + 1],face_node_ids[face_id*3 + 2],face_id);
    }
}
//---------------------------------------------------------------------------------------------------------


//---------------------------------------------------------------------------------------------------------
//Convert a mesh object into a cell object
//The mesh object structure is the following:
//mesh.face_point_ids = {{f0_n0,f0_n1,f0_n2}, {f1_n0,f1_n1,f1_n2},...,{fn_n0,fn_n1,fn_n2}]
//mesh.node_pos_lst =   {x0,y0,z0, x1,y1,z1, ..., xn,yn,zn}               
cell::cell(
            const mesh& m,
            unsigned cell_id,
            cell_type_param_ptr cell_type_parameters
        ) noexcept: cell_id_(cell_id), cell_type_(cell_type_parameters){


    //Make sure the number of node coordinates make sense
    assert(m.node_pos_lst.size() % 3 == 0);

    //Make sure all the faces are triangulated
    assert(
        std::all_of(m.face_point_ids.begin(), m.face_point_ids.end(), [](const auto& f) -> bool {return f.size() == 3;})
    );

    //Get the nb of nodes and faces
    const auto nb_nodes = static_cast<size_t>(m.node_pos_lst.size() / 3);
    const auto nb_faces = static_cast<size_t>(m.face_point_ids.size());

    //Check that the number of nodes or faces do not exceed the limits allowed by the identifier system
    assert(nb_nodes < std::numeric_limits<unsigned>::max());
    assert(nb_faces < std::numeric_limits<unsigned>::max());

    //Starts by reserving enough size in the 2 vectors
    node_lst_.reserve(nb_nodes);
    face_lst_.reserve(nb_faces);
 
    //Generate the node objects
    for(unsigned node_id = 0; node_id < nb_nodes; node_id++){
        node_lst_.emplace_back(m.node_pos_lst[node_id*3    ],m.node_pos_lst[node_id*3 + 1],m.node_pos_lst[node_id*3 + 2], node_id);
    }


    //Now generate the face objects
    for(unsigned face_id = 0; face_id < nb_faces; face_id++){
        face_lst_.emplace_back(m.face_point_ids[face_id][0], m.face_point_ids[face_id][1], m.face_point_ids[face_id][2], face_id);
    }

}
//---------------------------------------------------------------------------------------------------------



//---------------------------------------------------------------------------------------------------------
//Initialize all the cell data e.g. cell edges, face areas, cell volume etc
void cell::initialize_cell_properties(bool check_cell_integrity) noexcept(false){ //throw mesh_integrity_exception
    assert(node_lst_.size() != 0);
    assert(face_lst_.size() != 0);

    //Set the ids of faces and nodes of the cell
    set_local_ids();

    //Make sure that all the nodes that are not used by any face are marked as free
    remove_unused_nodes(); 

    //Indicate to all the faces that they belong to this cell
    set_face_owner_cell();

    //Generate the edges of the cell
    if(check_cell_integrity){

        generate_edge_set(); //throw mesh_integrity_exception

        if(!is_manifold()) throw initial_triangulation_exception(std::string("The surface of cell ") + std::to_string(cell_id_) + std::string(" is not properly triangulated."));

        //Orient the face normals outward
        check_face_normal_orientation();
    } 

    //Compute the area of each face
    update_all_face_normals_and_areas();

    //Update the area of the cell
    area_ = compute_area();
    target_area_ = area_;

    //Update the cell volume
    volume_ = compute_volume();
    target_volume_ = volume_;

    //Initialize the cell growth rate and division volumes that are drawn from a normal distribution
    if(cell_type_ != nullptr) initialize_random_properties();

}
//---------------------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------------------
//Initialize the cell growth rate and division volume
void cell::initialize_random_properties() noexcept{
    assert(cell_type_ != nullptr);

    //The growth rate and division volume are drawn from a normal distribution
    //The mean and standard deviation are defined in the cell type parameters

    //Set the growth rate of the cell
    if(cell_type_->std_growth_rate_ != 0.){
        std::minstd_rand gen(std::chrono::system_clock::now().time_since_epoch().count());
        std::normal_distribution<double> distribution(cell_type_->avg_growth_rate_, cell_type_->std_growth_rate_);
        growth_rate_ = distribution(gen);

        //Cap the value between -3x std and +3x std
        if(growth_rate_ > cell_type_->avg_growth_rate_ + 3*cell_type_->std_growth_rate_) growth_rate_ = cell_type_->avg_growth_rate_ + 3*cell_type_->std_growth_rate_;
        if(growth_rate_ < cell_type_->avg_growth_rate_ - 3*cell_type_->std_growth_rate_) growth_rate_ = cell_type_->avg_growth_rate_ - 3*cell_type_->std_growth_rate_;

    }else{
        growth_rate_ = cell_type_->avg_growth_rate_;
    }

    //Set the division volume of the cell
    if(cell_type_->std_division_vol_ != 0. && !std::isinf(cell_type_->avg_division_vol_)){
        std::minstd_rand gen(std::chrono::system_clock::now().time_since_epoch().count());
        std::normal_distribution<double> distribution(cell_type_->avg_division_vol_, cell_type_->std_division_vol_);
        division_volume_ = distribution(gen);

        //Cap the value between -3x std and +3x std
        if(division_volume_ > cell_type_->avg_division_vol_ + 3*cell_type_->std_division_vol_) division_volume_ = cell_type_->avg_division_vol_ + 3*cell_type_->std_division_vol_;
        if(division_volume_ < cell_type_->avg_division_vol_ - 3*cell_type_->std_division_vol_) division_volume_ = cell_type_->avg_division_vol_ - 3*cell_type_->std_division_vol_;

    }else{
        division_volume_ = cell_type_->avg_division_vol_;
    }
}



//---------------------------------------------------------------------------------------------------------
//Remove all the nodes and faces of the cell to avoid cyclic dependies betweeen the faces and the cell
void cell::clear_data() noexcept{
    node_lst_.clear();
    face_lst_.clear();
    edge_set_.clear();
    free_node_queue_.clear();
    free_face_queue_.clear();
    area_ = 0.0;
    volume_ = 0.0;
    target_volume_ = 0.0;
}
//---------------------------------------------------------------------------------------------------------


//---------------------------------------------------------------------------------------------------------
//Make sure that all the nodes that are not used by any face are marked as free
void cell::remove_unused_nodes() noexcept {
    assert(node_lst_.size() !=0);
    assert(face_lst_.size() !=0);

    //Make sure no nodes are marked as free
    assert(free_node_queue_.size() == 0);

    //Create a vector of bool to keep track of the nodes that are used by a face
    std::vector<bool> is_node_used(node_lst_.size(), false);

    //Find all the nodes that are not used by any face
    for(const face& f: face_lst_){ 
        if(f.is_used()){
            const auto f_node_ids = f.get_node_ids();
            for(unsigned node_id: f_node_ids){is_node_used[node_id] = true;}
        }
    }

    //Mark all the nodes that are not used by any face as free
    for(unsigned node_id = 0; node_id < node_lst_.size(); node_id++){
        if(!is_node_used[node_id]){
            node_lst_[node_id].reset();
            free_node_queue_.push_back(node_id);
        }
    }


}
//---------------------------------------------------------------------------------------------------------


#if CONTACT_MODEL_INDEX == 1 || CONTACT_MODEL_INDEX == 2
    //---------------------------------------------------------------------------------------------------------
    //Use the cotan-laplace operator to compute the mean curvature of the surface at each node
    void cell::compute_node_curvature_and_normals() noexcept{

        //Reset all the mean curvatures to 0, and the normals to 0
        std::for_each(node_lst_.begin(), node_lst_.end(), [=](node& n) -> void {
            n.curvature_ = 0;
            n.normal_.reset(0., 0., 0.);
        });


        //First use the faces to get the orientation of the nodes normals
        for(const face& f: face_lst_){
            if(f.is_used()){
            
                //Get the ids of the nodes composing the face
                const auto [n_id1, n_id2, n_id3] = f.get_node_ids();
                assert(n_id1 < node_lst_.size() && n_id2 < node_lst_.size() && n_id3 < node_lst_.size());

                //Add the face normal to the nodes
                node_lst_[n_id1].normal_.translate(f.get_normal() * f.get_area());
                node_lst_[n_id2].normal_.translate(f.get_normal() * f.get_area());
                node_lst_[n_id3].normal_.translate(f.get_normal() * f.get_area());
            }
        }


        //Then use the mean curvature normal used via the laplace beltrami operator as the normal of the surface at the node
        //The normal obtained via the laplace beltami operator is usually way more accurate than averaging the normals of the faces.
        //Howerver, when the mean curvature of the surface is small, the normal obtained via the laplace beltrami operator is not accurate. 
        //Therefore, we use the normal obtained via the laplace beltrami operator only when the mean curvature is large enough.

        //Compute a curvature threshold below which the normal obtained via the laplace beltrami operator is not used
        const double curvature_threshold = 1. / std::cbrt(3. * volume_ / (4. * M_PI));


        //For each node we need to calculate it's incident area (it's the area of the voronoi region around the node
        std::vector<double> incident_area(node_lst_.size(), 0.);

        //Keep track of the mean curvature normal in this vector
        std::vector<vec3> mean_curvature_normal(node_lst_.size(), vec3(0., 0., 0.));
   
        //Loop over the edges
        for(const edge& e: edge_set_){


            /*
                The whole arrangement looks something like that    
                
                n3
                / \     
               /   \      
              /  f1 \   
           n1/_______\n2
             \       /   
              \  f2 /   
               \   /     
                \ /      
                n4

            */


            //Extract the ids of the 2 nodes and the 2 faces stored in the edge object
            const auto [n1_id, n2_id] = e.get_node_ids();
            const auto [f1_id, f2_id] = e.get_face_ids();

            //Get the nodes and faces objects 
            assert(n1_id < node_lst_.size() && n2_id < node_lst_.size());
            assert(f1_id < face_lst_.size() && f2_id < face_lst_.size());

            node& n1 = node_lst_[n1_id];
            node& n2 = node_lst_[n2_id];

            assert(n1.is_used() && n2.is_used());
            const face& f1 = face_lst_[f1_id];
            const face& f2 = face_lst_[f2_id];

            //Add the incident areas of the 2 nodes
            incident_area[n1_id] += f1.get_area() / 3.;
            incident_area[n2_id] += f2.get_area() / 3.;

            //Get the nodes opposed to the edge in the face f1 and f2
            const unsigned n3_id = f1.get_opposite_node(n1_id, n2_id);
            const unsigned n4_id = f2.get_opposite_node(n1_id, n2_id);

            assert(n3_id < node_lst_.size() && n4_id < node_lst_.size());
            const node& n3 = node_lst_[n3_id];
            const node& n4 = node_lst_[n4_id];

            //Calculate the angle between the edges n4--n1 and n4--n2
            const vec3 vec_n4_to_n1 = n1.pos_ - n4.pos_;
            const vec3 vec_n4_to_n2 = n2.pos_ - n4.pos_;
            const double angle_1 = vec_n4_to_n1.get_angle_with(vec_n4_to_n2);

            //Calculate the angle between the edges n3--n1 and n3--n2
            const vec3 vec_n3_to_n1 = n1.pos_ - n3.pos_;
            const vec3 vec_n3_to_n2 = n2.pos_ - n3.pos_;
            const double angle_2 = vec_n3_to_n1.get_angle_with(vec_n3_to_n2);

            //Compute a constat that will be reused 
            const double constant = cot(angle_1) + cot(angle_2);
            if(!std::isfinite(constant)) continue;

            mean_curvature_normal[n1_id] = mean_curvature_normal[n1_id] + (n2.pos_ - n1.pos_) * constant;
            mean_curvature_normal[n2_id] = mean_curvature_normal[n2_id] + (n1.pos_ - n2.pos_) * constant;
        }

        //Now compute the mean curvature vector at the node
        for(size_t node_id = 0; node_id < node_lst_.size(); node_id++){
            node& n=  node_lst_[node_id];
            if(n.is_used()){

                //Compute the curvature of the node
                n.curvature_ = (incident_area[node_id] == 0.) ? 0. : mean_curvature_normal[node_id].norm() / (4. * incident_area[node_id]);
                
                //Compute the normal at the node
                if(n.curvature_> curvature_threshold){
                    n.normal_ = (n.normal_.dot(mean_curvature_normal[node_id]) >=0) ? mean_curvature_normal[node_id] : mean_curvature_normal[node_id] * (-1.);
                }

                //Normalize the normal
                n.normal_ = n.normal_.normalize();
            }
        }
    }





#endif
//---------------------------------------------------------------------------------------------------------




//---------------------------------------------------------------------------------------------------------
//A free node is a node that is not used by any face of the cell mesh. These nodes are not immediately removed
//from the cell meshes because this would cause frequent resizing of the node_lst_ vector.
void cell::add_free_node(const unsigned node_id) noexcept {
    assert(node_id < node_lst_.size());
    assert(node_lst_[node_id].is_used() == false);
    free_node_queue_.push_back(node_id);
}

//Same for the faces, when some faces are not used by the mesh of the cell, they are marked as unused/free
void cell::add_free_face(const unsigned face_id) noexcept {
    assert(face_id < face_lst_.size());
    assert(face_lst_[face_id].is_used() == false);
    free_face_queue_.push_back(face_id);
}
//---------------------------------------------------------------------------------------------------------







//---------------------------------------------------------------------------------------------------------
//Generate the set of edges of the cell from the faces
void cell::generate_edge_set() noexcept(false){

    //Check that the cell is not empty
    assert(node_lst_.size() != 0);
    assert(face_lst_.size() != 0);

    //Loop over the faces and used their node_ids to construct the set of edges of the cell 
    for(const face& f: face_lst_){

        //Get the ids of the nodes composing the face
        const auto [n_id1, n_id2, n_id3] = f.get_node_ids();

        auto [edge_it_1, insert_bool_1] = edge_set_.emplace(n_id1, n_id2);
        auto [edge_it_2, insert_bool_2] = edge_set_.emplace(n_id2, n_id3);
        auto [edge_it_3, insert_bool_3] = edge_set_.emplace(n_id3, n_id1);


        //The cost qualifier hqs to be removed in order to modify the edge objects
        edge& edge_1 = const_cast<edge&>(*edge_it_1);
        edge& edge_2 = const_cast<edge&>(*edge_it_2);
        edge& edge_3 = const_cast<edge&>(*edge_it_3);

        edge_1.add_face(f.get_local_id());
        edge_2.add_face(f.get_local_id());
        edge_3.add_face(f.get_local_id());
    }
}

//---------------------------------------------------------------------------------------------------------
//Returns a copy of the edge if it exists, otherwise returns an empty optional
std::optional<edge> cell::get_edge(const unsigned n1_id, const unsigned n2_id) const noexcept{

    //This is fast since we are using a set complexity is O(log(n))
    auto edge_it = edge_set_.find(edge(n1_id, n2_id));

    if(edge_it == edge_set_.end()) return std::nullopt;

    return *(edge_it);
}
//---------------------------------------------------------------------------------------------------------



//---------------------------------------------------------------------------------------------------------
//Check if the cell has a hole in its surface
bool cell::is_manifold() const noexcept{
    assert(edge_set_.size() > 0);

    //If not all the edges are connected to 2 faces
    if(!(std::all_of(edge_set_.begin(), edge_set_.end(), [](const auto& e) -> bool {return e.is_manifold();}))) return false;

    //The cell mesh should follow the euler formula V - E + F = 2

    //Get the number of nodes, edges and faces
    const int nb_nodes = static_cast<int>(get_nb_of_nodes());
    const int nb_edges = static_cast<int>(edge_set_.size());
    const int nb_faces = static_cast<int>(get_nb_of_faces());
    if(nb_nodes - nb_edges + nb_faces != 2) return false;

    return true;
}
//---------------------------------------------------------------------------------------------------------



//---------------------------------------------------------------------------------------------------------
//Remove a face from the cell mesh, the face is not deleted per see but marked as unused
void cell::delete_face(const unsigned face_id) noexcept(false){
    assert(face_id < face_lst_.size());
    face& f = face_lst_[face_id];
    assert(f.is_used() == true);

    //Remove the face from the edges it is composed of
    auto edge_it_1 = edge_set_.find(edge(f.n1_id_, f.n2_id_));
    auto edge_it_2 = edge_set_.find(edge(f.n2_id_, f.n3_id_));
    auto edge_it_3 = edge_set_.find(edge(f.n3_id_, f.n1_id_));

    //Make sure the edges are present in the set
    assert(edge_it_1 != edge_set_.end());
    assert(edge_it_2 != edge_set_.end());
    assert(edge_it_3 != edge_set_.end());

    
    //It's possible that the face is not contained in the edge, this can happen during the
    //edge merger operation. The 2 faces that are deleted by this operation at this stage 
    //are not stored in the edge objects anymore.
    if(edge_it_1->f1() == f.get_local_id() || edge_it_1->f2() == f.get_local_id()){
        //If the edge is just connected to the face that is about to be deleted
        if(!edge_it_1->is_manifold()){edge_set_.erase(edge_it_1);}
        else{const_cast<edge&>(*edge_it_1).delete_face(f.get_local_id());}
    }
    if(edge_it_2->f1() == f.get_local_id() || edge_it_2->f2() == f.get_local_id()){
        if(!edge_it_2->is_manifold()){edge_set_.erase(edge_it_2);}
        else{const_cast<edge&>(*edge_it_2).delete_face(f.get_local_id());}
    }
    if(edge_it_3->f1() == f.get_local_id() || edge_it_3->f2() == f.get_local_id()){
        if(!edge_it_3->is_manifold()){edge_set_.erase(edge_it_3);}
        else{const_cast<edge&>(*edge_it_3).delete_face(f.get_local_id());}
    }


    const unsigned f_id = f.get_local_id();

    //Mark the face as unused and remove its node data
    f.reset();

    //Add the face to the free face queue
    add_free_face(f_id);

}

void cell::delete_face(face& f) noexcept(false){
    delete_face(f.get_local_id());
}
//---------------------------------------------------------------------------------------------------------



//---------------------------------------------------------------------------------------------------------
//Remove a node from the cell mesh, the not is not deleted per see but marked as unused
//Be very careful when using this function, the faces using the nodes have to be a
//deleted before the node itself
void cell::delete_node(const unsigned node_id) noexcept{
    assert(node_id < node_lst_.size());
    node& n = node_lst_[node_id];
    assert(n.is_used() == true);

    //Mark the node as unused and remove its node data
    n.reset();

    //Add the face to the free face queue
    add_free_node(n.get_local_id());
}

void cell::delete_node(node& n) noexcept{
    delete_node(node_lst_[n.get_local_id()]);
}
//---------------------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------------------
//Add a node to the cell mesh and returns its id
unsigned cell::create_node(const vec3& pos) noexcept{
    node n(pos, 0);
    return add_node(n);
}

unsigned cell::create_node(const double x, const double y, const double z) noexcept{
    node n(x, y, z, 0);
    return add_node(n);
}
//---------------------------------------------------------------------------------------------------------


//---------------------------------------------------------------------------------------------------------
//Add a new node or face to the cell mesh
unsigned cell::add_node(const node& n) noexcept{

    //If there are free slots in the node_lst_ vector, we use them
    if(!free_node_queue_.empty()){

        //Get the id of the first free slot
        const unsigned free_node_id = free_node_queue_.back();

        //Remove the id from the free node queue
        free_node_queue_.pop_back();

        //Add the node to the free slot
        node_lst_[free_node_id] = n;

        //Update the id of the node
        node_lst_[free_node_id].set_local_id(free_node_id);
        node_lst_[free_node_id].set_is_used(true);

        return free_node_id;

    }
    //Else we add the node at the end of the node_lst_ vector
    else{
        node_lst_.push_back(n);

        node_lst_.back().set_local_id(node_lst_.size() - 1);

        node_lst_.back().set_is_used(true);

        return node_lst_.size() - 1;
    }
}
//---------------------------------------------------------------------------------------------------------


//---------------------------------------------------------------------------------------------------------
//Replace the old node with the new node in all its faces and edges, return a vector of the edges that have been deleted
std::pair<std::vector<edge>, std::vector<edge>> cell::replace_node(const edge& start_edge, const unsigned old_node_id, const unsigned new_node_id) noexcept{

    //A buch of tests to make sure the function is called correctly
    assert(old_node_id < node_lst_.size() && new_node_id < node_lst_.size());
    assert(old_node_id != new_node_id);
    assert(start_edge.n1() == old_node_id || start_edge.n2() == old_node_id);
    assert(start_edge.is_manifold());
    assert(node_lst_[old_node_id].is_used());

    //Get the ids of the 2 faces connected to the starting edge
    unsigned face_id =  start_edge.f1();
    unsigned old_face_id;

    //Get an iterator of the edge
    edge_set::const_iterator edge_it = edge_set_.find(edge(start_edge.n1(), start_edge.n2()));

    //Keep track of the edges that have been deleted and created by this operation
    std::vector<edge> deleted_edges;
    std::vector<edge> created_edges;

    //Loop over all the edges and faces connected to the old node
    do{

        assert(edge_it != edge_set_.end());
        assert(edge_it->is_manifold());

        //Get the next face
        old_face_id = face_id;
        face_id = edge_it->f1() == face_id ? edge_it->f2() : edge_it->f1();

        assert(face_id < face_lst_.size());
        assert(face_lst_[face_id].is_used());


        //Replace the old node with the new node in the face
        const vec3 old_normal = face_lst_[face_id].get_normal();

        face_lst_[face_id].replace_node(old_node_id, new_node_id);

        //Update the face normal and area
        update_face_normal_and_area(face_id);

        const double face_area = face_lst_[face_id].get_area();
        const vec3 new_normal = face_lst_[face_id].get_normal();

        //if(new_normal.dot(old_normal) < 0.0) face_lst_[face_id].swap_nodes();

        //We don't use structure bindings here to make sure that the 2 variables below
        //are always in the scope of the whole function and not just the if statement
        edge_set::iterator new_edge_it;
        bool insert_bool;


        //We cannot replace the node in the edge since it would change its hash value
        //Therefore, we have to delete it from the set and create a new one with the updated  node ids
        //Create the new edge, where the old node is replaced by the new node
        std::tie(new_edge_it, insert_bool) = edge_set_.emplace(
            edge_it->n1() == old_node_id ? new_node_id : edge_it->n1(), 
            edge_it->n2() == old_node_id ? new_node_id : edge_it->n2(), 
            edge_it->f1(), edge_it->f2()
        );

        
        //If the new edge already exists, (which happens at the second call of this method by lmr::merge_edge()) 
        //we have to update its faces such that the new edge will not contain one of the 2 faces that is going to be deleted.
        if(!insert_bool){

            const unsigned new_face = (old_face_id == start_edge.f1() || old_face_id == start_edge.f2()) ? face_id : old_face_id;

            edge e(new_edge_it->n1(), new_edge_it->n2(),
                (new_edge_it->f1() == start_edge.f1() || new_edge_it->f1() == start_edge.f2()) ? new_face :  new_edge_it->f1(),
                (new_edge_it->f2() == start_edge.f1() || new_edge_it->f2() == start_edge.f2()) ? new_face :  new_edge_it->f2()
            );
            
            edge_set_.erase(new_edge_it);

            std::tie(new_edge_it, insert_bool)  = edge_set_.insert(e);

            assert(insert_bool);
        }
        
        
        //Delete old edge
        created_edges.push_back(*new_edge_it);
        deleted_edges.push_back(*edge_it);
        edge_set_.erase(edge_it);
       
        //Use this node to get the next edge
        const unsigned opposite_node = face_lst_[face_id].get_opposite_node(new_edge_it->n1(), new_edge_it->n2());

        //Get the next edge
        edge_it = edge_set_.find(edge(old_node_id, opposite_node));

    }
    //This conditions is true when we have looped over all the edges and faces connected to the old node
    while (edge_it != edge_set_.end());

    //Delete the old node
    delete_node(old_node_id);

    return std::make_pair(deleted_edges, created_edges);

}
//---------------------------------------------------------------------------------------------------------



//---------------------------------------------------------------------------------------------------------
//Add a face to the mesh of the cell and returns the id of the newly inserted face
unsigned cell::add_face(const face& f) noexcept(false){

    //Check that the 3 nodes are already in the mesh and that they are not the same
    assert(f.n1_id_ < node_lst_.size() && f.n2_id_ < node_lst_.size() && f.n3_id_ < node_lst_.size());
    assert(f.n1_id_ != f.n2_id_ && f.n2_id_ != f.n3_id_ && f.n1_id_ != f.n3_id_);

    
    //If there are free slots in the face_lst_ vector, we use them
    if(!free_face_queue_.empty()){

        //Get the id of the first free slot
        const unsigned free_face_id = free_face_queue_.back();

        //Remove the id from the free node queue
        free_face_queue_.pop_back();

        //Add the face to the free slot
        face_lst_[free_face_id] = f;
        face& f_ = face_lst_[free_face_id];

        f_.set_local_id(free_face_id);
        f_.set_is_used(true);

        //We add the edges of the face to the edge set
        auto [edge_it_1, insert_bool_1] = edge_set_.emplace(f_.n1_id_, f_.n2_id_);
        auto [edge_it_2, insert_bool_2] = edge_set_.emplace(f_.n2_id_, f_.n3_id_);
        auto [edge_it_3, insert_bool_3] = edge_set_.emplace(f_.n3_id_, f_.n1_id_);

        //Add the face to the 3 edges
        const_cast<edge&>(*edge_it_1).add_face(f_.get_local_id());
        const_cast<edge&>(*edge_it_2).add_face(f_.get_local_id());
        const_cast<edge&>(*edge_it_3).add_face(f_.get_local_id());

        //Update the normal and area of the face
        update_face_normal_and_area(f_);

        f_.owner_cell_ = shared_from_this();

        return f_.get_local_id();
    }
    

    
    //Else we add the face at the end of the face_lst_ vector
    else{
        assert(face_lst_.size() + 1 < std::numeric_limits<unsigned>::max());

        face_lst_.push_back(f);
        face& f_ = face_lst_.back();

        f_.set_local_id(face_lst_.size() - 1);
        f_.set_is_used(true);

        auto [edge_it_1, insert_bool_1] = edge_set_.emplace(f_.n1_id_, f_.n2_id_);
        auto [edge_it_2, insert_bool_2] = edge_set_.emplace(f_.n2_id_, f_.n3_id_);
        auto [edge_it_3, insert_bool_3] = edge_set_.emplace(f_.n3_id_, f_.n1_id_);

        //Add the face to the 3 edges
        const_cast<edge&>(*edge_it_1).add_face(f_.get_local_id());
        const_cast<edge&>(*edge_it_2).add_face(f_.get_local_id());
        const_cast<edge&>(*edge_it_3).add_face(f_.get_local_id());

        //Update the normal and area of the face
        update_face_normal_and_area(f_);

        f_.owner_cell_ = shared_from_this();
        return f_.get_local_id();
    }
}


//---------------------------------------------------------------------------------------------------------



//---------------------------------------------------------------------------------------------------------
//Create a new face from the node ids
unsigned cell::create_face(
            const unsigned n1_id, 
            const unsigned n2_id, 
            const unsigned n3_id  
        ) noexcept(false){

    //The id of the face will anyway be reset when the face is added to the cell
    face f(n1_id, n2_id, n3_id, 0);
    return add_face(f);
}
//---------------------------------------------------------------------------------------------------------


//---------------------------------------------------------------------------------------------------------
//The cell can have some slots of its node_lst and face_lst vectors that are actually
//not used by the mesh. This is done to prevent these 2 vectors from being resized too frequently
//and thereby improve the computational performances. These ununsed slots however must be removed when
//the cell is written in a file or divided. Calling this method is really expensive.

void cell::rebase() noexcept(false){
    //Start by deleting the unsused faces
    //Sort the ids of the faces that can be deleted in reverse order
    const bool regenerate_edge_set = (free_face_queue_.empty() && free_node_queue_.empty()) ? false : true;

    //If some faces are marked as free
    if(free_face_queue_.size() != 0){
        std::sort(free_face_queue_.begin(), free_face_queue_.end());

        //Delete all the faces that are marked as free
        remove_index(face_lst_, free_face_queue_);
        free_face_queue_.clear();

        //Reset the IDs of the remaining faces
        for(unsigned face_id = 0; face_id < face_lst_.size(); face_id++){

            const unsigned old_face_id = face_lst_[face_id].get_local_id();
            face_lst_[face_id].set_local_id(face_id);
        }
    }


    //Do the exact same operation with the nodes
    if(free_node_queue_.size() != 0){
  
        std::sort(free_node_queue_.begin(), free_node_queue_.end());

        remove_index(node_lst_, free_node_queue_);

        free_node_queue_.clear();

        //Create a map to store the correspondence between the old and the new node ids
        std::map<unsigned, unsigned> node_id_correspondence;

        for(unsigned node_id = 0; node_id < node_lst_.size(); node_id++){

            node& n = node_lst_[node_id];

            //                             -------- old node id ------       new node id
            node_id_correspondence.insert({n.get_local_id(), node_id});
            //Relabel the remaining nodes with their correct ids
            n.set_local_id(node_id);

        }

        //Loop over the faces and correct their node ids
        std::for_each(face_lst_.begin(), face_lst_.end(), [&node_id_correspondence](auto& f) -> void {f.update_node_ids(node_id_correspondence);});
    }   


    //Rebuilding the whole edge set is quite expensive and this could be optimized
    if(regenerate_edge_set){
        edge_set_.clear();
        generate_edge_set();
    }


}
//---------------------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------------------
//Return in a flattened 1D vector the coordinates of the nodes
std::vector<double> cell::get_flat_node_coord_lst() const noexcept{
    std::vector<double> node_coord_lst;
    node_coord_lst.reserve(node_lst_.size() * 3);

    for(const node& n: node_lst_){

        //Get the position of the node
        node_coord_lst.push_back(n.pos_.dx());
        node_coord_lst.push_back(n.pos_.dy());
        node_coord_lst.push_back(n.pos_.dz());
    }

    return node_coord_lst;
}
//---------------------------------------------------------------------------------------------------------



//---------------------------------------------------------------------------------------------------------
//A static method to check the winding order of the nodes of a given face is correct compare to an adjacent 
//reference face. If the winding order is not correct the nodes of the face are swapped
void cell::check_face_winding_order(const face& ref_face, face& f) noexcept{

    //Get the local IDs of the 3 nodes constituting the 2 faces
    const auto [ref_n1_id, ref_n2_id, ref_n3_id] = ref_face.get_node_ids();
    const auto [che_n1_id, che_n2_id, che_n3_id] = f.get_node_ids();

    //Get the local ID in the faces of the 2 nodes that they have in common. The local IDs in the face
    //can only be 0, 1 or 2
    std::vector<unsigned> ref_face_local_common_nodes_id;
    std::vector<unsigned> che_face_local_common_nodes_id;

    if(ref_n1_id == che_n1_id){ref_face_local_common_nodes_id.push_back(0); che_face_local_common_nodes_id.push_back(0);}
    if(ref_n1_id == che_n2_id){ref_face_local_common_nodes_id.push_back(0); che_face_local_common_nodes_id.push_back(1);}
    if(ref_n1_id == che_n3_id){ref_face_local_common_nodes_id.push_back(0); che_face_local_common_nodes_id.push_back(2);}

    if(ref_n2_id == che_n1_id){ref_face_local_common_nodes_id.push_back(1); che_face_local_common_nodes_id.push_back(0);}
    if(ref_n2_id == che_n2_id){ref_face_local_common_nodes_id.push_back(1); che_face_local_common_nodes_id.push_back(1);}
    if(ref_n2_id == che_n3_id){ref_face_local_common_nodes_id.push_back(1); che_face_local_common_nodes_id.push_back(2);}

    if(ref_n3_id == che_n1_id){ref_face_local_common_nodes_id.push_back(2); che_face_local_common_nodes_id.push_back(0);}
    if(ref_n3_id == che_n2_id){ref_face_local_common_nodes_id.push_back(2); che_face_local_common_nodes_id.push_back(1);}
    if(ref_n3_id == che_n3_id){ref_face_local_common_nodes_id.push_back(2); che_face_local_common_nodes_id.push_back(2);} 
    assert(ref_face_local_common_nodes_id.size() == 2);

    //Check if these 2 common nodes appear in the same order in the 2 faces
    const bool same_order = ((ref_face_local_common_nodes_id[0] + 1) % 3 == ref_face_local_common_nodes_id[1]) == ((che_face_local_common_nodes_id[0] + 1) % 3 == che_face_local_common_nodes_id[1]);
    
    //If yes, the winding order of the nodes of the face is correct
    if(same_order){std::swap(f.n1_id_, f.n3_id_);}
}
//---------------------------------------------------------------------------------------------------------




//---------------------------------------------------------------------------------------------------------
//Checked that all the faces have their normals oriented outward. 
void cell::check_face_normal_orientation() noexcept{
    assert(face_lst_.size() >= 4); assert(node_lst_.size() >= 4);

    //Use a seed reference face to orient all the facees in the same direction. Then 
    //compute the signed volume of the cell. If the volume is negative, then the faces are oriented inward 
    //and we have to flip the normal of all the faces

    //Store the indices of all the faces that have been already checked
    std::vector<bool> face_checked(face_lst_.size(), false);

    //Store the pairs of faces. The first face is always the reference face, and the second face is the face that is being checked
    typedef std::pair<unsigned, unsigned> face_pair;

    //Store in this vector all the pairs to check
    std::list<face_pair> face_pair_to_check_lst;

    //Loop over the face_lst until a face is_used is found and use it as the reference face
    const auto seed_face_it = std::find_if(face_lst_.begin(), face_lst_.end(), [](const face& f) -> bool {return f.is_used();});
    assert(seed_face_it != face_lst_.end());
    const face& seed_face = *seed_face_it;

    //Set the seed face as checked
    face_checked[seed_face.get_local_id()] = true;

    //Get the local IDs of the 3 nodes constituting the reference face
    const auto [n1_id, n2_id, n3_id] = seed_face.get_node_ids();

    //Get the 3 edges of the face
    const auto edge_1_opt = get_edge(n1_id, n2_id);
    const auto edge_2_opt = get_edge(n2_id, n3_id);
    const auto edge_3_opt = get_edge(n3_id, n1_id);

    assert (edge_1_opt.has_value() && edge_2_opt.has_value() && edge_3_opt.has_value());

    const edge& edge_1 = edge_1_opt.value();
    const edge& edge_2 = edge_2_opt.value();
    const edge& edge_3 = edge_3_opt.value();

    //Get the local IDs of the 3 neighboring faces
    const unsigned f1_id = edge_1.f1() == seed_face.get_local_id() ? edge_1.f2() : edge_1.f1();
    const unsigned f2_id = edge_2.f1() == seed_face.get_local_id() ? edge_2.f2() : edge_2.f1();
    const unsigned f3_id = edge_3.f1() == seed_face.get_local_id() ? edge_3.f2() : edge_3.f1();

    //Make sure that these faces are in the list of faces and that they are used
    assert(f1_id < face_lst_.size() && f2_id < face_lst_.size() && f3_id < face_lst_.size());
    assert(face_lst_[f1_id].is_used() && face_lst_[f2_id].is_used() && face_lst_[f3_id].is_used());

    //Add to the set of pairs of faces to check
    face_pair_to_check_lst.push_back({seed_face.get_local_id(), f1_id});
    face_pair_to_check_lst.push_back({seed_face.get_local_id(), f2_id});
    face_pair_to_check_lst.push_back({seed_face.get_local_id(), f3_id});


    while(!face_pair_to_check_lst.empty()){

        //Get the 2 faces to check
        const auto [ref_face_id, face_to_check_id] = face_pair_to_check_lst.front();

        //Remove the pair from the list
        face_pair_to_check_lst.pop_front();

        //Check if the face_to_check has not been already checked
        if(face_checked[face_to_check_id]) continue;
        else face_checked[face_to_check_id] = true;

        //Get the two faces
        assert(ref_face_id < face_lst_.size() && face_to_check_id < face_lst_.size());
        const face& ref_face = face_lst_[ref_face_id];
        face& face_to_check = face_lst_[face_to_check_id];
        assert (ref_face.is_used() && face_to_check.is_used());

        //Check the winding order of the nodes of the face_to_check
        check_face_winding_order(ref_face, face_to_check);

        //We now need to add the 3 neighboring faces of the face_to_check to the list of pairs to check
        const auto [che_n1_id, che_n2_id, che_n3_id] = face_to_check.get_node_ids();
        const auto edge_1_opt = get_edge(che_n1_id, che_n2_id);
        const auto edge_2_opt = get_edge(che_n2_id, che_n3_id);
        const auto edge_3_opt = get_edge(che_n3_id, che_n1_id);

        assert (edge_1_opt.has_value() && edge_2_opt.has_value() && edge_3_opt.has_value());
        const edge& edge_1 = edge_1_opt.value();
        const edge& edge_2 = edge_2_opt.value();
        const edge& edge_3 = edge_3_opt.value();

        const unsigned new_f1_id = edge_1.f1() == face_to_check_id ? edge_1.f2() : edge_1.f1();
        const unsigned new_f2_id = edge_2.f1() == face_to_check_id ? edge_2.f2() : edge_2.f1();
        const unsigned new_f3_id = edge_3.f1() == face_to_check_id ? edge_3.f2() : edge_3.f1();

        //Check that these faces are in the list of faces
        assert(new_f1_id < face_lst_.size() && new_f2_id < face_lst_.size() && new_f3_id < face_lst_.size());

        //If these neighboring faces have not been checked, we add them to the list of pairs to check
        if(!face_checked[new_f1_id]) face_pair_to_check_lst.push_back({face_to_check_id, new_f1_id});
        if(!face_checked[new_f2_id]) face_pair_to_check_lst.push_back({face_to_check_id, new_f2_id});
        if(!face_checked[new_f3_id]) face_pair_to_check_lst.push_back({face_to_check_id, new_f3_id});
    }

    //At this stage all the faces have been checked and their normals are oriented consistently
    //We can now compute the signed volume of the cell
    double signed_volume = 0.0;
    for(const auto& f: face_lst_){
        if(f.is_used()){
            const auto [n1_id, n2_id, n3_id] = f.get_node_ids();
            const node& n1 = node_lst_[n1_id];
            const node& n2 = node_lst_[n2_id];
            const node& n3 = node_lst_[n3_id];

            signed_volume += (n1.pos_.dx() * n2.pos_.dy() * n3.pos_.dz() - n1.pos_.dx() * n3.pos_.dy() * n2.pos_.dz() - n2.pos_.dx() * n1.pos_.dy() * n3.pos_.dz() + n2.pos_.dx() * n3.pos_.dy() * n1.pos_.dz() + n3.pos_.dx() * n1.pos_.dy() * n2.pos_.dz() - n3.pos_.dx() * n2.pos_.dy() * n1.pos_.dz());
        }
    }

    //If the volume is negative, we have to flip the normal of all the faces
    if(signed_volume < 0.0){
        for(auto& f: face_lst_){
            if(f.is_used()) f.swap_nodes();
        }
    }

}
//---------------------------------------------------------------------------------------------------------




//---------------------------------------------------------------------------------------------------------
//Compute the normal of the given face. WARNING: the normal is not a unit vector!!
vec3 cell::get_face_normal(const unsigned face_id) const noexcept{

    //Make sure the face is in the cell
    assert(face_id < face_lst_.size());
    return face_lst_[face_id].get_normal();
}



//---------------------------------------------------------------------------------------------------------
//Update the area of each face
void cell::update_all_face_normals_and_areas() noexcept{

    std::for_each(face_lst_.begin(), face_lst_.end(), 
        [this](face& f) -> void {if(f.is_used()) update_face_normal_and_area(f);}
    );

}
//---------------------------------------------------------------------------------------------------------




//---------------------------------------------------------------------------------------------------------
//Update at the ssame time the normal and the area of the given face
void cell::update_face_normal_and_area(const unsigned face_id) noexcept{
    assert(face_id < face_lst_.size());
    update_face_normal_and_area(face_lst_[face_id]);
}


void cell::update_face_normal_and_area(face& f) noexcept{

    //Make sure the face is used by nodes
    assert(f.is_used());

    //Get the ids of the face's nodes, and make sure they are part of the cell 
    const unsigned n1_id = f.n1_id_, n2_id = f.n2_id_, n3_id = f.n3_id_; 

    //Make sure that the nodes are part of the cell
    assert(n1_id < node_lst_.size() && n2_id < node_lst_.size() && n3_id < node_lst_.size());

    const node& n1 = node_lst_[n1_id], n2 = node_lst_[n2_id], n3 = node_lst_[n3_id]; 

    //Use the edges of the face to compute the normal
    vec3 face_normal =  (n2 - n1).cross(n3 - n1);
    const double face_normal_norm = face_normal.norm();
    assert(std::isfinite(face_normal_norm));
    assert(face_normal_norm >=0.);

    //Update the face area
    f.set_area(0.5 * face_normal_norm);

    face_normal = face_normal_norm == 0.0 ? vec3(0.0, 0.0, 0.0) : face_normal / face_normal_norm;

    //Update the face normal
    f.set_normal(face_normal);

}
//---------------------------------------------------------------------------------------------------------




//-------------------------------------------------------------------------------------------
//Return the volume of the cell, if signed is set to true return the signed volume
double cell::compute_volume() const noexcept{

    assert(face_lst_.size() >= 4);
    assert(node_lst_.size() >= 4);

    double vol = 0.0;
    for(const auto f : face_lst_){
        //Make sure that all the faces are used by nodes at this stage
        if(f.is_used()){

            //Get the nodes of the face
            const auto [n1_id, n2_id, n3_id] = f.get_node_ids();
            assert(n1_id < node_lst_.size() && n2_id < node_lst_.size() && n3_id < node_lst_.size());
            const node& n1 = node_lst_[n1_id];
            const node& n2 = node_lst_[n2_id];
            const node& n3 = node_lst_[n3_id];

            //Get the coordinates of the nodes
            const auto [x1, y1, z1] = n1.pos_.to_array();
            const auto [x2, y2, z2] = n2.pos_.to_array(); 
            const auto [x3, y3, z3] = n3.pos_.to_array(); 

            //Add the contribution of the face to the signed volume
            vol += -x3*y2*z1 + x2*y3*z1 + x3*y1*z2 - x1*y3*z2 - x2*y1*z3 + x1*y2*z3;

        }                       
    }

    vol /= 6.;
    vol = std::abs(vol);


    //Make sure that the function doesn't return non sense, which can happen very easily
    //if one of the node coordinates equals nan
    assert(std::isfinite(vol));

    //If the volume is not signed, make sure that the returned volume is not equal to 0 
    assert(vol != 0.);

    return vol;
}
//---------------------------------------------------------------------------------------------------------




//---------------------------------------------------------------------------------------------------------
//The name is explicit enough
vec3 cell::compute_centroid() const noexcept{
    
    assert(area_ != 0.);

    vec3 cell_centroid(0., 0., 0.);

    //Use the centroid of each face to compute the centroid of the cell
    for(const auto& f : face_lst_){
        if(f.is_used()){

            //Compute the centroid of the face
            auto [n1_id, n2_id, n3_id] = f.get_node_ids();
            assert(n1_id < node_lst_.size() && n2_id < node_lst_.size() && n3_id < node_lst_.size());

            const node& n1 = node_lst_[n1_id];
            const node& n2 = node_lst_[n2_id];
            const node& n3 = node_lst_[n3_id];
            const vec3 face_centroid  = (n1.pos_ + n2.pos_ + n3.pos_) / 3.;

            //Add the contribution of the face to the cell centroid, weighted by the face area
            cell_centroid.translate(face_centroid * f.get_area());
        }
    }

    //Divide by the total area of the cell
    cell_centroid = cell_centroid / area_;

    return cell_centroid;
}
//---------------------------------------------------------------------------------------------------------



//---------------------------------------------------------------------------------------------------------
//Return the area of the cell const noexcept;
double cell::compute_area() const noexcept{


    const double cell_area = std::accumulate(face_lst_.begin(), face_lst_.end(), 0.,

        [](double sum_area, const face& f) -> double {return sum_area + ((f.is_used()) ? f.get_area() : 0.);}
    
    );
    return cell_area;
}
//---------------------------------------------------------------------------------------------------------




//---------------------------------------------------------------------------------------------------------
const node& cell::get_const_ref_node(const unsigned node_id) const noexcept{
    assert(node_id < node_lst_.size());
    assert(node_lst_[node_id].is_used());
    return node_lst_[node_id];
}
//---------------------------------------------------------------------------------------------------------


//---------------------------------------------------------------------------------------------------------
//Given its ID this method returns a const copy of the face 
const face& cell::get_const_ref_face(const unsigned face_id) const noexcept{
    assert(face_id < face_lst_.size());
    assert(face_lst_[face_id].is_used());
    return face_lst_[face_id];
}
//---------------------------------------------------------------------------------------------------------



//---------------------------------------------------------------------------------------------------------
//Get the axis aligned bounding box of the cell
std::array<double, 6> cell::get_aabb() const noexcept{

    double min_x =   std::numeric_limits<double>::infinity();
    double min_y =   std::numeric_limits<double>::infinity();
    double min_z =   std::numeric_limits<double>::infinity();

    double max_x = - std::numeric_limits<double>::infinity();
    double max_y = - std::numeric_limits<double>::infinity();
    double max_z = - std::numeric_limits<double>::infinity();

    //Loop over the nodes
    for(const node& n: node_lst_){

        if(n.is_used()){

            //Get the position of the node
            const auto [n_x, n_y, n_z] = n.pos_.to_array();

            if(n_x < min_x) min_x = n_x;
            if(n_y < min_y) min_y = n_y;
            if(n_z < min_z) min_z = n_z;

            if(n_x > max_x) max_x = n_x;
            if(n_y > max_y) max_y = n_y;
            if(n_z > max_z) max_z = n_z;
        }
    }

    return {min_x, min_y, min_z, max_x, max_y, max_z};
}
//---------------------------------------------------------------------------------------------------------


//---------------------------------------------------------------------------------------------------------
void cell::update_target_volume(const double time_step) noexcept{
    assert(cell_type_ != nullptr);
    target_volume_ += time_step * growth_rate_;

    //The target volume cannot be smaller than the minimum allowed volume
    if(target_volume_ < cell_type_->min_vol_) target_volume_ = cell_type_->min_vol_;
}
//---------------------------------------------------------------------------------------------------------


//---------------------------------------------------------------------------------------------------------
//Compute the pressure difference between the cell and the environment
void cell::update_pressure() noexcept{
    assert(cell_type_ != nullptr);

    pressure_ = - cell_type_->bulk_modulus_ * std::log(volume_ / target_volume_);
    assert(std::isfinite(pressure_)); //throw unstable_simulation_exception("The pressure is not finite. Cell volume = " + format_number(volume_, "%.2e") + ", target volume = " + format_number(target_volume_, "%.2e"));

    //The pressure can be capped to a maximum value
    if(pressure_ > cell_type_->max_pressure_) pressure_ = cell_type_->max_pressure_;

    //Update the pressure energy of the cell
    pressure_energy_ = std::abs(cell_type_->bulk_modulus_ * volume_ * (std::log(volume_ / target_volume_) - 1.));
}
//---------------------------------------------------------------------------------------------------------



//---------------------------------------------------------------------------------------------------------
//Apply the pressure on the nodes of the cell mesh
void cell::apply_pressure_on_surface() noexcept{

    //Loop over the faces
    for(face& f : face_lst_){

        //Make sure that the face is used
        if(f.is_used()){

            //Get the nodes of the face
            assert(f.n1_id_ < node_lst_.size() && f.n2_id_ < node_lst_.size() && f.n3_id_ < node_lst_.size());

            node& n1 = node_lst_[f.n1_id_];
            node& n2 = node_lst_[f.n2_id_];
            node& n3 = node_lst_[f.n3_id_];

            //Make sure the nodes are not marked as free
            assert(n1.is_used() && n2.is_used() && n3.is_used());

            //Get the normal of the face and its area
            const vec3 face_normal = f.get_normal();
            const double face_area = f.get_area();

            //The pressure that will be applied on the nodes of the face, We divide by 3 because the pressure is distributed on 
            //the 3 nodes of the face
            const vec3 face_pressure = face_normal * pressure_ * face_area / 3.;

            //Apply the pressure on the nodes of the face
            n1.add_force(face_pressure);
            n2.add_force(face_pressure);
            n3.add_force(face_pressure);
        }
    }
}
//---------------------------------------------------------------------------------------------------------





//---------------------------------------------------------------------------------------------------------
//Apply at the same time the surface tension forces and the membrane elasticity forces
void cell::apply_surface_tension_and_membrane_elasticity() noexcept{
    assert(cell_type_ != nullptr);

    //Compute the taget area of the cell based on its isoperimetric ratio
    target_area_ = std::cbrt(cell_type_->target_isoperimetric_ratio_ * volume_ * volume_);
    const double membrane_elasticity_factor =  -(cell_type_->area_elasticity_modulus_ / target_area_) * ((area_ / target_area_) - 1.);

    //Reset the membrane elasticity and surface tension energy  
    membrane_elasticity_energy_ = 0.;
    surface_tension_energy_ = 0.;

    //Loop over the faces
    for(face& f : face_lst_){

        //Make sure that the face is used
        if(f.is_used()){

            if(f.get_area() == 0.) continue;

            //Get the nodes of the face
            assert(f.n1_id_ < node_lst_.size() && f.n2_id_ < node_lst_.size() && f.n3_id_ < node_lst_.size());

            //Get the type of the face
            assert(f.type_id_ < cell_type_->face_types_.size());
            const face_type_parameters& face_type = cell_type_->face_types_[f.type_id_];

            node& n1 = node_lst_[f.n1_id_];
            node& n2 = node_lst_[f.n2_id_];
            node& n3 = node_lst_[f.n3_id_];

            //Compute the 3 area gradients of the face wrt its 3 nodes positions
            const vec3 grad_pos_n1 = f.get_normal().cross(n2.pos_ - n3.pos_) * (-0.5);
            const vec3 grad_pos_n2 = f.get_normal().cross(n3.pos_ - n1.pos_) * (-0.5);
            const vec3 grad_pos_n3 = f.get_normal().cross(n1.pos_ - n2.pos_) * (-0.5);

            const double membrane_elasticity_factor =  -(cell_type_->area_elasticity_modulus_ / target_area_) * ((area_/ target_area_) - 1.);

            
            //Compute the membrane elasticity force factor
            const double force_factor = -face_type.surface_tension_ + membrane_elasticity_factor;
            assert(std::isfinite(force_factor));

            //Compute the forces applied on the 3 nodes of the face
            const vec3 force_n1 = grad_pos_n1 * force_factor;
            const vec3 force_n2 = grad_pos_n2 * force_factor;
            const vec3 force_n3 = grad_pos_n3 * force_factor;

            //Apply the forces on the nodes of the face
            n1.add_force(force_n1);
            n2.add_force(force_n2);
            n3.add_force(force_n3);

            //Update the membrane elasticity energy and surface tension energy
            //membrane_elasticity_energy_ += 0.5 * cell_type_->area_elasticity_modulus_ * std::pow((f.get_area() / target_area_face) - 1., 2.);
            surface_tension_energy_     += 0.5 * face_type.surface_tension_ * f.get_area();
        }
    }
}
//---------------------------------------------------------------------------------------------------------



void cell::apply_bending_forces() noexcept{
    /*For a comprehensive inderstanding of the bending model used here see:
    M. Wardetzky, M. Bergou, D. Harmon, D. Zorin, and E. Grinspun, Computer Aided Geometric Design 24, 499 (2007).
    */
   assert(cell_type_ != nullptr);

   //Angle threshold below which we don't apply any bending force on the nodes of the hinge, to prevent numerical instabilities
   constexpr double max_angle_threshold = 135. * M_PI / 180.;

   //Only apply if some of the face types have a bending rigidicty greater than 0
   if(std::all_of(cell_type_->face_types_.begin(), cell_type_->face_types_.end(), [](const face_type_parameters& f){return f.bending_modulus_ ==  0.;})) return;

    //Reset the bending energy to 0
    bending_energy_ = 0.;

   //Loop over the edges
   for(const edge& e: edge_set_){

        //Make sure the edge is connected to 2 faces
        assert(e.is_manifold());

        //Get the nodes of the edge
        assert(e.n1() < node_lst_.size() && e.n2() < node_lst_.size());
        node& n1 = node_lst_[e.n1()];
        node& n2 = node_lst_[e.n2()];

        //Get the faces of the edge
        assert(e.f1() < face_lst_.size() && e.f2() < face_lst_.size());
        face& f1 = face_lst_[e.f1()];
        face& f2 = face_lst_[e.f2()];

        assert(n1.is_used() && n2.is_used());
        assert(f1.is_used() && f2.is_used());

        //Get the type of the faces
        assert(f1.type_id_ < cell_type_->face_types_.size() && f2.type_id_ < cell_type_->face_types_.size());
        const face_type_parameters& face_type_1 = cell_type_->face_types_[f1.type_id_];
        const face_type_parameters& face_type_2 = cell_type_->face_types_[f2.type_id_];

        //Get the average bending elasticity of the 2 faces
        const double avg_bending_stiffness = (face_type_1.bending_modulus_ + face_type_2.bending_modulus_) / 2.;

        //Get the sum of the 2 faces areas
        const double sum_face_areas = f1.get_area() + f2.get_area();

        //Get the normals of the 2 faces and compute the angle of the hinge
        const vec3& normal_1 = f1.get_normal();
        const vec3& normal_2 = f2.get_normal();

        //Get the nodes opposite to the edge in each face
        const unsigned n3_id = f1.get_opposite_node(e.n1(), e.n2());
        const unsigned n4_id = f2.get_opposite_node(e.n1(), e.n2());
        
        assert(n3_id < node_lst_.size() && n4_id < node_lst_.size());
        node& n3 = node_lst_[n3_id];
        node& n4 = node_lst_[n4_id];

        //Compute the different edge vectors
        const vec3 e0 = n2 - n1;
        const vec3 e1 = n3 - n1;
        const vec3 e2 = n4 - n1;
        const vec3 e3 = n3 - n2;
        const vec3 e4 = n4 - n2;


        //Compute the different angles between the vectors of the diamond region
        const double alpha_1 = e0.get_angle_with(e1);
        const double alpha_2 = e0.get_angle_with(e2);
        const double alpha_3 = e3.get_angle_with(e0 * (-1.));
        const double alpha_4 = e4.get_angle_with(e0 * (-1.));

        //Compute the dihedrral angle between the 2 faces
        double dot = normal_1.dot(normal_2);
        double theta;
        if(dot >= 1.0)          theta = 0.0;
        else if(dot <= -1.0)    theta = M_PI;
        else                    theta = std::acos(dot);

        //If the angle is below a certain thresold, we don't apply any bending force
        if(theta > max_angle_threshold) continue;

        //If the edge is concave
        if (e2.dot(normal_1) > 0) theta = (2* M_PI) - theta;
        theta = M_PI - theta;

        //Compute various factors
        const double edge_length = e0.norm();
        const double prefactor_1 = (-3.*(1. + std::cos(theta)))*avg_bending_stiffness;
        const double prefactor_2 = (3. * edge_length * edge_length / sum_face_areas) * std::sin(theta) * avg_bending_stiffness;
        
        if(edge_length == 0. || f1.get_area() == 0 || f2.get_area() == 0.) continue;
        if(alpha_1 == 0. || alpha_2 == 0. || alpha_3 == 0. || alpha_4 == 0.) continue;

        //Compute the gradient of the angle theta wrt the position of each node position
        const vec3 grad_x0_theta = (normal_1 * cot(alpha_3)   + normal_1 * cot(alpha_4)) * (-1. / edge_length);
        const vec3 grad_x1_theta = (normal_1 * cot(alpha_1)   + normal_1 * cot(alpha_2)) * (-1. / edge_length);
        const vec3 grad_x2_theta = normal_1 * (edge_length / (2. * f1.get_area()));
        const vec3 grad_x3_theta = normal_2 * (edge_length / (2. * f2.get_area()));

        //Compute the tangent of the edge vectors
        const vec3 t1   = e1.rotate_around_axis(normal_1,  M_PI / 2.);
        const vec3 t2   = e2.rotate_around_axis(normal_2, -M_PI / 2.);
        const vec3 t3   = e3.rotate_around_axis(normal_1, -M_PI / 2.);
        const vec3 t4   = e4.rotate_around_axis(normal_2,  M_PI / 2.);
        const vec3 t0_0 = e0.rotate_around_axis(normal_1, -M_PI / 2.);
        const vec3 t0_1 = e0.rotate_around_axis(normal_2,  M_PI / 2.);

        //Compute the gradients of the in plane bending forces
        const double prefactor_3 = edge_length * edge_length / (2. * sum_face_areas * sum_face_areas);
        const vec3 grad_x0_ip = e0 * (-2. / sum_face_areas) + (t3 + t4) * prefactor_3;
        const vec3 grad_x1_ip = e0 * ( 2. / sum_face_areas) + (t1 + t2) * prefactor_3;
        const vec3 grad_x2_ip = t0_0 * prefactor_3;
        const vec3 grad_x3_ip = t0_1 * prefactor_3;

        //Assemble everything together to get the bending forces applied on the nodes of the hinge
        const vec3 f_bend_x0 = grad_x0_ip * prefactor_1 +  grad_x0_theta * prefactor_2;
        const vec3 f_bend_x1 = grad_x1_ip * prefactor_1 +  grad_x1_theta * prefactor_2;
        const vec3 f_bend_x2 = grad_x2_ip * prefactor_1 +  grad_x2_theta * prefactor_2;
        const vec3 f_bend_x3 = grad_x3_ip * prefactor_1 +  grad_x3_theta * prefactor_2;

        //Apply these bending forces to the edges
        n1.add_force(f_bend_x0);
        n2.add_force(f_bend_x1);
        n3.add_force(f_bend_x2);
        n4.add_force(f_bend_x3);

        //Compute the bending energy associated with this edge
        bending_energy_ += (avg_bending_stiffness * edge_length * edge_length / sum_face_areas) * std::pow(2. * std::cos(theta / 2.), 2) ;
   }
}
//---------------------------------------------------------------------------------------------------------



//---------------------------------------------------------------------------------------------------------
void cell::apply_internal_forces(const double time_step) noexcept{
    assert(cell_type_ != nullptr);
    assert(node_lst_.size() >= 4);
    assert(face_lst_.size() >= 4);

    //Update the normals and the areas of all the faces
    update_all_face_normals_and_areas();

    //Update the cell area
    area_ = compute_area();

    //Update the cell volume
    volume_ = compute_volume(); 

    //Update the target volume of the cell
    update_target_volume(time_step);

    //Update the pressure of the cell
    update_pressure();

    //Apply the pressure forces on the nodes of the surface
    apply_pressure_on_surface();

    //Apply the surface tension and membrane elasticity forces
    apply_surface_tension_and_membrane_elasticity();

    //Apply the bending forces
    apply_bending_forces();

    //Regularize the angles of each face to be approximately 60 degrees
    regularize_all_face_angles();

    //If the surfaces of adjacent epithelial cells are mechanically coupled
    #if CONTACT_MODEL_INDEX == 1  || CONTACT_MODEL_INDEX == 2
        //Calculate the mean curvature at the node as well as the surface normal
        compute_node_curvature_and_normals();


    #endif

}
//---------------------------------------------------------------------------------------------------------




//---------------------------------------------------------------------------------------------------------
//Returns the ids of all the nodes that are connected to the given node by an edge
std::vector<unsigned> cell::get_connected_nodes(const unsigned node_id, const edge& start_edge) const noexcept{
    assert(node_id < node_lst_.size());
    assert(start_edge.n1() == node_id || start_edge.n2() == node_id);
    assert(node_lst_[node_id].is_used());


    const unsigned stop_face = start_edge.f1();
    
    //Store the connected nodes in this vector
    std::vector<unsigned> node_ids{(node_id == start_edge.n1() ? start_edge.n2() : start_edge.n1())};


    edge e = start_edge;
    unsigned next_face = start_edge.f2();

    while(next_face != stop_face){


        //Get the next face
        const face& f = get_const_ref_face(next_face);

        //Get the opposite node in this face to the current edge
        const unsigned opp_node = f.get_opposite_node(e.n1(), e.n2());

        //Get the next edge
        auto e_opt = get_edge(opp_node, node_id);
        assert(e_opt.has_value());
        e = e_opt.value();

        node_ids.push_back(node_id == e.n1() ? e.n2() : e.n1());
        

        //Get the next face
        next_face = next_face == e.f1() ? e.f2() : e.f1();
    }

    return node_ids;
}
//---------------------------------------------------------------------------------------------------------


//---------------------------------------------------------------------------------------------------------
//Returns a mesh object of the cell
mesh cell::get_mesh() const noexcept{
    mesh m; 
    m.node_pos_lst = get_flat_node_coord_lst();

    m.face_point_ids = get_face_point_ids();

    return m;
}
//---------------------------------------------------------------------------------------------------------


//---------------------------------------------------------------------------------------------------------
//Return the eigen vector associated with the largest eigen value of the mass covariance matrix
vec3 cell::get_cell_longest_axis() const noexcept{

    //First get the centroid of the cell
    const vec3 centroid = compute_centroid();

    //Compute the number of points that are used in the cell
    const double nb_points = static_cast<double>(get_nb_of_nodes());

    //Compute the covariance matrix of the cell point cloud
    double cov_xx = 0.;
    double cov_xy = 0.;
    double cov_xz = 0.;
    double cov_yy = 0.;
    double cov_yz = 0.;
    double cov_zz = 0.;

    for(const node& n: node_lst_){
        if(n.is_used()){
            cov_xx += (n.pos().dx() - centroid.dx()) * (n.pos().dx() - centroid.dx());
            cov_xy += (n.pos().dx() - centroid.dx()) * (n.pos().dy() - centroid.dy());
            cov_xz += (n.pos().dx() - centroid.dx()) * (n.pos().dz() - centroid.dz());
            cov_yy += (n.pos().dy() - centroid.dy()) * (n.pos().dy() - centroid.dy());
            cov_yz += (n.pos().dy() - centroid.dy()) * (n.pos().dz() - centroid.dz());
            cov_zz += (n.pos().dz() - centroid.dz()) * (n.pos().dz() - centroid.dz());
        }
    }

    cov_xx /= nb_points;
    cov_xy /= nb_points;
    cov_xz /= nb_points;
    cov_yy /= nb_points;
    cov_yz /= nb_points;
    cov_zz /= nb_points;


    //Create the covariance matrix
    const mat33 covariance_matrix(
        {cov_xx, cov_xy, cov_xz},
        {cov_xy, cov_yy, cov_yz},
        {cov_xz, cov_yz, cov_zz}
    );
    
    //Get the eigen values and eigen vectors of the covariance matrix
    const auto [eigen_values, eigen_vectors] = covariance_matrix.eigen_decomposition();

    vec3 longest_axis;
    if(     std::abs(eigen_values.dx()) > std::abs(eigen_values.dy()) && std::abs(eigen_values.dx()) > std::abs(eigen_values.dz())){longest_axis = eigen_vectors.get_col(0); }
    else if(std::abs(eigen_values.dy()) > std::abs(eigen_values.dx()) && std::abs(eigen_values.dy()) > std::abs(eigen_values.dz())){longest_axis = eigen_vectors.get_col(1);}
    else{longest_axis = eigen_vectors.get_col(2);}
    return longest_axis.normalize();
}
//---------------------------------------------------------------------------------------------------------



//---------------------------------------------------------------------------------------------------------
//Given 3 nodes, that form a face with edge vectors a = j - i, b = k - i, 
//return the gradient of the angle between a and b with respect to the coordinates of the 3 nodes
std::array<vec3, 3> cell::get_angle_gradient(const vec3& i, const vec3& j, const vec3& k) noexcept{

    //Construct the edge vectors a and b
    const vec3 a = j - i;
    const vec3 b = k - i;

    //Precompute some dot products
    const double d_ab = a.dot(b);
    const double d_aa = a.dot(a);
    const double d_bb = b.dot(b);

    //Precompute the norm of the vectors a and b
    const double norm_a = std::sqrt(d_aa);  
    const double norm_b = std::sqrt(d_bb);    

    //Precompute some factors that are used repeatedly
    const double denominator = std::sqrt(1 - d_ab * d_ab / (d_aa*d_bb));
    const double factor_1 = norm_b*std::pow(d_aa,1.5);
    const double factor_2 = norm_a*std::pow(d_bb,1.5);


    if(almost_equal(denominator, 0.) || !std::isfinite(denominator) || almost_equal(norm_a, 0.) || almost_equal(norm_b, 0.)) return {vec3(0., 0., 0.), vec3(0., 0., 0.), vec3(0., 0., 0.)};
    if(almost_equal(factor_1, 0.) || !std::isfinite(factor_1)) return {vec3(0., 0., 0.), vec3(0., 0., 0.), vec3(0., 0., 0.)};
    if(almost_equal(factor_2, 0.) || !std::isfinite(factor_2)) return {vec3(0., 0., 0.), vec3(0., 0., 0.), vec3(0., 0., 0.)};


    //Compute the gradients
    const vec3 grad_theta_i(
        -((2*i.dx() - j.dx() - k.dx())/(norm_a*norm_b) + (d_ab*a.dx())/ factor_1 + (d_ab*b.dx())/ factor_2) / denominator,
        -((2*i.dy() - j.dy() - k.dy())/(norm_a*norm_b) + (d_ab*a.dy())/ factor_1 + (d_ab*b.dy())/ factor_2) / denominator,
        -((2*i.dz() - j.dz() - k.dz())/(norm_a*norm_b) + (d_ab*a.dz())/ factor_1 + (d_ab*b.dz())/ factor_2) / denominator
    );

    const vec3 grad_theta_j(
        -((-i.dx() + k.dx()) / (norm_a*norm_b) - (d_ab*a.dx()) / factor_1) / denominator,
        -((-i.dy() + k.dy()) / (norm_a*norm_b) - (d_ab*a.dy()) / factor_1) / denominator,
        -((-i.dz() + k.dz()) / (norm_a*norm_b) - (d_ab*a.dz()) / factor_1) / denominator
    );

    const vec3 grad_theta_k(
        -((-i.dx() + j.dx()) / (norm_a*norm_b) - (d_ab*b.dx()) / factor_2) / denominator,
        -((-i.dy() + j.dy()) / (norm_a*norm_b) - (d_ab*b.dy()) / factor_2) / denominator,
        -((-i.dz() + j.dz()) / (norm_a*norm_b) - (d_ab*b.dz()) / factor_2) / denominator
    );

    return {grad_theta_i, grad_theta_j, grad_theta_k};
}
//---------------------------------------------------------------------------------------------------------


//---------------------------------------------------------------------------------------------------------
void cell::regularize_all_face_angles() noexcept{
    std::for_each(face_lst_.begin(), face_lst_.end(), 
        [&](const face& f) -> void{
            if(f.is_used()) regularize_face_angles(f);
        }
    );
}
//---------------------------------------------------------------------------------------------------------




//---------------------------------------------------------------------------------------------------------
//Apply a force to make sure that the angle of the face are all aproximately 60 degrees
void cell::regularize_face_angles(const face& f) noexcept{
    assert(cell_type_ != nullptr);
    assert(f.is_used());

    if(cell_type_->angle_regularization_factor_ == 0.) return;

    //Get the ids of the nodes of the face
    const auto [n1_id, n2_id, n3_id] = f.get_node_ids();
    assert(n1_id < node_lst_.size() && n2_id < node_lst_.size() && n3_id < node_lst_.size());
    node& n1 = node_lst_[n1_id];
    node& n2 = node_lst_[n2_id];
    node& n3 = node_lst_[n3_id];

    //There are three angles in a triangle, we want them all to be 60 degrees
    const double angle_1 = (n2 - n1).get_angle_with(n3 - n1);
    const double angle_2 = (n1 - n2).get_angle_with(n3 - n2);
    const double angle_3 = (n1 - n3).get_angle_with(n2 - n3);

    //If one of the angle is below 20 deg we don't apply any regularization to the face angles
    constexpr double min_angle = 10  * M_PI / 180.0;
    constexpr double max_angle = 170 * M_PI / 180.0;

    if(angle_1 < min_angle || angle_2 < min_angle || angle_3 < min_angle) return;
    if(angle_1 > max_angle || angle_2 > max_angle || angle_3 > max_angle) return;


    //Compute the gradient of the angles with respect to the coordinates of the nodes
    const auto [grad_angle_1_i, grad_angle_1_j, grad_angle_1_k] = get_angle_gradient(n1.pos(), n2.pos(), n3.pos());
    const auto [grad_angle_2_i, grad_angle_2_j, grad_angle_2_k] = get_angle_gradient(n2.pos(), n1.pos(), n3.pos());
    const auto [grad_angle_3_i, grad_angle_3_j, grad_angle_3_k] = get_angle_gradient(n3.pos(), n1.pos(), n2.pos());

    //Check that there is no nan in the gradients
    if(
        std::isfinite(grad_angle_1_i.dx()) && std::isfinite(grad_angle_1_i.dy()) && std::isfinite(grad_angle_1_i.dz()) &&
        std::isfinite(grad_angle_2_i.dx()) && std::isfinite(grad_angle_2_i.dy()) && std::isfinite(grad_angle_2_i.dz()) &&
        std::isfinite(grad_angle_3_i.dx()) && std::isfinite(grad_angle_3_i.dy()) && std::isfinite(grad_angle_3_i.dz())
    ){
        //Compute the force that is needed to make the angles 60 degrees
        const vec3 force_n1 = ((grad_angle_1_i * ((M_PI / 3.) - angle_1))  + (grad_angle_2_j * ((M_PI / 3.) - angle_2)) + (grad_angle_3_j * ((M_PI / 3.) - angle_3))) * cell_type_->angle_regularization_factor_;
        const vec3 force_n2 = ((grad_angle_1_j * ((M_PI / 3.) - angle_1))  + (grad_angle_2_i * ((M_PI / 3.) - angle_2)) + (grad_angle_3_k * ((M_PI / 3.) - angle_3))) * cell_type_->angle_regularization_factor_;
        const vec3 force_n3 = ((grad_angle_1_k * ((M_PI / 3.) - angle_1))  + (grad_angle_2_k * ((M_PI / 3.) - angle_2)) + (grad_angle_3_i * ((M_PI / 3.) - angle_3))) * cell_type_->angle_regularization_factor_;

        //Apply the force to the nodes
        n1.add_force(force_n1);
        n2.add_force(force_n2);
        n3.add_force(force_n3);
    }
}
//---------------------------------------------------------------------------------------------------------


//---------------------------------------------------------------------------------------------------------
//Translate the cell given a translation vector
void cell::translate(const vec3& translation) noexcept{
    std::for_each(node_lst_.begin(), node_lst_.end(), 
        [&](node& n) -> void{
            if(n.is_used()) n.pos_.translate(translation);
        }
    );
}
//---------------------------------------------------------------------------------------------------------





//---------------------------------------------------------------------------------------------------------
//Indicate on which side of a face a point lies
bool cell::point_is_on_positive_side_of_face(const unsigned face_id, const vec3& point) const noexcept{
    assert(face_id < face_lst_.size());
    assert(face_lst_[face_id].is_used());

    //Get the first node of the face as well as its normal
    const node& n1 = get_const_ref_node(face_lst_[face_id].n1_id_);
    const vec3 face_normal = face_lst_[face_id].get_normal();

    //Compute the vectors from the point to the first node of the face 
    const vec3 vec_n1_to_p = point - n1.pos();

    return face_normal.dot(vec_n1_to_p) > 0.;        
}
//---------------------------------------------------------------------------------------------------------


//---------------------------------------------------------------------------------------------------------
//Get the face made up of the 3 nodes. If the face doesn't exist return a nullptr
face* cell::get_face(const unsigned n1_id, const unsigned n2_id, const unsigned n3_id) noexcept{

    //Get the 3 edges between the nodes
    const auto e1_opt = get_edge(n1_id, n2_id);
    const auto e2_opt = get_edge(n2_id, n3_id);
    const auto e3_opt = get_edge(n3_id, n1_id);

    //If one of the edges doesn't exist, the face doesn't exist
    if(!e1_opt.has_value() || !e2_opt.has_value() || !e3_opt.has_value()) return nullptr;

    //Find the common face between the 3 edges
    const edge& e1 = e1_opt.value();
    const edge& e2 = e2_opt.value();
    const edge& e3 = e3_opt.value();

    assert(e1.f1() < face_lst_.size() && e1.f2() < face_lst_.size());

    if(
        
        (e1.f1() == e2.f1() || e1.f1() == e2.f2()) &&
        (e1.f1() == e3.f1() || e1.f1() == e3.f2())

    ){
        return &(face_lst_[e1.f1()]);

    }
    

    if(
        
        (e1.f2() == e2.f1() || e1.f2() == e2.f2()) &&
        (e1.f2() == e3.f1() || e1.f2() == e3.f2())

    ){
        return &(face_lst_[e1.f2()]);

    }
    
    return nullptr;
}
//---------------------------------------------------------------------------------------------------------
