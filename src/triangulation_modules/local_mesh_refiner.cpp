#include "local_mesh_refiner.hpp"



/*
    This class makes sure that all the edges of the cell meshes have lengths comprised in the range [l_min, l_max]. 
    When an edge is too long, it is subdivided into two edges. When an edge is too short, it is merged with into one node.
    Triangles with very high isoperimetric ratio are also subdivided into two triangles.
*/


//-----------------------------------------------------------------------------------------------
//Trivial constructor
local_mesh_refiner::local_mesh_refiner(const double l_min, const double l_max, const bool enable_edge_swap_operation /*=true*/ ) noexcept: 
    l_min_(l_min), l_max_(l_max), 
    l_min_squared_(l_min * l_min), l_max_squared_(l_max * l_max),
    enable_edge_swap_operation_(enable_edge_swap_operation) 
{
    assert(l_min_ > 0.);
    assert(l_max_ > l_min_);    
}
//-----------------------------------------------------------------------------------------------



//-----------------------------------------------------------------------------------------------
//Refine the meshes of all the cells in the list in parallel
void local_mesh_refiner::refine_meshes(const std::vector<cell_ptr> cell_lst) const noexcept(false){


    //Refine the mesh of the cells in parallel. If an exception occurs, the parallel_exception_handler will rethrow it
    parallel_exception_handler(cell_lst, refine_mesh_func_); //noexcept(false)
}
//-----------------------------------------------------------------------------------------------


//-----------------------------------------------------------------------------------------------
void local_mesh_refiner::refine_mesh(cell_ptr c) const noexcept(false){
    //This method ensures a high mesh quality of the cell surface. The elongated 
    //triangles are removed via edge swaps. The edges longer than l_max are subdivided and 
    //the edges shorter than l_min are merged into a node. 

    //The method cell::rebase() must not be called during this method. Do not call the mesh writer for instance

    //Start by updating the centroid of the cell
    c->update_centroid();

    //Remove all the elongated triangles
    if(enable_edge_swap_operation_) remove_elongated_triangles(c);

    //Create a copy of the edge set
    edge_set edge_to_check_set = c->get_edge_set();

    //Keep track of the number of iterations
    unsigned iteration = 0;

    //Loop over all the edges of the cell until all edges have been checked
    while(edge_to_check_set.size() > 0 && iteration < c->get_edge_set().size()){

        //Pop the first element of the set
        auto e_ab = *edge_to_check_set.begin();
        edge_to_check_set.erase(edge_to_check_set.begin());

        //Get the 2 nodes of the edge
        const node& n_a = c->get_node(e_ab.n1());
        const node& n_b = c->get_node(e_ab.n2());

        //Compute the squared length of the edge
        const double l_ab_squared = (n_a - n_b).squared_norm();
        
        //Check if the edge is too long
        if(l_ab_squared > l_max_squared_){
            
            split_edge(e_ab, c, edge_to_check_set);
            iteration++;
        }

        //If the edge is too short
        else if(l_ab_squared < l_min_squared_){

            //Check that the edge can be merged
            if(can_be_merged(e_ab, c)){

                //Merge the edge into a node
                merge_edge(e_ab, c, edge_to_check_set);
                iteration++;

            }
        }
    }
    if(iteration == c->get_edge_set().size()){
        throw mesh_integrity_exception(std::string("The refinement of the mesh of cell ") + std::to_string(c->get_id()) + std::string(" failed.") +
        std::string(" The simulation is unstable. Reducing the time step might help")
        );
    }



}
//-----------------------------------------------------------------------------------------------





//-----------------------------------------------------------------------------------------------
//Compute the quality score of the triangle and returns the longest edge in the triangle
std::pair<double, edge> local_mesh_refiner::get_triangle_score(cell_ptr c, const face& f) const noexcept{
    assert(f.is_used());

    //Get the nodes composing the face
    auto [n_a_id, n_b_id, n_c_id] = f.get_node_ids();

    const node& n_a = c->get_node(n_a_id);
    const node& n_b = c->get_node(n_b_id);
    const node& n_c = c->get_node(n_c_id);

    //Compute the sqared length of the edges
    const double l_ab = (n_a - n_b).norm();
    const double l_bc = (n_b - n_c).norm();
    const double l_ca = (n_c - n_a).norm();


    //Find the longest edge of the triangle
    std::optional<edge> longest_edge_opt;
    if(l_ab > l_bc){
        if(l_ab > l_ca){
            longest_edge_opt = c->get_edge(n_a_id, n_b_id);
        }
        else{
            longest_edge_opt = c->get_edge(n_c_id, n_a_id);
        }
    }
    else{
        if(l_bc > l_ca){
            longest_edge_opt = c->get_edge(n_b_id, n_c_id);
        }
        else{
            longest_edge_opt = c->get_edge(n_c_id, n_a_id);
        }
    }

    assert(longest_edge_opt.has_value());
    edge longest_edge = longest_edge_opt.value();

    const double squared_perimeter = std::pow(l_ab + l_bc + l_ca, 2);

    //Get the area of the face
    const double area = f.get_area();

    //Compare the isoperimetric ratio of the triangle to the minimum possible isoperimetric ratio
    //to get the quality score of the face
    const double quality_score = q_min_ *  area / squared_perimeter;

    return std::make_pair(quality_score, longest_edge);
}
//-----------------------------------------------------------------------------------------------



//-----------------------------------------------------------------------------------------------
//Remove all the triangles with high isoeperimetric ratio
void local_mesh_refiner::remove_elongated_triangles(cell_ptr c) const noexcept(false){
    if(c->get_face_lst().size() == 0){return;}

    //Loop over the faces
    for(unsigned int i = 0; i < c->get_face_lst().size(); i++){

        //Take and remove the last face of the list
        const face& f_1 = c->get_face_lst()[i];

        if(!f_1.is_used()){continue;}

        //Get the quality score of the face and its longest edge
        auto [triangle_score, longest_edge] = get_triangle_score(c, f_1);

        //If the quality score of the face is too low, we perform an edge swap operation
        if(triangle_score < triangle_score_min_){swap_edge(longest_edge, c);}
    }
}
//-----------------------------------------------------------------------------------------------





//-----------------------------------------------------------------------------------------------
void local_mesh_refiner::swap_edge(edge& e_ab, cell_ptr c) const noexcept(false){
   /*Swap the edge to prevent the formation of triangles with high isoperimetric ratios.
    Those triangles can be a source of numerical instabilities.
    
    BEFORE SWAP :                  
            A 
           /|\     
       5  / | \ 6    
         /  |  \   
       C/ 1 | 2 \D
        \   |   /   
         \  |  /   
        8 \ | / 7    
           \|/      
            B

    AFTER SWAP :                  
            A
           / \     
       5  /   \  6     
         /  3  \   
       C/_______\D
        \       /   
         \  4  /   
          \   /     
       8   \ /   7   
            B  
    */


    //Make sure the node is connected to at least 2 faces
    assert(e_ab.is_manifold());
    

    //Extract the nodes composing the edge
    const node& n_a = c->get_node(e_ab.n1());
    const node& n_b = c->get_node(e_ab.n2());
    
    const unsigned f_1_id = e_ab.f1();
    const unsigned f_2_id = e_ab.f2();

    //Extract the faces connected to the edge
    const face& f_1 = c->get_face(e_ab.f1());
    const face& f_2 = c->get_face(e_ab.f2());

    assert(n_a.is_used() && n_b.is_used());
    assert(f_1.is_used() && f_2.is_used());
    
    //Get the opposite nodes of the edge in the two faces
    const unsigned id_n_c = f_1.get_opposite_node(n_a.get_local_id(), n_b.get_local_id());
    const unsigned id_n_d = f_2.get_opposite_node(n_a.get_local_id(), n_b.get_local_id());

    const node& n_c = c->get_node(id_n_c);
    const node& n_d = c->get_node(id_n_d);

    //Delete the face 1 and 2 which should have for effect to delete the 
    //edge e_ab at the same time

    //Get a copy of the edges AC, CB, BD, and DA
    auto e_ac_opt = c->get_edge(n_a.get_local_id(), n_c.get_local_id());
    auto e_cb_opt = c->get_edge(n_c.get_local_id(), n_b.get_local_id());
    auto e_bd_opt = c->get_edge(n_b.get_local_id(), n_d.get_local_id());
    auto e_da_opt = c->get_edge(n_d.get_local_id(), n_a.get_local_id());

    assert(e_ac_opt.has_value() && e_cb_opt.has_value() && e_bd_opt.has_value() && e_da_opt.has_value());
    edge e_ac = e_ac_opt.value();
    edge e_cb = e_cb_opt.value();
    edge e_bd = e_bd_opt.value();
    edge e_da = e_da_opt.value();

    const unsigned f_5_id = (e_ac.f1() == f_1_id) ? e_ac.f2() : e_ac.f1();
    const unsigned f_8_id = (e_cb.f1() == f_1_id) ? e_cb.f2() : e_cb.f1();
    const unsigned f_7_id = (e_bd.f1() == f_2_id) ? e_bd.f2() : e_bd.f1();
    const unsigned f_6_id = (e_da.f1() == f_2_id) ? e_da.f2() : e_da.f1();

    //If we are in a pathologic configuration, we do not swap the edge
    if(f_5_id == f_6_id || f_7_id == f_8_id){return;}

    //Make sure the edge does not exist already
    if(c->get_edge(n_c.get_local_id(), n_d.get_local_id()).has_value() == true){return;}

    //Make sure they are connected to the correct faces
    assert(e_ac.has_face(f_1_id) && e_ac.has_face(f_5_id));
    assert(e_cb.has_face(f_1_id) && e_cb.has_face(f_8_id));
    assert(e_bd.has_face(f_2_id) && e_bd.has_face(f_7_id));
    assert(e_da.has_face(f_2_id) && e_da.has_face(f_6_id));

    //Delete the faces 1 and 2 from the cell
    c->delete_face(f_1.get_local_id());
    c->delete_face(f_2.get_local_id());

    //Make sure the edge AB does not exist anymore
    assert(!c->get_edge(n_a.get_local_id(), n_b.get_local_id()).has_value());

    //Get a copy of the edges AC, CB, BD, and DA and check that they are non manifold
    e_ac_opt = c->get_edge(n_a.get_local_id(), n_c.get_local_id());
    e_cb_opt = c->get_edge(n_c.get_local_id(), n_b.get_local_id());
    e_bd_opt = c->get_edge(n_b.get_local_id(), n_d.get_local_id());
    e_da_opt = c->get_edge(n_d.get_local_id(), n_a.get_local_id());

    assert(e_ac_opt.has_value() && e_cb_opt.has_value() && e_bd_opt.has_value() && e_da_opt.has_value());
    e_ac = e_ac_opt.value();
    e_cb = e_cb_opt.value();
    e_bd = e_bd_opt.value();
    e_da = e_da_opt.value();

    //Make sure they are connected to the correct faces
    assert(e_ac.is_manifold() == false);
    assert(e_cb.is_manifold() == false);
    assert(e_bd.is_manifold() == false);
    assert(e_da.is_manifold() == false);

    //The local ids of the new faces
    unsigned f_3_id, f_4_id;

    //Swap the edges such that the normals are pointing in the same direction
    f_3_id = c->create_face(n_a.get_local_id(), n_d.get_local_id(), n_c.get_local_id());
    f_4_id = c->create_face(n_b.get_local_id(), n_c.get_local_id(), n_d.get_local_id());

    //We need to make sure that the winding order of the face 3 and 4 are correct.
    face& f_3 = c->face_lst_[f_3_id];
    face& f_4 = c->face_lst_[f_4_id];
    const face& ref_f_5 = c->get_face(f_5_id);
    const face& ref_f_8 = c->get_face(f_8_id);
    c->check_face_winding_order(ref_f_5, f_3);
    c->check_face_winding_order(ref_f_8, f_4);

    //Get a copy of the edges AC, CB, BD, and DA and check that they are connected to the correct faces
    e_ac_opt = c->get_edge(n_a.get_local_id(), n_c.get_local_id());
    e_cb_opt = c->get_edge(n_c.get_local_id(), n_b.get_local_id());
    e_bd_opt = c->get_edge(n_b.get_local_id(), n_d.get_local_id());
    e_da_opt = c->get_edge(n_d.get_local_id(), n_a.get_local_id());

    assert(e_ac_opt.has_value() && e_cb_opt.has_value() && e_bd_opt.has_value() && e_da_opt.has_value());
    e_ac = e_ac_opt.value();
    e_cb = e_cb_opt.value();
    e_bd = e_bd_opt.value();
    e_da = e_da_opt.value();

    //Make sure they are connected to the correct faces
    assert(e_ac.has_face(f_3_id) && e_ac.has_face(f_5_id));
    assert(e_cb.has_face(f_4_id) && e_cb.has_face(f_8_id));
    assert(e_bd.has_face(f_4_id) && e_bd.has_face(f_7_id));
    assert(e_da.has_face(f_3_id) && e_da.has_face(f_6_id));



}
//-----------------------------------------------------------------------------------------------


//-----------------------------------------------------------------------------------------------
void local_mesh_refiner::split_edge(edge& e_ab, cell_ptr c, edge_set& edge_to_check_set) const noexcept(false){

    /*When an edge is too long, it is split in two.

        BEFORE SPLIT :                  
            A
           /|\     
          / | \     
         /  |  \   
       C/ 1 | 2 \D
        \   |   /   
         \  |  /   
          \ | /     
           \|/      
            B

    AFTER SPlIT :                  
            A
           /|\     
          / | \     
         / 3|4 \   
       C/___|E__\D
        \   |   /   
         \ 5|6 /   
          \ | /     
           \|/      
            B  
    */


    //Make sure the node is connected to at least 2 faces
    assert(e_ab.is_manifold());

    //Extract the nodes composing the edge
    node& n_a = c->get_node(e_ab.n1());
    node& n_b = c->get_node(e_ab.n2());

    const unsigned id_n_a = n_a.get_local_id();
    const unsigned id_n_b = n_b.get_local_id();

    //Extract the faces connected to the edge
    const  unsigned f_1_id = e_ab.f1();
    const  unsigned f_2_id = e_ab.f2();

    const face& f_1 = c->get_face(e_ab.f1());
    const face& f_2 = c->get_face(e_ab.f2());

    //Get the normals of the 2 faces
    const vec3 f1_normal = f_1.get_normal();
    const vec3 f2_normal = f_2.get_normal();

    //Get the opposite nodes of the edge in the two faces
    const unsigned id_n_c = f_1.get_opposite_node(n_a.get_local_id(), n_b.get_local_id());
    const unsigned id_n_d = f_2.get_opposite_node(n_a.get_local_id(), n_b.get_local_id());
    const node& n_c = c->get_node(id_n_c);
    const node& n_d = c->get_node(id_n_d);

    //Create a node in the middle of the edge
    const vec3 n_e_pos = (n_b.pos() + n_a.pos()) * 0.5;
    node n_e(n_e_pos, 0);


    //Full equations of motion solved with the semi implicit euler scheme
    #if DYNAMIC_MODEL_INDEX == 0 

        //The nodes a and b distribute a third of their momentum to the new node e
        const vec3 n_a_momentum = n_a.momentum();
        const vec3 n_b_momentum = n_b.momentum();
        n_a.set_momentum(n_a_momentum * (2./3.));
        n_b.set_momentum(n_b_momentum * (2./3.));
        n_e.set_momentum((n_a_momentum + n_b_momentum) / 3.0);

    #elif DYNAMIC_MODEL_INDEX == 2 //Overdamped equations of motion solved with the improved euler scheme

        //The nodes a and b distribute their previous force to the new node e
        const vec3 n_a_previous_force = n_a.previous_force();
        const vec3 n_b_previous_force = n_b.previous_force();
        n_a.set_previous_force(n_a_previous_force * (2./3.));
        n_b.set_previous_force(n_b_previous_force * (2./3.));
        n_e.set_previous_force((n_a_previous_force + n_b_previous_force) / 3.0);

        //The previous position of the node e is the average of the previous positions of the nodes a and b
        n_e.set_previous_pos((n_a.previous_pos() + n_b.previous_pos()) / 2.0);
    #endif 


    //Add the new node to the cell
    const unsigned id_n_e = c->add_node(n_e);

    //Use the vector from the cell centroid to the new node to determine the normal of the new faces

    //Get the 2 face normals before their deletion
    const vec3 ref_normal = n_e.pos() - c->get_centroid();

    //Get the local face type id of the face 1 and 2
    const unsigned short f_1_type_id = f_1.get_local_face_type_id();
    const unsigned short f_2_type_id = f_2.get_local_face_type_id();


    //Delete the faces 1 and 2 which should have for effect to delete the edge e_ab at the same time
    c->delete_face(f_1.get_local_id());
    c->delete_face(f_2.get_local_id());

    //The local ids of the new faces that will be created
    unsigned f_3_id, f_4_id, f_5_id, f_6_id;


    //Create the faces 3 and 5 such that they have the normal pointing in
    //the same direction as the normal of f1
    if ((n_a - n_c).cross(n_b - n_c).dot(f1_normal) >= 0.){
        f_3_id = c->create_face(id_n_c, id_n_a, id_n_e);
        f_5_id = c->create_face(id_n_c, id_n_e, id_n_b);
    }else{
        f_3_id = c->create_face(id_n_c, id_n_e, id_n_a);
        f_5_id = c->create_face(id_n_c, id_n_b, id_n_e);

    }


    //Create the faces 4 and 6 such that they have the normal pointing in
    //the same direction as the normal of f2
    if ((n_a - n_d).cross(n_b - n_d).dot(f2_normal) >= 0.){
        f_4_id = c->create_face(id_n_d, id_n_a, id_n_e);
        f_6_id = c->create_face(id_n_d, id_n_e, id_n_b);
    }else{
        f_4_id = c->create_face(id_n_d, id_n_e, id_n_a);
        f_6_id = c->create_face(id_n_d, id_n_b, id_n_e);
    }


    //Set the face type id of the new faces to the same as the old ones
    c->get_face(f_3_id).set_face_type_id(f_1_type_id);
    c->get_face(f_4_id).set_face_type_id(f_2_type_id);
    c->get_face(f_5_id).set_face_type_id(f_1_type_id);
    c->get_face(f_6_id).set_face_type_id(f_2_type_id);

    //Get the newly created edge 
    auto edge_opt_ea = c->get_edge(id_n_e, id_n_a);
    auto edge_opt_eb = c->get_edge(id_n_e, id_n_b);
    auto edge_opt_ec = c->get_edge(id_n_e, id_n_c);
    auto edge_opt_ed = c->get_edge(id_n_e, id_n_d);
    assert(edge_opt_ea.has_value() && edge_opt_eb.has_value() && edge_opt_ec.has_value() && edge_opt_ed.has_value());

    const edge& e_ea = edge_opt_ea.value();
    const edge& e_eb = edge_opt_eb.value();
    const edge& e_ec = edge_opt_ec.value();
    const edge& e_ed = edge_opt_ed.value();
    assert(e_ea.is_manifold() && e_eb.is_manifold() && e_ec.is_manifold() && e_ed.is_manifold());

    //Add the newly created edge to the set of edges that have to be checked
    edge_to_check_set.insert({e_ea, e_eb, e_ec, e_ed});

    //We need to update the edge AC, AD, BC, and BD if they are still in the set of edges to check
    auto edge_i_ac = edge_to_check_set.find(edge(id_n_a, id_n_c));
    auto edge_i_ad = edge_to_check_set.find(edge(id_n_a, id_n_d));
    auto edge_i_bc = edge_to_check_set.find(edge(id_n_b, id_n_c));
    auto edge_i_bd = edge_to_check_set.find(edge(id_n_b, id_n_d));

    //Replace the face f1 by the faces f3 and f5
    //Replace the face f2 by the faces f4 and f6
    if(edge_i_ac != edge_to_check_set.end()){const_cast<edge&>(*edge_i_ac).replace_face(f_1_id, f_3_id);}
    if(edge_i_bc != edge_to_check_set.end()){const_cast<edge&>(*edge_i_bc).replace_face(f_1_id, f_5_id);}
    if(edge_i_ad != edge_to_check_set.end()){const_cast<edge&>(*edge_i_ad).replace_face(f_2_id, f_4_id);}
    if(edge_i_bd != edge_to_check_set.end()){const_cast<edge&>(*edge_i_bd).replace_face(f_2_id, f_6_id);}


}
//-----------------------------------------------------------------------------------------------


//-----------------------------------------------------------------------------------------------
bool local_mesh_refiner::can_be_merged(edge& e_ab, cell_ptr c) const noexcept(false){

/*
    In a rare cases, the arrangement of the nodes prevents the merger operation. These specific regions looks like that
       G___ A_
       |   /|\＼     
       | 6/ | \ ＼    
       | /  |  \ 3＼  
       E/ 1 | 2 \D__＼C
       |\   |   /   ／
       | \  |  / 4／  
       |  \ | / ／    
       | 5 \|/／      
      F|___ B
    
    Merging the edge between node A and B would form a new node X and also join together the triangles (ACD) and (BCD). 
    The new edge X-C would be shared by more than 2 faces which is a pathological case leading to a crash. 

    The following code detects this type of configurations and return false if it is the case.
 */   

    //Make sure the node is connected to at least 2 faces
    assert(e_ab.is_manifold());

    //Extract the nodes composing the edge
    const node& n_a = c->get_node(e_ab.n1());
    const node& n_b = c->get_node(e_ab.n2());

    const face& f_a = c->get_face(e_ab.f1());
    const face& f_b = c->get_face(e_ab.f2());

    //Get all the faces connected to the nodes a and b
    std::vector<unsigned> node_lst_a = c->get_connected_nodes(n_a.get_local_id(), e_ab);
    std::vector<unsigned> node_lst_b = c->get_connected_nodes(n_b.get_local_id(), e_ab);

    std::sort(node_lst_a.begin(), node_lst_a.end());
    std::sort(node_lst_b.begin(), node_lst_b.end());

    std::set<unsigned> node_set_ab;
    std::set_intersection(node_lst_a.begin(),node_lst_a.end(),node_lst_b.begin(),node_lst_b.end(),
    std::inserter(node_set_ab, node_set_ab.begin()));

    assert(node_set_ab.size() >= 2);
    return node_set_ab.size() == 2;
}

//-----------------------------------------------------------------------------------------------


//-----------------------------------------------------------------------------------------------
void local_mesh_refiner::merge_edge(edge& e_ab, cell_ptr c, edge_set& edge_to_check_set) const noexcept(false){

    /*When an edge is too short, it is merged into a node

    BEFORE MERGER: 

       E____C___F
       |   /\   |  
       |  /  \  |   
       | /    \ |  
       A/______\B
       |\      /|   
       | \    / |  
       |  \  /  |   
       |___\/___|   
       H    D    G

    AFTER MERGER:

       E____C____F
        \   |   /
         \  |  /
          \ | /
           \|/
            I
           /|\
          / | \      
         /  |  \    
        /   |   \    
       /___ |____\   
       H    D    G

    */


    //Make sure the node is connected to at least 2 faces
    assert(e_ab.is_manifold());

    //Extract the nodes composing the edge
    const unsigned id_n_a = e_ab.n1();
    const unsigned id_n_b = e_ab.n2();

    node& n_a = c->get_node(id_n_a);
    node& n_b = c->get_node(id_n_b);
    
    //Get the ids of the 2 faces that are connected to the edge
    const unsigned id_f_1 = e_ab.f1();
    const unsigned id_f_2 = e_ab.f2();


    //The node i will be inserted at the middle of the edge
    const vec3 n_i_pos = (n_b.pos() + n_a.pos()) * 0.5;

    //Create the new node
    node n_i(n_i_pos, 0);
    
    //Full equations of motion solved with the semi implicit euler scheme
    #if DYNAMIC_MODEL_INDEX == 0 

        //The momentum of the node A and B will be distributed entirely to the new node I
        const vec3 n_i_momentum = n_a.momentum() + n_b.momentum();
        n_i.set_momentum(n_i_momentum);

    //Overdamped equations of motion solved with the improved euler scheme
    #elif DYNAMIC_MODEL_INDEX == 2 

        //The node i inherit from the previous force of the nodes a and b
        n_i.set_previous_pos((n_a.previous_pos() + n_b.previous_pos()) * 0.5);
        n_i.set_previous_force(n_a.previous_force() + n_b.previous_force());
    #endif 



    //Add the new node to the cell
    const unsigned id_n_i = c->add_node(n_i);

    //Replace the node A and B by the new node I in the edges
    auto [deleted_edge_a_lst, created_edge_a_lst]  = c->replace_node(e_ab, id_n_a, id_n_i);

    //The edge e_ab has been deleted, so we need to get the new edge e_bi
    auto edge_opt_bi = c->get_edge(id_n_b, id_n_i);
    assert(edge_opt_bi.has_value());
    const edge& e_bi = edge_opt_bi.value();

    auto [deleted_edge_b_lst, created_edge_b_lst] = c->replace_node(e_bi, id_n_b, id_n_i);

    //The nodes A and B have already been deleted by the replace node operation
    //Delete the faces 1 and 2 which should have for effect to delete the edge e_ii which was created
    //after the 2 consecutive cell::replace_node() calls 
    c->delete_face(id_f_1);
    c->delete_face(id_f_2);


    //Remove all the edges that have been deleted by the cell::replace_node() methods from edge_to_check_set
    std::for_each(deleted_edge_a_lst.begin(), deleted_edge_a_lst.end(), [&](const edge& e){
        auto edge_to_erase_it = edge_to_check_set.find(e);
        if(edge_to_erase_it != edge_to_check_set.end()){edge_to_check_set.erase(e);}
    });

    std::for_each(deleted_edge_b_lst.begin(), deleted_edge_b_lst.end(), [&](const edge& e){
        auto edge_to_erase_it = edge_to_check_set.find(e);
        if(edge_to_erase_it != edge_to_check_set.end()){edge_to_check_set.erase(e);}
    });

    //Insert all the edges that have been created by the cell::replace_node() methods into edge_to_check_set
    std::for_each(created_edge_a_lst.begin(), created_edge_a_lst.end(), [&](const edge& e){

        if(!e.has_node(id_n_a) && !e.has_node(id_n_b)  && e.n1() != e.n2() &&!e.has_face(id_f_1) && !e.has_face(id_f_2)){
            edge_to_check_set.insert(e);
        }
    });

    std::for_each(created_edge_b_lst.begin(), created_edge_b_lst.end(), [&](const edge& e){
        if(!e.has_node(id_n_a) && !e.has_node(id_n_b)  && e.n1() != e.n2() &&!e.has_face(id_f_1) && !e.has_face(id_f_2)){
            edge_to_check_set.insert(e);
        }
    });

}
//-----------------------------------------------------------------------------------------------


