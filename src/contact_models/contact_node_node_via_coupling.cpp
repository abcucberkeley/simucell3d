#include "contact_node_node_via_coupling.hpp"




#if CONTACT_MODEL_INDEX == 1

/*
This class is a derived class of the contact_model_node_face_ and has the advantage that it also works when the cell meshes are not locally refined 
during the simulations. The mechanic is exactly the same, but this version of the face-node contact model is less optimized than its
base class since the size of the triangles is unpredictable when the local mesh refiner is turned off.
*/

//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//The trivial constructor of the contact model
contact_node_node_via_coupling::contact_node_node_via_coupling(const global_simulation_parameters& sim_parameters) noexcept(false): contact_model_abstract(sim_parameters){
}
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------





//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//Find the contacts between the faces of the mesh and apply the contact forces
void contact_node_node_via_coupling::run(const std::vector<cell_ptr>& cell_lst) noexcept{

        //Store all the faces of the mesh in a vector for easier access
        const size_t nb_faces = std::accumulate(cell_lst.begin(), cell_lst.end(), 0, [](size_t acc, const cell_ptr& c){return acc + c->get_nb_of_faces();});
        face_lst_.clear();
        face_lst_.reserve(nb_faces);
        
        //Give to each face a global id which is the index of the face in the face_lst_ vector
        size_t face_global_id = 0;
        for(cell_ptr c : cell_lst){
            for(auto& f : c->face_lst_){
                if(f.is_used()){
                    face_lst_.push_back(&f);
                    
                    //Set the global id of all the faces
                    f.global_face_id_ = face_global_id++;

                    //Only compiles this section if the faces store their contact energies
                    #if FACE_STORE_CONTACT_ENERGY
                        //Reset the adhesion and repulsion energies of the face
                        f.adhesion_energy_ = 0.0;
                        f.repulsion_energy_ = 0.0;
                    #endif
                }
            }

            for(node& n: c->node_lst_){
                if(n.is_used()){
                    //Reset the coupling of the node
                    n.coupled_node_ = std::nullopt;
                    n.squared_distance_to_closest_node_ = std::numeric_limits<double>::max();
                }
            }
        }

        //Compute the axis aligned bounding box of ech face and store it in the face_aabb_lst_ vector and update the grid dimensions
        update_face_aabbs(cell_lst);

        //Store all the faces of the mesh in a space partitionning grid
        store_face_in_uspg();
    
        //Find the faces that are within a distance below the contact cutoff and apply adhesive or repulsive forces
        resolve_all_contacts(cell_lst);
}
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------





//-----------------------------------------------------------------------------------------------
void contact_node_node_via_coupling::resolve_all_contacts(const std::vector<cell_ptr>& cell_lst) noexcept{

    //Loop over the cells in parallel
    #pragma omp parallel for
    for(size_t cell_id = 0; cell_id < cell_lst.size(); cell_id++){
        cell_ptr c1 = cell_lst[cell_id];

        //Get the maximum mean curvature of the cell, above which on eof its point cannot be coupled to a face
        //on another cell surface
        const double surface_coupling_max_curvature = c1->get_cell_type()->surface_coupling_max_curvature_;

        //Loop over the nodes of the cell
        for(node& n: c1->node_lst_){

            if(n.is_used() && n.curvature_ < surface_coupling_max_curvature){

                //Get the position of the node in the space partitioning grid
                const unsigned voxel_1_x = std::floor((n.pos().dx() - grid_.min_x_) / grid_.voxel_size_);
                const unsigned voxel_1_y = std::floor((n.pos().dy() - grid_.min_y_) / grid_.voxel_size_);
                const unsigned voxel_1_z = std::floor((n.pos().dz() - grid_.min_z_) / grid_.voxel_size_);

                //Get the ID of the voxel in the space partitionning grid
                const size_t voxel_id = grid_.get_voxel_index(voxel_1_x, voxel_1_y, voxel_1_z);
                assert(voxel_id < grid_.voxel_lst_.size());

                //Get the list of faces stored in this voxel
                for(face* f: grid_.voxel_lst_[voxel_id]){
                    assert(f != nullptr);
                    cell_ptr c2 = f->get_owner_cell();

                    if(c1->get_id() != c2->get_id()){

                        //Check if the node is located in the AABB of the face
                        if(
                            
                            aabb_intersection_check(f->global_face_id_ * 6, n.pos()) &&
                            n.normal_.dot(f->normal_) < max_dot_product_repulsion_
                        ){
                            resolve_contact(c1, c2, n, f);
                        }
                    }
                }
            }  
        }
    }



    for(size_t c1_id = 0; c1_id < cell_lst.size(); c1_id++){
        cell_ptr c1 = cell_lst[c1_id];

        //Loop over the nodes of the cell
        for(node& n1: c1->node_lst_){

                if(n1.is_used()){

                //If the node has been coupled to another node
                if(n1.coupled_node_.has_value()){
                    const auto [c2_id, n2_id] = n1.coupled_node_.value();

                    if(c1_id > c2_id){
                        assert(c2_id < cell_lst.size());

                        //Get the coupled node
                        cell_ptr c2 = cell_lst[c2_id];
                        
                        assert(n2_id < c2->node_lst_.size());
                        node& n2 = c2->node_lst_[n2_id];

                        //Find the center point position between the 2 nodes and move them at this location
                        const vec3 center_point = (n1.pos() + n2.pos()) * 0.5;

                        //Move the node n1 to the position of the coupled node
                        n1.pos_.reset(center_point);
                        n2.pos_.reset(center_point);
                    }
                }
            }   
        }
    }
}
//-----------------------------------------------------------------------------------------------










//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//Prevent the surfaces of the cells to interpenetrate by applying repulsive forces on the surfaces
void contact_node_node_via_coupling::resolve_contact(cell_ptr c1, cell_ptr c2, node& n1, face* f) const noexcept{

    //Make sure that the nodes of the face are part of the cell number 2
    assert(f->n1_id_ < c2->node_lst_.size() && f->n2_id_ < c2->node_lst_.size() && f->n3_id_ < c2->node_lst_.size());

    //Get the node of the face
    node& f_n1 = c2->node_lst_[f->n1_id_];
    node& f_n2 = c2->node_lst_[f->n2_id_];
    node& f_n3 = c2->node_lst_[f->n3_id_];

    assert(f_n1.is_used() && f_n2.is_used() && f_n3.is_used());
    
    //Extract the position of the nodes of the face
    const vec3& f_n1_pos = f_n1.pos();
    const vec3& f_n2_pos = f_n2.pos();
    const vec3& f_n3_pos = f_n3.pos();


    //If the 2 cells are epithelial cell
    if(
        c1->get_cell_type_id() == 0 && c2->get_cell_type_id() == 0
    ){

        //Compute the squared distances between the node n1 and the nodes of the face
        double squared_distance_f_n1 = (n1.pos() - f_n1_pos).squared_norm();
        double squared_distance_f_n2 = (n1.pos() - f_n2_pos).squared_norm();
        double squared_distance_f_n3 = (n1.pos() - f_n3_pos).squared_norm();

        //If the normal at the node n1 is pointing in an opposite direction than the normal of the face's nodes, we set the squared distance to infinity
        if(n1.normal_.dot(f_n1.normal_) > max_dot_product_adhesion_){squared_distance_f_n1 = std::numeric_limits<double>::max();}
        if(n1.normal_.dot(f_n2.normal_) > max_dot_product_adhesion_){squared_distance_f_n2 = std::numeric_limits<double>::max();}
        if(n1.normal_.dot(f_n3.normal_) > max_dot_product_adhesion_){squared_distance_f_n3 = std::numeric_limits<double>::max();}

        //If the node in the face are already coupled to another node, and their distance with this other coupled 
        //node is smaller than the distance with the node n1, we also set the squared distance to infinity
        if(f_n1.coupled_node_.has_value() && f_n1.squared_distance_to_closest_node_ < squared_distance_f_n1){squared_distance_f_n1 = std::numeric_limits<double>::max();}
        if(f_n2.coupled_node_.has_value() && f_n2.squared_distance_to_closest_node_ < squared_distance_f_n2){squared_distance_f_n2 = std::numeric_limits<double>::max();}
        if(f_n3.coupled_node_.has_value() && f_n3.squared_distance_to_closest_node_ < squared_distance_f_n3){squared_distance_f_n3 = std::numeric_limits<double>::max();}


        //Make sure the curvature at the other point is below the threshold
        const double surface_coupling_max_curvature = c1->get_cell_type()->surface_coupling_max_curvature_;
        if(f_n1.curvature_ > surface_coupling_max_curvature){squared_distance_f_n1 = std::numeric_limits<double>::max();}
        if(f_n2.curvature_ > surface_coupling_max_curvature){squared_distance_f_n2 = std::numeric_limits<double>::max();}
        if(f_n3.curvature_ > surface_coupling_max_curvature){squared_distance_f_n3 = std::numeric_limits<double>::max();}


        //Find the node of the face that is the closest to the node n1
        double min_squared_distance_to_node_face;
        node* n2 = nullptr;

        if(squared_distance_f_n1 < squared_distance_f_n2 && squared_distance_f_n1 < squared_distance_f_n3){
            min_squared_distance_to_node_face = squared_distance_f_n1;
            n2 = &f_n1;
        }
        else if(squared_distance_f_n2 < squared_distance_f_n1 && squared_distance_f_n2 < squared_distance_f_n3){
            min_squared_distance_to_node_face = squared_distance_f_n2;
            n2 = &f_n2;
        }
        else{
            min_squared_distance_to_node_face = squared_distance_f_n3;
            n2 = &f_n3;
        }


        //If this closesst node in the face (n2) is closer to the node n1 than any previous closest node
        //and if the distance is below the adhesion cutoff, then we mechanically couple the node n1 to the node n2
        if(                                                                              
            min_squared_distance_to_node_face < interaction_cutoff_square_adhesion_ &&      //If the distance is below the adhesion cutoff
            min_squared_distance_to_node_face < n1.squared_distance_to_closest_node_       //If the distance is smaller than the previous closest distance
        ){
            assert(n2->get_local_id() < c2->get_node_lst().size());
            assert(n2->is_used());
            assert(c2->get_node_lst()[n2->get_local_id()].is_used());

            //This section is thread safe
            n1.set_coupled_node_and_min_distance(std::make_pair(c2->get_local_id(), n2->get_local_id()), min_squared_distance_to_node_face);
            n2->set_coupled_node_and_min_distance(std::make_pair(c1->get_local_id(), n1.get_local_id()), min_squared_distance_to_node_face);

            //The contact has been resolved, we can stop here
            return;
        }
    }

    //If the 2 cells are not epithelial cells, or if the node n1 has not been coupled to one of the node of the face
    //We compute the minimum distance between the node and the face


    //Get the minimal squared distance between the node and the face
    const auto [min_squared_distance, bary_pos] = compute_node_triangle_distance(n1.pos(), f_n1_pos, f_n2_pos, f_n3_pos);


    //If the node is located within the interaction cutoff distance of the face, we have to check if a repulsion force need to be applied.
    if(
        min_squared_distance < max_interaction_cutoff_square_
    ){

        //The contact force vector goes from the closest point of approach on the face to the node
        const vec3 cpa_pos = f_n1_pos * bary_pos.dx() + f_n2_pos * bary_pos.dy() + f_n3_pos * bary_pos.dz();
        const vec3 cpa_f_to_node = n1.pos() - cpa_pos;

        //Get a reference to the types of the face
        const auto& face_type = c2->get_face_type(f->local_face_id_);

        //Determine if the contact is of adhesive or repulsive type based on the normal of the face
        bool contact_is_repulsive = (cpa_f_to_node.dot(f->normal_) < 0.0);

        //Reverse the direction of contact if the cell to which the face belongs is an ECM cell
        if (c1->get_cell_type_id() == 0 && c2->get_cell_type_id() == 1){contact_is_repulsive = !contact_is_repulsive;}

        //Do the same of the node belongs to a nucleus and the face to an epithelial cell
        if (c1->get_cell_type_id() == 3 && c2->get_cell_type_id() == 0){contact_is_repulsive = !contact_is_repulsive;}


        //If the 2 membranes are interpenetrating, we apply a repulsive force
        if(contact_is_repulsive){
            
            //Compute the distance
            const double min_distance = std::sqrt(min_squared_distance);
            const double integration_region = f->get_area();

            //Get a reference to the types of the face
            const auto& face_type = c2->get_face_type(f->local_face_id_);

            #if FACE_STORE_CONTACT_ENERGY
                const double repulsion_energy = 0.5 * face_type.repulsion_strength_ * integration_region * min_distance * min_distance;
                f->add_repulsion_energy(repulsion_energy);
            #endif

            //The force that will be distributed on the nodes of the 2 faces
            const vec3 repulsion_force_vector =  cpa_f_to_node  *  face_type.repulsion_strength_ * integration_region;

            //Apply the repulsion forces on the nodes of the 2 faces
            f_n1.add_force(repulsion_force_vector * bary_pos.dx());
            f_n2.add_force(repulsion_force_vector * bary_pos.dy());
            f_n3.add_force(repulsion_force_vector * bary_pos.dz());
            n1.add_force(repulsion_force_vector * -1.0);   
        }
    }
}
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------




#endif