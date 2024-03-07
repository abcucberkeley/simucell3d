#include "contact_node_face_via_spring.hpp"


#if CONTACT_MODEL_INDEX == 0


/*
    This class contains the methods to compute the contact forces between pairs of nodes and faces belonging to adjacent cells. 
    The contact forces are calculated in the following way: 

    1 -> All the axias-aligned bounding boxes of the faces are computed and stored in a vector
    2 -> The faces are stored in an uniform space partitioning grid (USPG) based on their AABB
    3 -> For each point, we used the USPG to find the faces that are close to it
    4 -> We check if the point is within the AABB of the face
    5 -> If so, we compute the minimal distance between the point and the face
    6 -> The contact force is the computed based on the minimum distance and the contact parameters of the face

*/



//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//The trivial constructor of the contact model
contact_node_face_via_spring::contact_node_face_via_spring(const global_simulation_parameters& sim_parameters) noexcept(false): contact_model_abstract(sim_parameters){

    //There are 2 interaction cutoffs, one for the adhesion and one for the repulsion. We use the max of the two to setup
    //the grid and the AABBs
    interaction_cutoff_ = std::max(sim_parameters.contact_cutoff_adhesion_, sim_parameters.contact_cutoff_repulsion_);
    interaction_cutoff_square_ = interaction_cutoff_*interaction_cutoff_;

    //The distance between 2 faces where the contact force is maximal, by default it's equal to interaction_cutoff_ / 2
    hardening_distance_ = sim_parameters.contact_cutoff_adhesion_ / 2.0;

}
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//Find the contacts between the faces of the mesh and apply the contact forces
void contact_node_face_via_spring::run(const std::vector<cell_ptr>& cell_lst) noexcept{

        //Store all the faces of the mesh in a vector for easier access
        const size_t nb_faces = std::accumulate(cell_lst.begin(), cell_lst.end(), 0, [](size_t acc, const cell_ptr& c){return acc + c->get_nb_of_faces();});
        face_lst_.clear();
        face_lst_.reserve(nb_faces);
        
        //Give to each face a global id which is the index of the face in the face_lst_ vector
        size_t face_global_id = 0;
        for(cell_ptr c : cell_lst){
            for(auto& f : c->face_lst_){if(f.is_used()){
                face_lst_.push_back(&f);
                
                //Set the global id of all the faces
                f.global_face_id_ = face_global_id++;

                //Only compiles this section if the faces store their contact energies
                #if FACE_STORE_CONTACT_ENERGY
                    //Reset the adhesion and repulsion energies of the face
                    f.adhesion_energy_ = 0.0;
                    f.repulsion_energy_ = 0.0;
                #endif
    
            }}
        }

        //Compute the axis aligned bounding box of ech face and store it in the face_aabb_lst_ vector and update the grid dimensions
        update_face_aabbs(cell_lst);

        //Store all the faces of the mesh in a space partitionning grid
        store_face_in_uspg();
    
        //Find the faces that are within a distance below the contact cutoff and apply adhesive or repulsive forces
        resolve_contacts(cell_lst);
}
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------









//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//Store all the faces in an unform space partitioning grid
void contact_node_face_via_spring::resolve_contacts(const std::vector<cell_ptr>& cell_lst) noexcept{


    //Loop over the cells in parallel
    #pragma omp parallel for
    for(size_t cell_id = 0; cell_id < cell_lst.size(); cell_id++){
        cell_ptr c1 = cell_lst[cell_id];

        //Loop over the nodes of the cell
        for(node& n: c1->node_lst_){

            if(n.is_used()){
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

                    if(c1->get_id() != f->get_owner_cell()->get_id()){
 
                        //Check if the node is located in the AABB of the face
                        const size_t face_aabb_pos = f->global_face_id_ * 6;
                        if(aabb_intersection_check(face_aabb_pos, n.pos())){apply_contact_forces(c1, n, f);}
                    }
                }

            }  
        }
    }   
}
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------







//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//Compute the distance between the 2 faces and apply the contact forces
void contact_node_face_via_spring::apply_contact_forces(cell_ptr c1, node& n, face* f) const noexcept{

    //Get the cell owning the face
    cell_ptr c2 = f->get_owner_cell();

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

    //Get the minimal squared distance between the node and the face
    const auto [min_squared_distance, bary_pos] = contact_model_abstract::compute_node_triangle_distance(n.pos(), f_n1_pos, f_n2_pos, f_n3_pos);

  
    //Now that we know the minimal distance, we compute the resulting contact forces
    //If the distance between the 2 faces is below the interaction cutoff
    if(min_squared_distance < interaction_cutoff_square_ && min_squared_distance != 0.0){

        //The contact force vector goes from the closest point of approach on the face to the node
        const vec3 cpa_pos = f_n1_pos * bary_pos.dx() + f_n2_pos * bary_pos.dy() + f_n3_pos * bary_pos.dz();
        const vec3 cpa_f_to_node = n.pos() - cpa_pos;

        //Only compile this section if the faces store their contact energies
        #if FACE_STORE_CONTACT_ENERGY
            double adhesion_energy = 0.,    repulsion_energy = 0.;
        #endif

        //Get a reference to the types of the face
        const auto& face_type = c2->get_face_type(f->local_face_id_);

        //Determine if the contact is of adhesive or repulsive type based on the normal of the face
        bool adhesive_contact = (cpa_f_to_node.dot(f->normal_) > 0.0);

        //Reverse the direction of contact if the cell to which the face belongs is an ECM cell
        if (c1->get_cell_type_id() == 0 && c2->get_cell_type_id() == 1){adhesive_contact = !adhesive_contact;}

        //Do the same of the node belongs to a nucleus and the face to an epithelial cell
        if (c1->get_cell_type_id() == 3 && c2->get_cell_type_id() == 0){adhesive_contact = !adhesive_contact;}


        //Check that the faces have normals pointing in opposite directions
        if(adhesive_contact && min_squared_distance < interaction_cutoff_square_adhesion_){

            #if POLARIZATION_MODE_INDEX == 1
                //Indicate to the face that it is in contact with another cell
                c2->face_is_in_contact(f->local_face_id_, c1);
            #endif

            //Compute the distance
            double min_distance = std::sqrt(min_squared_distance);
            const double integration_region = f->area_;

            //Might change based on the hardening or softening regime
            double force_amplitude;

            //Hardening regime
            if(min_distance >= hardening_distance_){
                force_amplitude = face_type.adherence_strength_ * (interaction_cutoff_adhesion_ /min_distance - 1) * integration_region;

                #if FACE_STORE_CONTACT_ENERGY
                    adhesion_energy = integration_region * 0.5 * face_type.adherence_strength_ * (0.5 * interaction_cutoff_adhesion_ * interaction_cutoff_adhesion_ - pow(interaction_cutoff_adhesion_ - min_distance,2));
                #endif
            }

            //Softening regime
            else {
                force_amplitude =  face_type.adherence_strength_ * integration_region;

                #if FACE_STORE_CONTACT_ENERGY
                    adhesion_energy = 0.5 * face_type.adherence_strength_ * integration_region * min_distance * min_distance;
                #endif

            }

            #if FACE_STORE_CONTACT_ENERGY
                f->add_adhesion_energy(adhesion_energy);
            #endif

            //Apply the adhesion forces 
            const vec3 adhesion_force_vector =  cpa_f_to_node *  force_amplitude;

            //Linearly distribute the adhesion force onto the 3 nodes of the face
            f_n1.add_force(adhesion_force_vector * bary_pos.dx());
            f_n2.add_force(adhesion_force_vector * bary_pos.dy());
            f_n3.add_force(adhesion_force_vector * bary_pos.dz());

            //Apply the opposite force on the node
            n.add_force(adhesion_force_vector * -1.0);
        }
        


        //Repulsion regime
        if(!adhesive_contact && min_squared_distance < interaction_cutoff_square_repulsion_){

            #if POLARIZATION_MODE_INDEX == 1
            //Indicate to the face that it is in contact with another cell
            c2->face_is_in_contact(f->local_face_id_, c1);
            #endif

            //Compute the distance
            const double min_distance = std::sqrt(min_squared_distance);
            const double integration_region = f->area_;

            //Get a reference to the types of the face
            const auto& face_type = c2->get_face_type(f->local_face_id_);

            #if FACE_STORE_CONTACT_ENERGY
                repulsion_energy = 0.5 * face_type.repulsion_strength_ * integration_region * min_distance * min_distance;
                f->add_repulsion_energy(repulsion_energy);
            #endif

            //The force that will be distributed on the nodes of the 2 faces
            const vec3 repulsion_force_vector =  cpa_f_to_node  *  face_type.repulsion_strength_ * integration_region;

            //Apply the repulsion forces on the nodes of the 2 faces
            f_n1.add_force(repulsion_force_vector * bary_pos.dx());
            f_n2.add_force(repulsion_force_vector * bary_pos.dy());
            f_n3.add_force(repulsion_force_vector * bary_pos.dz());

            n.add_force(repulsion_force_vector * -1.0);   
        }
    }
}
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------




#endif

