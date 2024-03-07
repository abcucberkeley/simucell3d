#include "time_integration.hpp"

/*
SimuCell3D possesses three methods to solve the equations of motions:
    - DYNAMIC_MODEL_INDEX = 0 : Full equations of motions solved with semi-implicit euler  (1st order integration scheme)
    - DYNAMIC_MODEL_INDEX = 1 : Overdamped equations of motions solved with forward euler  (1st order integration scheme)
    - DYNAMIC_MODEL_INDEX = 2 : Overdamped equations of motions solved with improved euler (2nd order integration scheme)
//The DYNAMIC_MODEL_INDEX can be changed in the ./include/global_configuration.hpp file.

For more information, please refer to the paper associated to SimuCell3D.
*/ 



void time_integration_scheme::update_nodes_positions(const std::vector<cell_ptr>& cell_lst){

    //Regardless of the scheme used to update the node positions, we first need to reset the cell kinetic energies
    std::for_each(cell_lst.begin(), cell_lst.end(), [](cell_ptr c){c->kinetic_energy_ = 0.0;});

    //----------------------------------------------------------------------------------------------------------------------
    //If we use the contact model that uses springs between faces and vertices
    #if CONTACT_MODEL_INDEX == 0

        #pragma omp parallel for schedule(static)
        for(size_t c1_id = 0; c1_id < cell_lst.size(); c1_id++){
            cell_ptr c1 = cell_lst[c1_id];

            if(c1->is_static_){continue;}

            const double c1_node_mass = c1->get_node_mass();

            for(node& n1: c1->node_lst_){

                if(n1.is_used()){

                    //If we use the dynamic model that solves the full equations of motion
                    #if DYNAMIC_MODEL_INDEX == 0

                        //First update the momentum of the node
                        n1.momentum_.translate((n1.force_ - (n1.momentum_ * (damping_coeff_ / c1_node_mass))) * dt_);

                        //Compute the node velocity on the fly and update the node position
                        n1.pos_.translate(n1.momentum_ * (dt_ / c1_node_mass)); 

                        //Update the cell kinetic energy
                        c1->kinetic_energy_ += 0.5 * n1.momentum_.squared_norm() / c1_node_mass;

                    //If we use the dynamic model that solves the overdamped equations of motion
                    #elif DYNAMIC_MODEL_INDEX == 1

                        n1.pos_.translate(n1.force_ * (dt_ / damping_coeff_));

                        //Update the cell kinetic energy
                        c1->kinetic_energy_ += 0.5 * (n1.force_ * (c1_node_mass / damping_coeff_)).squared_norm() ;

                        //Reset the force of the node to the 0 vector
                        n1.force_.reset();

                    #else
                        std::runtime_error("The combination of DYNAMIC_MODEL_INDEX = " + std::to_string(DYNAMIC_MODEL_INDEX) + " and CONTACT_MODEL_INDEX = " + std::to_string(CONTACT_MODEL_INDEX) + " is not supported by the time integration scheme. Please choose another combination.");
                    #endif

                    //Reset the force of the node to the 0 vector
                    n1.force_.reset();
                }
            }
        }
    //----------------------------------------------------------------------------------------------------------------------


    //----------------------------------------------------------------------------------------------------------------------
    //If we use the contact model that mechanically couples nodes of adjacent cells
    #elif CONTACT_MODEL_INDEX == 1

    #pragma omp parallel for schedule(static)
    for(size_t c1_id = 0; c1_id < cell_lst.size(); c1_id++){
        cell_ptr c1 = cell_lst[c1_id];
        if(c1->is_static_){continue;}

        const double c1_node_mass = c1->get_node_mass();

        for(node& n1: c1->node_lst_){

            if(n1.is_used()){

                //If the node is coupled to another node
                if(n1.coupled_node_.has_value()){

                    //Get the coupled node
                    const auto [c2_id, n2_id] = n1.coupled_node_.value();

                    if(c1->get_local_id() > c2_id){
                        assert(c2_id < cell_lst.size());

                        //Get the coupled node
                        cell_ptr c2 = cell_lst[c2_id];
                        const double c2_node_mass = c2->get_node_mass();
                        
                        assert(n2_id < c2->node_lst_.size());
                        node& n2 = c2->node_lst_[n2_id];

                        //If we use the dynamic model that solves the full equations of motion
                        #if DYNAMIC_MODEL_INDEX == 0
                            //Sum the momentum and forces of the 2 nodes
                            const vec3 avg_momentum = (n1.momentum_ + n2.momentum_) * (0.5);
                            n1.momentum_.reset(avg_momentum);
                            n2.momentum_.reset(avg_momentum);
                        #endif


                        const vec3 avg_force = (n1.force_ + n2.force_) * (0.5);
                        n1.force_.reset(avg_force);
                        n2.force_.reset(avg_force);

                        //Sum the mass of the 2 nodes
                        const double avg_node_mass = (c1_node_mass + c2_node_mass) *0.5;

                        //If we use the dynamic model that solves the full equations of motion
                        #if DYNAMIC_MODEL_INDEX == 0

                            //Update simulataaneously the momentum of the 2 nodes
                            n1.momentum_.translate((n1.force_ - (n1.momentum_ * (damping_coeff_ / avg_node_mass))) * dt_);
                            n2.momentum_.translate((n2.force_ - (n2.momentum_ * (damping_coeff_ / avg_node_mass))) * dt_);

                            //Do the same for the position of the 2 nodes
                            n1.pos_.translate(n1.momentum_ * (dt_ / avg_node_mass));
                            n2.pos_.translate(n2.momentum_ * (dt_ / avg_node_mass));

                            //Update the kinetic energy of the 2 cells
                            const double kinetic_energy_node = 0.5 * n1.momentum_.squared_norm() / avg_node_mass;
                            c1->kinetic_energy_ += kinetic_energy_node;
                            c2->kinetic_energy_ += kinetic_energy_node;


                        //If we use the dynamic model that solves the overdamped equations of motion
                        #elif DYNAMIC_MODEL_INDEX == 1

                            //Move the nodes
                            n1.pos_.translate(avg_force * (dt_ / damping_coeff_));
                            n2.pos_.translate(avg_force * (dt_ / damping_coeff_));

                            //Compute the kinetic energy of the nodes
                            const double node_kinetic_energy = 0.5 * (avg_force * (dt_ / damping_coeff_)).squared_norm() / avg_node_mass;
                            c1->kinetic_energy_ += node_kinetic_energy;
                            c2->kinetic_energy_ += node_kinetic_energy;

                        #else
                            std::runtime_error("The combination of DYNAMIC_MODEL_INDEX = " + std::to_string(DYNAMIC_MODEL_INDEX) + " and CONTACT_MODEL_INDEX = " + std::to_string(CONTACT_MODEL_INDEX) + " is not supported by the time integration scheme. Please choose another combination.");
                        #endif


                        //Reset the forces of the nodes
                        n1.force_.reset();
                        n2.force_.reset();
                    }
                }

          
                //If the node is not coupled
                else{
                    
                    //If we use the dynamic model that solves the full equations of motion
                    #if DYNAMIC_MODEL_INDEX == 0

                        //First update the momentum of the node
                        n1.momentum_.translate((n1.force_ - (n1.momentum_ * (damping_coeff_ / c1_node_mass))) * dt_);

                        //Compute the node velocity on the fly and update the node position
                        n1.pos_.translate(n1.momentum_ * (dt_ / c1_node_mass)); 

                        //Update the cell kinetic energy
                        c1->kinetic_energy_ += 0.5 * n1.momentum_.squared_norm() / c1_node_mass;

                    //If we use the dynamic model that solves the overdamped equations of motion
                    #elif DYNAMIC_MODEL_INDEX == 1
                            
                        n1.pos_.translate(n1.force_ * (dt_ / damping_coeff_));

                        //Update the cell kinetic energy
                        c1->kinetic_energy_ += 0.5 * (n1.force_ * (c1_node_mass / damping_coeff_)).squared_norm() ;

                        //Reset the force of the node to the 0 vector
                        n1.force_.reset();

                    #else

                        std::runtime_error("The combination of DYNAMIC_MODEL_INDEX = " + std::to_string(DYNAMIC_MODEL_INDEX) + " and CONTACT_MODEL_INDEX = " + std::to_string(CONTACT_MODEL_INDEX) + " is not supported by the time integration scheme. Please choose another combination.");

                    #endif

                    //Reset the force of the node to the 0 vector
                    n1.force_.reset();
                }
            }
        }
    }
    //----------------------------------------------------------------------------------------------------------------------




    //----------------------------------------------------------------------------------------------------------------------
    //If we use the contact model that mechanically couples faces of adjacent cells
    #elif CONTACT_MODEL_INDEX == 2

        for(size_t c1_id = 0; c1_id < cell_lst.size(); c1_id++){
            cell_ptr c1 = cell_lst[c1_id];
            if(c1->is_static_){continue;}

            for(node& n1: c1->node_lst_){

                if(n1.is_used() == false){continue;}
                
                //If this node is coupled to other nodes on other cells, and this cell has the greatest id
                bool c1_has_greatest_id = std::all_of(n1.coupled_nodes_map_.begin(), n1.coupled_nodes_map_.end(), 
                    [&](const auto& coupled_node_data){
                        return c1->get_local_id() > coupled_node_data.first;
                    }
                );

                //If this node is coupled to other nodes on other cells, and one of them has a greater id than this cell
                if(c1_has_greatest_id == false){continue;}

                //Loop over the coupled nodes and calculate the average force that should be applied to the nodes
                vec3 avg_force = n1.force_;
                double avg_node_mass = c1->get_node_mass();

                //If we use the dynamic model that solves the full equations of motion
                #if DYNAMIC_MODEL_INDEX == 0
                    vec3 avg_momentum = n1.momentum_;
                #endif


                for(auto it = n1.coupled_nodes_map_.begin() ; it != n1.coupled_nodes_map_.end(); it++){

                    const unsigned c2_local_id = it->first;
                    const unsigned n2_local_id = it->second.first;

                    //Get the adjacent cell
                    assert(c2_local_id < cell_lst.size());
                    cell_ptr c2 = cell_lst[c2_local_id];

                    //Get the coupled node in the adjacent cell
                    assert(n2_local_id < c2->node_lst_.size());
                    node& n2 = c2->node_lst_[n2_local_id];

                    //Add the force of the second node to the average force
                    avg_force = avg_force + n2.force_;
                    avg_node_mass = avg_node_mass + c2->get_node_mass();
                    
                    #if DYNAMIC_MODEL_INDEX == 0
                        avg_momentum = avg_momentum + n2.momentum_;
                    #endif 
                }

                //Divide the average force by the number of coupled nodes
                avg_force = avg_force  / (n1.get_nb_coupled_nodes() + 1);
                avg_node_mass = avg_node_mass  / (n1.get_nb_coupled_nodes() + 1);


                //If we use the dynamic model that solves the full equations of motion
                #if DYNAMIC_MODEL_INDEX == 0
                    avg_momentum = avg_momentum  / (n1.get_nb_coupled_nodes() + 1);
                    const double kinetic_energy_node = 0.5 * avg_momentum.squared_norm() / avg_node_mass;
                
                    //Update the position and momentum of the node n1
                    n1.momentum_.translate((avg_force - (avg_momentum * (damping_coeff_ / avg_node_mass))) * dt_);
                    n1.pos_.translate(avg_momentum* (dt_ / avg_node_mass)); 

                #elif DYNAMIC_MODEL_INDEX == 1

                    const double kinetic_energy_node = 0.5 * (avg_force * (avg_node_mass / damping_coeff_)).squared_norm() ;
                    
                    //Update the position of the node n1
                    n1.pos_.translate(avg_force * (dt_ / damping_coeff_));
                    

                #else
                    std::runtime_error("The combination of DYNAMIC_MODEL_INDEX = " + std::to_string(DYNAMIC_MODEL_INDEX) + " and CONTACT_MODEL_INDEX = " + std::to_string(CONTACT_MODEL_INDEX) + " is not supported by the time integration scheme. Please choose another combination.");
                #endif 

                c1->kinetic_energy_ += kinetic_energy_node;
                n1.force_.reset();


                //Apply the average force to all the nodes coupled to this node
                for(auto it = n1.coupled_nodes_map_.begin() ; it != n1.coupled_nodes_map_.end(); it++){

                    const unsigned c2_local_id = it->first;
                    const unsigned n2_local_id = it->second.first;

                    //Get the adjacent cell
                    assert(c2_local_id < cell_lst.size());
                    cell_ptr c2 = cell_lst[c2_local_id];

                    //Get the coupled node in the adjacent cell
                    assert(n2_local_id < c2->node_lst_.size());
                    node& n2 = c2->node_lst_[n2_local_id];
                    
                    #if DYNAMIC_MODEL_INDEX == 0
                        n2.momentum_.translate((avg_force- (avg_momentum * (damping_coeff_ / avg_node_mass))) * dt_);
                        n2.pos_.translate(avg_momentum * (dt_ / avg_node_mass)); 

                    #elif DYNAMIC_MODEL_INDEX == 1
                        n2.pos_.translate(avg_force * (dt_ / damping_coeff_));
                    #endif

                    c2->kinetic_energy_ += kinetic_energy_node;
                    n2.force_.reset();
                }
            }
        }

    #else
        std::runtime_error("The contact model corresponding to the INDEX " + std::to_string(CONTACT_MODEL_INDEX) + " is not supported by the time integration scheme. Please choose another contact model.");
    #endif


    //Update the simulation time
    simulation_time_ += dt_;
    
}










