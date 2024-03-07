#ifndef DEF_ECM_CELL
#define DEF_ECM_CELL

#include <cassert>

#include "../cell.hpp"


/*
    Use polymorphism to create the different cell types. Even though using IF statements is more efficient,
    polymorphism creates code that is more readable and maintainable. The IF statements are also more error prone.
*/



class ecm_cell: public cell {

    private: 
        friend class cell_divider;

    protected:


    public:

        //Make sure the ecm_cell cannot be default instantiated
        ecm_cell() = delete;                        //default constructor
        ecm_cell(const ecm_cell& c) = default;           //copy constructor
        ecm_cell(ecm_cell&& c) = delete;                //move constructor
        ecm_cell& operator=(const ecm_cell& c) = delete;//copy assignment operator
        ecm_cell& operator=(ecm_cell&& c) = delete;     //move assignment operator 


        
        //Construct the ecm_cell with a node_lst and a face_lst
        ecm_cell(
            const std::vector<node>& node_lst,
            const std::vector<face>& face_lst,
            unsigned cell_id,
            cell_type_param_ptr cell_type_parameters
        ) noexcept : cell(node_lst, face_lst, cell_id, cell_type_parameters){
            
            //This cell should not move
            is_static_ = true;
        }


        //Construct the ecm_cell based on the node position and the node  ids of the faces
        ecm_cell(
            const std::vector<double>&  node_position,
            const std::vector<unsigned>& face_node_ids,
            unsigned cell_id,
            cell_type_param_ptr cell_type_parameters
        ) noexcept : cell(node_position, face_node_ids, cell_id, cell_type_parameters){

            //This cell should not move
            is_static_ = true;
        }


        //Convert a mesh object into a ecm_cell object
        ecm_cell(
            const mesh& m,
            unsigned cell_id,
            cell_type_param_ptr cell_type_parameters
        ) noexcept : cell(m, cell_id, cell_type_parameters){

            //This cell should not move
            is_static_ = true;
        }


        void apply_internal_forces(const double time_step) noexcept override{
            //Only if the cell is not static we apply the internal forces such as the cell pressure, sureface tension etc
            if(!is_static_){cell::apply_internal_forces(time_step);}
        }


        cell_ptr get_cell_same_type(const mesh& m) noexcept(false) override{
            //Return a copy of the cell
            return std::make_shared<ecm_cell>(m, cell_id_, cell_type_);
        }

};

#endif