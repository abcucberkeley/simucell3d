#ifndef DEF_STATIC_CELL
#define DEF_STATIC_CELL

#include <cassert>

#include "../cell.hpp"


/*
    Use polymorphism to create the different cell types. Even though using IF statements is more efficient,
    polymorphism creates code that is more readable and maintainable. The IF statements are also more error prone.
*/


class static_cell: public cell {

    private: 
        friend class cell_divider;

    protected:


    public:

        //Make sure the static_cell cannot be default instantiated
        static_cell() = delete;                        //default constructor
        static_cell(const static_cell& c) = default;           //copy constructor
        static_cell(static_cell&& c) = delete;                //move constructor
        static_cell& operator=(const static_cell& c) = delete;//copy assignment operator
        static_cell& operator=(static_cell&& c) = delete;     //move assignment operator 


        
        //Construct the static_cell with a node_lst and a face_lst
        static_cell(
            const std::vector<node>& node_lst,
            const std::vector<face>& face_lst,
            unsigned cell_id,
            cell_type_param_ptr cell_type_parameters
        ) noexcept : cell(node_lst, face_lst, cell_id, cell_type_parameters){

            //This cell should not move
            is_static_ = true;
        }

        //Construct the static_cell based on the node position and the node  ids of the faces
        static_cell(
            const std::vector<double>&  node_position,
            const std::vector<unsigned>& face_node_ids,
            unsigned cell_id,
            cell_type_param_ptr cell_type_parameters
        ) noexcept : cell(node_position, face_node_ids, cell_id, cell_type_parameters){
                
                //This cell should not move
                is_static_ = true;
        }


        //Convert a mesh object into a static_cell object
        static_cell(
            const mesh& m,
            unsigned cell_id,
            cell_type_param_ptr cell_type_parameters
        ) noexcept : cell(m, cell_id, cell_type_parameters){
                
                //This cell should not move
                is_static_ = true;
        }


        cell_ptr get_cell_same_type(const mesh& m) noexcept(false) override{
            //Return a copy of the cell
            return std::make_shared<static_cell>(m, cell_id_, cell_type_);
        }
};

#endif