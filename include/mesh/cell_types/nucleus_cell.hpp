#ifndef DEF_NUCLEUS_CELL
#define DEF_NUCLEUS_CELL

#include <cassert>

#include "../cell.hpp"


/*
    Use polymorphism to create the different cell types. Even though using IF statements is more efficient,
    polymorphism creates code that is more readable and maintainable. The IF statements are also more error prone.
*/



class nucleus_cell: public cell {

    private: 
        friend class cell_divider;

    protected:


    public:

        //Make sure the nucleus_cell cannot be default instantiated
        nucleus_cell() = delete;                        //default constructor
        nucleus_cell(const nucleus_cell& c) = default;           //copy constructor
        nucleus_cell(nucleus_cell&& c) = delete;                //move constructor
        nucleus_cell& operator=(const nucleus_cell& c) = delete;//copy assignment operator
        nucleus_cell& operator=(nucleus_cell&& c) = delete;     //move assignment operator 


        //Construct the nucleus_cell with a node_lst and a face_lst
        nucleus_cell(
            const std::vector<node>& node_lst,
            const std::vector<face>& face_lst,
            unsigned cell_id,
            cell_type_param_ptr cell_type_parameters
        ) noexcept : cell(node_lst, face_lst, cell_id, cell_type_parameters){}


        //Construct the nucleus_cell based on the node position and the node  ids of the faces
        nucleus_cell(
            const std::vector<double>&  node_position,
            const std::vector<unsigned>& face_node_ids,
            unsigned cell_id,
            cell_type_param_ptr cell_type_parameters
        ) noexcept : cell(node_position, face_node_ids, cell_id, cell_type_parameters){}


        //Convert a mesh object into a nucleus_cell object
        nucleus_cell(
            const mesh& m,
            unsigned cell_id,
            cell_type_param_ptr cell_type_parameters
        ) noexcept : cell(m, cell_id, cell_type_parameters){}

        cell_ptr get_cell_same_type(const mesh& m) noexcept(false) override{
            //Return a copy of the cell
            return std::make_shared<nucleus_cell>(m, cell_id_, cell_type_);
        }


};

#endif