#ifndef DEF_LUMEN_CELL
#define DEF_LUMEN_CELL

#include <cassert>

#include "../cell.hpp"


/*
    Use polymorphism to create the different cell types. Even though using IF statements is more efficient,
    polymorphism creates code that is more readable and maintainable. The IF statements are also more error prone.
*/



class lumen_cell: public cell {

    private: 
        friend class cell_divider;

    protected:


    public:

        //Make sure the lumen_cell cannot be default instantiated
        lumen_cell() = delete;                        //default constructor
        lumen_cell(const lumen_cell& c) = default;           //copy constructor
        lumen_cell(lumen_cell&& c) = delete;                //move constructor
        lumen_cell& operator=(const lumen_cell& c) = delete;//copy assignment operator
        lumen_cell& operator=(lumen_cell&& c) = delete;     //move assignment operator 


        //Construct the lumen_cell with a node_lst and a face_lst
        lumen_cell(
            const std::vector<node>& node_lst,
            const std::vector<face>& face_lst,
            unsigned cell_id,
            cell_type_param_ptr cell_type_parameters
        ) noexcept : cell(node_lst, face_lst, cell_id, cell_type_parameters){}


        //Construct the lumen_cell based on the node position and the node  ids of the faces
        lumen_cell(
            const std::vector<double>&  node_position,
            const std::vector<unsigned>& face_node_ids,
            unsigned cell_id,
            cell_type_param_ptr cell_type_parameters
        ) noexcept : cell(node_position, face_node_ids, cell_id, cell_type_parameters){}


        //Convert a mesh object into a lumen_cell object
        lumen_cell(
            const mesh& m,
            unsigned cell_id,
            cell_type_param_ptr cell_type_parameters
        ) noexcept : cell(m, cell_id, cell_type_parameters){}


        cell_ptr get_cell_same_type(const mesh& m) noexcept(false) override{
            //Return a copy of the cell
            return std::make_shared<lumen_cell>(m, cell_id_, cell_type_);
        }
};

#endif