#include "node.hpp"



//Trivial constructor
//---------------------------------------------------------------------------------------------------------
node::node(const double dx, const double dy, const double dz, const unsigned node_id) noexcept :node_id_(node_id){

    //Check the coordinates of the node
    assert(std::isfinite(dx) && std::isfinite(dy) && std::isfinite(dz));

    //Set the position of the node
    pos_.translate(dx, dy, dz);


    //Overdamped equations of motion solved with the improved euler scheme
    #if DYNAMIC_MODEL_INDEX == 2

        //The previous position is the starting position
        previous_pos_.translate(dx, dy, dz);
    #endif 

    #if CONTACT_MODEL_INDEX == 1 || CONTACT_MODEL_INDEX == 2
        omp_init_lock(&lock_); // Initialize the lock in the constructor
    #endif

}

node::node(const vec3& pos, const unsigned node_id) noexcept{
    node_id_ = node_id;
    pos_ = pos;

    //Overdamped equations of motion solved with the improved euler scheme
    #if DYNAMIC_MODEL_INDEX == 2
     
        //The previous position is the starting position
        previous_pos_ = pos;
    #endif 

    #if CONTACT_MODEL_INDEX == 1 || CONTACT_MODEL_INDEX == 2
        omp_init_lock(&lock_); // Initialize the lock in the constructor
    #endif
}

//---------------------------------------------------------------------------------------------------------
//Reset the vectors owned by the node
void node::reset() noexcept {
    assert(is_used_ = true);

    pos_.reset();
    force_.reset();


    #if CONTACT_MODEL_INDEX == 1
        //Remove the coupling to the face
        coupled_node_.reset();

    #endif



    #if CONTACT_MODEL_INDEX == 2
        //Remove the coupling to the face
        coupled_nodes_map_.clear();

    #endif

    //Overdamped equations of motion solved with the improved euler scheme
    #if DYNAMIC_MODEL_INDEX == 0
        momentum_.reset();
    
    #elif DYNAMIC_MODEL_INDEX == 2
    
        previous_pos_.reset();
        previous_force_.reset();
    #endif 

    //Indicate that the node is not used by the cell anymore
    is_used_ = false;
}
//---------------------------------------------------------------------------------------------------------



//---------------------------------------------------------------------------------------------------------
//Check if the 2 nodes have the same id
bool node::operator==(const node& n) const noexcept {return node_id_ == n.get_local_id();}

bool node::operator!=(const node& n) const noexcept{return node_id_ != n.get_local_id();}
//---------------------------------------------------------------------------------------------------------

//The substraction of 2 nodes returns a vector
vec3 node::operator-(const node& n) const noexcept{return pos_ - n.pos();}

