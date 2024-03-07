#include "face.hpp"


/*
    Class for a triangular face of the cell mesh
*/


//---------------------------------------------------------------------------------------------------------
//Instantiate a face with the unique id of itss nodes
face::face(const unsigned node_id_1, const unsigned node_id_2, const unsigned node_id_3, const unsigned face_id) noexcept: local_face_id_(face_id){
    n1_id_ = node_id_1;
    n2_id_ = node_id_2;
    n3_id_ = node_id_3;
}
//---------------------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------------------

//Instantiate a face with the reference of its nodes 
face::face(const node& n1, const node& n2, const node& n3, const unsigned face_id) noexcept: local_face_id_(face_id){

    n1_id_ = n1.get_local_id();
    n2_id_ = n2.get_local_id();
    n3_id_ = n3.get_local_id();
}

//---------------------------------------------------------------------------------------------------------



//---------------------------------------------------------------------------------------------------------
//Returns the nodes of the face in an array
std::array<unsigned, 3> face::get_node_ids() const noexcept{

    //First check that the face is associated with nodes
    assert(is_used_);

    //Return the ids of the nodes
    return {n1_id_,  n2_id_,  n3_id_};
} 
//---------------------------------------------------------------------------------------------------------


//---------------------------------------------------------------------------------------------------------
//Replace one node of the face by another one
void face::replace_node(const unsigned old_node_id, const unsigned new_node_id) noexcept{
    //First check that the face is associated with nodes
    assert(is_used_);
    assert(old_node_id != new_node_id);
    assert(old_node_id == n1_id_ || old_node_id == n2_id_ || old_node_id == n3_id_);

    //Replace the old node id with the new one
    if      (n1_id_ == old_node_id) n1_id_ = new_node_id;
    else if (n2_id_ == old_node_id) n2_id_ = new_node_id;
    else                            n3_id_ = new_node_id;
}
//---------------------------------------------------------------------------------------------------------


//---------------------------------------------------------------------------------------------------------
//Only the cell can set the area of a face
void face::set_area(double area) noexcept{

    //Some checks before we accept the value
    assert(is_used_);
    assert(std::isfinite(area));
    assert(area >= 0.);

    area_ = area;
}
//---------------------------------------------------------------------------------------------------------


//---------------------------------------------------------------------------------------------------------
//Empties the face content
void face::reset() noexcept{

    //First check that the face is associated with nodes
    assert(is_used_);

    //The face is no longer used by the cell to store nodes
    is_used_ = false;

    //The face remains however connected to the cell

}
//---------------------------------------------------------------------------------------------------------


//---------------------------------------------------------------------------------------------------------
//The ids of the nodes might change, so we need to update them in each face
void face::update_node_ids(std::map<unsigned, unsigned>& node_id_correspondence) noexcept{

    //First check that the face is associated with nodes
    assert(is_used_);

    n1_id_ = node_id_correspondence[n1_id_];
    n2_id_ = node_id_correspondence[n2_id_];
    n3_id_ = node_id_correspondence[n3_id_];
}

//---------------------------------------------------------------------------------------------------------


//---------------------------------------------------------------------------------------------------------
unsigned face::get_opposite_node(const unsigned n1, const unsigned n2) const noexcept{

    //Make sure the face is used by nodes
    assert(is_used_);

    //Make sure the given nodes are in the face
    assert(n1 == n1_id_ || n1 == n2_id_ || n1 == n3_id_);
    assert(n2 == n1_id_ || n2 == n2_id_ || n2 == n3_id_);

    unsigned opposite_node_id;

    for(const auto node_id: get_node_ids()){
        if(node_id != n1 && node_id != n2){
            opposite_node_id = node_id;
            break;
        }
    }
    return opposite_node_id;
}
//---------------------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------------------
void face::swap_nodes() noexcept{
    assert( is_used_);
    std::swap(n2_id_, n3_id_);
}
//---------------------------------------------------------------------------------------------------------


//---------------------------------------------------------------------------------------------------------
//Check if they have the same adress in memory
bool face::operator==(const face& f) const noexcept{return local_face_id_ == f.get_local_id();}
bool face::operator!=(const face& f) const noexcept{return local_face_id_ != f.get_local_id();}
//---------------------------------------------------------------------------------------------------------
