#include "edge.hpp"


/*
Stores the edge information of the cell mesh. An edge is defined by 2 nodes and 2 faces.
*/

//----------------------------------------------------------------------------
//Constructors
edge::edge(const node& n1, const node& n2) noexcept{
    //assert(n1.get_local_id() != n2.get_local_id());

    //Make sure the node 1 is the one with the smallest id
    if(n1.get_local_id() < n2.get_local_id()) {
        n1_id_ = n1.get_local_id();
        n2_id_ = n2.get_local_id();
    }else{
        n1_id_ = n2.get_local_id();
        n2_id_ = n1.get_local_id();
    }

}

//Instantiate an edge just with its nodes ids
edge::edge(unsigned n1, unsigned n2) noexcept{

    //assert(n1 != n2);
    //Make sure the node 1 has the smallest id
    if(n1 < n2) {
        n1_id_ = n1;
        n2_id_ = n2;
    }else{
        n1_id_ = n2;
        n2_id_ = n1;
    }

}

//Instantiate an edge with the ids of its nodes and faces
edge::edge(unsigned n1, unsigned n2, unsigned f1, unsigned f2) noexcept: f1_id_(f1), f2_id_(f2){
    
    //assert(n1 != n2); 
    assert(f1 != f2);

    //Make sure the node 1 has the smallest id
    if(n1 < n2) {
        n1_id_ = n1;
        n2_id_ = n2;
    }else{
        n1_id_ = n2;
        n2_id_ = n1;
    }

}
//----------------------------------------------------------------------------




//----------------------------------------------------------------------------
//Add a new face to the edge
void edge::add_face(const unsigned f_id) noexcept(false){

    //Add the face id
    if(!f1_id_) {
        f1_id_ = f_id;
    }
    else if(!f2_id_){
        //Check that the same face is not added twice 
        assert(f1_id_.value() != f_id);
        f2_id_ = f_id;
    }
    //This means that the edge is connected to more than 2 faces
    else{
        throw mesh_integrity_exception("Edge ("+std::to_string(n1_id_)+"--"+std::to_string(n2_id_)+") already connected to the faces "+std::to_string(f1_id_.value())+" and "+std::to_string(f2_id_.value())+". Cannot add the face "+std::to_string(f_id)+".");
    }
}
//----------------------------------------------------------------------------


//----------------------------------------------------------------------------
//Swap the ids of the 2 faces
void edge::swap_face_ids() noexcept{
    assert(f1_id_ && f2_id_);
    std::swap(f1_id_, f2_id_);
}
//----------------------------------------------------------------------------



//----------------------------------------------------------------------------
//Returns wether or not the edge is associated with 2 faces
bool edge::is_manifold() const noexcept{
    return f1_id_.has_value() && f2_id_.has_value();

}
//----------------------------------------------------------------------------

//----------------------------------------------------------------------------
//Remove the face
void edge::delete_face(const unsigned f_id) noexcept(false){
    assert(f1_id_ || f2_id_);

    if(f1_id_ && f1_id_.value() == f_id) {

        //Delete the first face and then move the 2nd face to the first position
        if (f2_id_) swap_face_ids();
        f2_id_.reset();
    }
    else if(f2_id_ && f2_id_.value() == f_id){
        f2_id_.reset();
    }
    else{
        throw mesh_integrity_exception("Face not found in edge");
    }
}
//----------------------------------------------------------------------------



//----------------------------------------------------------------------------
//Returns the ids of the nodes composing this edge
std::pair<unsigned, unsigned> edge::get_node_ids() const noexcept{

    return {n1_id_, n2_id_};
}
//----------------------------------------------------------------------------



//----------------------------------------------------------------------------
//Returns the ids of the nodes composing this edge
std::pair<unsigned, unsigned> edge::get_face_ids() const noexcept{

    //Check that the nodes are defined in the edge
    assert(f1_id_ && f2_id_);

    return {f1_id_.value(), f2_id_.value()};
}

//----------------------------------------------------------------------------



//----------------------------------------------------------------------------

//Equal operator overload
bool edge::operator==(edge const& e2) const noexcept{
    const auto [e2_n1, e2_n2] = e2.get_node_ids();
    return n1_id_ == e2_n1 && n2_id_ == e2_n2;
}



bool edge::operator!=(edge const& e2) const noexcept{
    const auto [e2_n1, e2_n2] = e2.get_node_ids();
    return !(n1_id_ == e2_n1 && n2_id_ == e2_n2);
}
//----------------------------------------------------------------------------
