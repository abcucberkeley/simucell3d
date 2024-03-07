#ifndef DEF_EDGE
#define DEF_EDGE

#include <cassert>
#include <iostream>
#include <optional>
#include <memory>
#include <unordered_set>


class node;
class face;


#include "node.hpp"
#include "face.hpp"
#include "custom_exception.hpp"



/*
    Stores the edge information of the cell mesh. An edge is defined by 2 nodes and 2 faces.
*/

class edge 
{

    private: 

        //The ids of the 2 nodes owning this edge
        unsigned n1_id_;
        unsigned n2_id_;

        //The ids of the 2 edges owning this edge
        std::optional<unsigned> f1_id_ = std::nullopt;
        std::optional<unsigned> f2_id_ = std::nullopt;


    public:

        //Make sure a edge cannot be instantiated without its nodes
        edge() = default; // Have to exist in order to create empty arrayrs of edges
        edge(const edge& e) = default;           //copy constructor
        edge(edge&& e) = default;                //move constructor
        edge& operator=(const edge& e) = default;//copy assignment operator
        edge& operator=(edge&& e) = default;     //move assignment operator 

        //Instantiate the edge with the reference of its nodes
        edge(const node& n1, const node& n2) noexcept;

        //Instantiate an edge just with its nodes ids
        explicit edge(unsigned n1, unsigned n2) noexcept;

        //Instantiate an edge with the ids of its nodes and faces
        explicit edge(unsigned n1, unsigned n2, unsigned f1, unsigned f2) noexcept;

        //Returns the ids of the nodes
        std::pair<unsigned, unsigned> get_node_ids() const noexcept;
        std::pair<unsigned, unsigned> get_face_ids() const noexcept;

        //Swap the ids of the 2 faces
        void swap_face_ids() noexcept;

        //Remove the face
        void delete_face(const unsigned f_id) noexcept(false);

        unsigned f1() const noexcept {assert(f1_id_); return f1_id_.value();}
        unsigned f2() const noexcept {assert(f2_id_); return f2_id_.value();}
        unsigned n1() const noexcept {return n1_id_;}
        unsigned n2() const noexcept {return n2_id_;}

        //Add a face to the edge
        void add_face(const unsigned f_id) noexcept(false);
        
        //Returns wether or not the edge is associated with 2 faces
        bool is_manifold() const noexcept;

        //Create a hash function to uniquely identify each edge object
        size_t hash() const noexcept{
            return 0.5 * (n1_id_ + n2_id_) * (n1_id_ + n2_id_ + 1) + n2_id_;
        }

        //Overload the equality operator
        bool operator==(const edge& e) const noexcept;
        bool operator!=(const edge& e) const noexcept;

        //Returns wether or not the edge has a face with the given id
        bool has_face(const unsigned f_id) const noexcept{
            return (f1_id_ && f1_id_.value() == f_id) || (f2_id_ && f2_id_.value() == f_id);
        }
        
        bool has_node(const unsigned n_id) const noexcept{
            return n1_id_ == n_id || n2_id_ == n_id;
        }

        //Replace a face by another one, doesn't affect the hash of the edge
        //Can only be used if the edge is manifold
        void replace_face(const unsigned old_f_id, const unsigned new_f_id) noexcept(false){
            assert(f1_id_ && f2_id_);
            assert(old_f_id == f1_id_ || old_f_id == f2_id_);
            if(f1_id_ == old_f_id){f1_id_ = new_f_id;}
            else{f2_id_ = new_f_id;}
        }

        bool operator<(const edge& e2) const{
            return hash() < e2.hash();
        }

};







//-----------------------------------------------------------------------------------------------------
//Create a hash out of an edge
struct edge_hasher{
    int operator()(const edge& e) const{
        
        //Create a hash function for pairs of edges with the Cantor pairing function
        //The created hash will always be unique
        const auto [e_n1, e_n2] = e.get_node_ids();
        return 0.5 * (e_n1 + e_n2) * (e_n1 + e_n2 + 1) + e_n2;
    }
};
//-----------------------------------------------------------------------------------------------------


//-----------------------------------------------------------------------------------------------------
//Compare pairs of edges
struct edge_comparator{
    bool operator()(const edge& e1,const edge& e2) const{
        return e1 == e2;
    }
};
//----------------------------------------------------------------------------------------------------



#endif