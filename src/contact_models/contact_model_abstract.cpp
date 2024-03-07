#include "contact_model_abstract.hpp"










//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//The trivial constructor of the contact model
contact_model_abstract::contact_model_abstract(const global_simulation_parameters& sim_parameters) noexcept(false){

        //Store the cutoffs for the adhesion and repulsion
        interaction_cutoff_adhesion_ = sim_parameters.contact_cutoff_adhesion_;
        interaction_cutoff_square_adhesion_ = interaction_cutoff_adhesion_*interaction_cutoff_adhesion_;

        interaction_cutoff_repulsion_ = sim_parameters.contact_cutoff_repulsion_;
        interaction_cutoff_square_repulsion_ = interaction_cutoff_repulsion_*interaction_cutoff_repulsion_;
        max_interaction_cutoff_square_ = std::max(interaction_cutoff_square_repulsion_, interaction_cutoff_square_adhesion_);

        //The extra padding added to the AABB of the faces, to check for potential contacts
        aabb_padding_ = std::max(sim_parameters.contact_cutoff_repulsion_, sim_parameters.contact_cutoff_adhesion_);

        //The length size of each voxel used to instantiate the uniform space partiotionning grid
        grid_.voxel_size_ = sim_parameters.min_edge_len_ * 3.0 + 2. * aabb_padding_;

}
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------





//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//Check if the node is in the AABB of the face
bool contact_model_abstract::aabb_intersection_check(const size_t face_aabb_pos, const vec3& node_pos) const noexcept{
    
    if (node_pos.dx() < face_aabb_lst_[face_aabb_pos    ] || node_pos.dx() > face_aabb_lst_[face_aabb_pos + 3]) return false;
    if (node_pos.dy() < face_aabb_lst_[face_aabb_pos + 1] || node_pos.dy() > face_aabb_lst_[face_aabb_pos + 4]) return false;
    if (node_pos.dz() < face_aabb_lst_[face_aabb_pos + 2] || node_pos.dz() > face_aabb_lst_[face_aabb_pos + 5]) return false;
    return true;
}
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------





//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//Compute the axis aligned bounding box of ech face and store it in the face_aabb_lst_ vec3
void contact_model_abstract::update_face_aabbs(
    const std::vector<cell_ptr>& cell_lst) noexcept{

    //Get the number of faces in the mesh
    const size_t nb_faces = face_lst_.size(); 

    //Reset all the face AABBs
    face_aabb_lst_.clear();
    face_aabb_lst_.reserve(nb_faces * 6);

    //Reset the min and max coordinates of the mesh
    global_min_x_ = std::numeric_limits<double>::infinity();
    global_min_y_ = std::numeric_limits<double>::infinity();
    global_min_z_ = std::numeric_limits<double>::infinity();

    global_max_x_ = -std::numeric_limits<double>::infinity();
    global_max_y_ = -std::numeric_limits<double>::infinity();
    global_max_z_ = -std::numeric_limits<double>::infinity();

    //Loop over the faces in the mesh
    for(unsigned i = 0; i < nb_faces; i++){
        face* f =  face_lst_[i];
        assert(f != nullptr && f->is_used());

        //Get the cell that owns the face
        cell_ptr c = f->get_owner_cell();

        //Get the coordinates of the nodes composing the face
        const vec3& n1_pos = c->get_const_ref_node(f->n1_id_).pos();
        const vec3& n2_pos = c->get_const_ref_node(f->n2_id_).pos();
        const vec3& n3_pos = c->get_const_ref_node(f->n3_id_).pos();

        //Compute the min and max coordinates of the face
        const double face_min_x = std::min(n1_pos.dx(), std::min(n2_pos.dx(), n3_pos.dx())) - aabb_padding_;
        const double face_min_y = std::min(n1_pos.dy(), std::min(n2_pos.dy(), n3_pos.dy())) - aabb_padding_;
        const double face_min_z = std::min(n1_pos.dz(), std::min(n2_pos.dz(), n3_pos.dz())) - aabb_padding_;

        const double face_max_x = std::max(n1_pos.dx(), std::max(n2_pos.dx(), n3_pos.dx())) + aabb_padding_;
        const double face_max_y = std::max(n1_pos.dy(), std::max(n2_pos.dy(), n3_pos.dy())) + aabb_padding_;
        const double face_max_z = std::max(n1_pos.dz(), std::max(n2_pos.dz(), n3_pos.dz())) + aabb_padding_;
    
        //Update the global min and max coordinates of the mesh
        if(face_min_x < global_min_x_) global_min_x_ = face_min_x;
        if(face_min_y < global_min_y_) global_min_y_ = face_min_y;
        if(face_min_z < global_min_z_) global_min_z_ = face_min_z;

        if(face_max_x > global_max_x_) global_max_x_ = face_max_x;
        if(face_max_y > global_max_y_) global_max_y_ = face_max_y;
        if(face_max_z > global_max_z_) global_max_z_ = face_max_z;

        //Store the AABB of the face in the face_aabb_lst_ vec3        
        face_aabb_lst_.insert(face_aabb_lst_.end(), {face_min_x, face_min_y, face_min_z, face_max_x , face_max_y, face_max_z});
    } 
    
    //---------------------------------------------------------------------------------------------------------

    //Update the space partitioning grid dimensions
    global_min_x_ -= aabb_padding_;
    global_min_y_ -= aabb_padding_;
    global_min_z_ -= aabb_padding_;
}
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//Store all the faces in an uinform space partitioning grid
void contact_model_abstract::store_face_in_uspg() noexcept{

    //Empties the space partitioning grid and update its dimensions
    grid_.update_dimensions(face_lst_.size(), global_min_x_, global_min_y_, global_min_z_, global_max_x_, global_max_y_, global_max_z_);

    //Loop over all the faces in the mesh
    for(size_t i =  0; i < face_lst_.size(); i++){

        //Get the face
        face* f = face_lst_[i];
        assert(f != nullptr && f->is_used());

        //Get the position where the AABB of the face is stored in the face_aabb_lst_ vec3
        size_t face_aabb_pos = i * 6;

        //Get the coordinates of the min points of the face AABB
        const double face_min_x = face_aabb_lst_[face_aabb_pos    ]; 
        const double face_min_y = face_aabb_lst_[face_aabb_pos + 1];
        const double face_min_z = face_aabb_lst_[face_aabb_pos + 2];
        const double face_max_x = face_aabb_lst_[face_aabb_pos + 3];
        const double face_max_y = face_aabb_lst_[face_aabb_pos + 4];
        const double face_max_z = face_aabb_lst_[face_aabb_pos + 5];

        //Get in which voxels of the uniform space partitionning grid the face should be stored
        const unsigned voxel_x_start = static_cast<unsigned>(std::floor((face_min_x - grid_.min_x_) / grid_.voxel_size_));
        const unsigned voxel_y_start = static_cast<unsigned>(std::floor((face_min_y - grid_.min_y_) / grid_.voxel_size_));
        const unsigned voxel_z_start = static_cast<unsigned>(std::floor((face_min_z - grid_.min_z_) / grid_.voxel_size_));

        const unsigned voxel_x_stop  = static_cast<unsigned>(std::floor((face_max_x - grid_.min_x_) / grid_.voxel_size_));
        const unsigned voxel_y_stop  = static_cast<unsigned>(std::floor((face_max_y - grid_.min_y_) / grid_.voxel_size_));
        const unsigned voxel_z_stop  = static_cast<unsigned>(std::floor((face_max_z - grid_.min_z_) / grid_.voxel_size_));


        for(unsigned voxel_x = voxel_x_start; voxel_x <= voxel_x_stop; voxel_x++){
        for(unsigned voxel_y = voxel_y_start; voxel_y <= voxel_y_stop; voxel_y++){
        for(unsigned voxel_z = voxel_z_start; voxel_z <= voxel_z_stop; voxel_z++){

                //Insert the face in the space partitioning grid
                const size_t voxel_id = grid_.get_voxel_index(voxel_x, voxel_y, voxel_z);
                grid_.place_object(f, voxel_id);

        }
        }
        }
    }
}
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------





//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//Find the closest point on a triangle to a given point, returns the barycentric coordinates of the closest point as well as the closest point
// and the distance to the closest point
std::pair<double, vec3> contact_model_abstract::compute_node_triangle_distance(
    const vec3& p, 
    const vec3& a, 
    const vec3& b, 
    const vec3& c
) noexcept{


    // Check if P in vertex region outside A
    const vec3 ab = b - a;
    const vec3 ac = c - a;
    const vec3 ap = p - a;
    const double  d1 = ab.dot(ap);
    const double  d2 = ac.dot(ap);
    if (d1 <= 0.0 && d2 <= 0.0) return {ap.squared_norm(), vec3(1,0,0)}; // barycentric coordinates (1,0,0)


    // Check if P in vertex region outside B
    const vec3 bp = p - b;
    const double d3 = ab.dot(bp);
    const double d4 = ac.dot(bp);
    if (d3 >= 0.0 && d4 <= d3) return {bp.squared_norm(), vec3(0,1,0)}; // barycentric coordinates (0,1,0)



    // Check if P in edge region of AB, if so return projection of P onto AB
    const double vc = d1*d4 - d3*d2;
    if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) {
        const double v = d1 / (d1 - d3);
        const vec3 abv = a + ab * v;
        return {(abv - p).squared_norm(), vec3(1-v,v,0)}; // barycentric coordinates (1-v,v,0)
    }

    // Check if P in vertex region outside C
    const vec3 cp = p - c;
    const double d5 = ab.dot(cp);
    const double d6 = ac.dot(cp);
    if (d6 >= 0.0 && d5 <= d6){
        return {cp.squared_norm(), vec3(0,0,1)}; // barycentric coordinates (0,0,1)
    }

    // Check if P in edge region of AC, if so return projection of P onto AC
    const double vb = d5*d2 - d1*d6;
    if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) {
        const double w = d2 / (d2 - d6);
        const vec3 acw = a + ac * w;
        return {(p - acw).squared_norm(), vec3(1-w,0,w)}; // barycentric coordinates (1-w,0,w)
    }

    // Check if P in edge region of BC, if so return projection of P onto BC
    const double va = d3*d6 - d5*d4;
    if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0) {
        const double z = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        const vec3 bcz = b + (c - b) * z;
        return {(bcz - p).squared_norm(), vec3(0,1-z,z)}; // barycentric coordinates (0,1-w,w)
    }

    // P inside face region. Compute Q through its barycentric coordinates (u,v,w)
    const double denom = 1.0 / (va + vb + vc);
    const double v = vb * denom;
    const double w = vc * denom;
    const vec3 cpa = a + ab * v + a + ac * w;
    return {(p - cpa).squared_norm(), vec3(1.0-v-w, v, w)}; // = u*a + v*b + w*c, u = va * denom = 1.0-v-w
}
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



