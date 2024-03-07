#include "automatic_polarizer.hpp"



//----------------------------------------------------------------------------------------------------------------------------------------------------
void automatic_polarizer::polarize_faces(const std::vector<cell_ptr>& cell_lst) noexcept(false){

    //Start by updating the dimensions of the grid
    update_grid_dimensions(cell_lst);

    //Mark all the boundary voxels
    mark_boundary_voxels(cell_lst);

    //Mark all the voxels that are localized in a cytoplasm
    mark_cytoplasm_voxels(cell_lst);

    //Mark all the voxels that are localized in the exterior space and the lumen
    mark_exterior_space_voxels(cell_lst);

    //Mark all the voxels that belong to a luminal region
    mark_luminal_voxels(cell_lst);

    //As it name indicates, use the discretize space to update the type of the faces
    for(cell_ptr c: cell_lst){
        polarize_faces_based_on_contacts_with_discretized_space(c);
    }
}
//----------------------------------------------------------------------------------------------------------------------------------------------------


//----------------------------------------------------------------------------------------------------------------------------------------------------
//Update the dimensions of the discretization grid, that will be used to segment the space in different regions
void automatic_polarizer::update_grid_dimensions(const std::vector<cell_ptr>& cell_lst) noexcept{

    //Compute the dimensions of the grid
    //Get the global min and max coordinates of the mesh
    double min_x = std::numeric_limits<double>::infinity();
    double min_y = std::numeric_limits<double>::infinity();
    double min_z = std::numeric_limits<double>::infinity();
    double max_x = std::numeric_limits<double>::lowest();
    double max_y = std::numeric_limits<double>::lowest();
    double max_z = std::numeric_limits<double>::lowest();

    //Loop over all the nodes in the mesh
    for(cell_ptr c: cell_lst){
        for(const node& n: c->get_node_lst()){
            if(n.pos().dx() < min_x){min_x = n.pos().dx();}
            if(n.pos().dy() < min_y){min_y = n.pos().dy();}
            if(n.pos().dz() < min_z){min_z = n.pos().dz();}
            if(n.pos().dx() > max_x){max_x = n.pos().dx();}
            if(n.pos().dy() > max_y){max_y = n.pos().dy();}
            if(n.pos().dz() > max_z){max_z = n.pos().dz();}
        }
    }

    //Add an extra layer of voxels around the geometry
    min_x -= grid_voxel_size_ * 2.;
    min_y -= grid_voxel_size_ * 2.;
    min_z -= grid_voxel_size_ * 2.;
    max_x += grid_voxel_size_ * 2.;
    max_y += grid_voxel_size_ * 2.;
    max_z += grid_voxel_size_ * 2.;
   
    //Compute the number of voxels in the grid
    constexpr double delta = std::numeric_limits<double>::epsilon();
    const unsigned nb_voxels_x_ = static_cast<unsigned>(std::ceil((max_x + delta - min_x) / grid_voxel_size_));
    const unsigned nb_voxels_y_ = static_cast<unsigned>(std::ceil((max_y + delta - min_y) / grid_voxel_size_));
    const unsigned nb_voxels_z_ = static_cast<unsigned>(std::ceil((max_z + delta - min_z) / grid_voxel_size_));
    const size_t total_nb_voxels = nb_voxels_x_ * nb_voxels_y_ * nb_voxels_z_;

    //Empties the space partitioning grid and update its dimensions
    grid_ = uspg_3d<unsigned short>(min_x,min_y,min_z, max_x,max_y,max_z, grid_voxel_size_, total_nb_voxels);

    const auto [max_voxel_x_id, max_voxel_y_id, max_voxel_z_id] = grid_.get_nb_voxels();

    //Update the ray maximum length
    ray_maximum_length_ =  grid_voxel_size_  * std::sqrt(3) * static_cast<double>(std::max({nb_voxels_x_, nb_voxels_y_, nb_voxels_z_}) + 1);


}
//----------------------------------------------------------------------------------------------------------------------------------------------------


//----------------------------------------------------------------------------------------------------------------------------------------------------
//Mark all the boundary voxels
void automatic_polarizer::mark_boundary_voxels(const std::vector<cell_ptr>& cell_lst) noexcept{
    //If a node is located in the voxel then it is marked as boundary
    for(cell_ptr c: cell_lst){
        for(const node& n: c->get_node_lst()){
            grid_.place_object(boundary_region_id, n.pos());
        }
    }  
}
//----------------------------------------------------------------------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------------------------------------------------------------------
//Mark all the voxels that are localized in a cytoplasm
//This method will fail if the cells are strongly concave and their centroids are outside 
//of the cell volume/domain
void automatic_polarizer::mark_cytoplasm_voxels(const std::vector<cell_ptr>& cell_lst) noexcept{

    //Loop over the cells
    for(cell_ptr c: cell_lst){

        //Get the centroid of the cell
        const vec3 centroid = c->get_centroid();

        //Get the index of the voxel where the centroid is located
        const std::array<unsigned ,3> voxel_index = grid_.get_3d_voxel_index(centroid);

        //Make this voxel a cytoplasm voxel
        grid_.place_object(cytoplasm_region_id, centroid);

        //Expand the label to the surrounding voxels
        expand_voxel_label(voxel_index, cytoplasm_region_id);
    }
}
//----------------------------------------------------------------------------------------------------------------------------------------------------


//----------------------------------------------------------------------------------------------------------------------------------------------------
//Mark all the voxels that are localized in the exterior space 
void automatic_polarizer::mark_exterior_space_voxels(const std::vector<cell_ptr>& cell_lst){

    //We use the fact that we have added one layer of voxels around the geometry
    //therefore we can select one voxel that is in the exterior space e.g. voxel (0,0,0)
    //and then expand the label to the surrounding voxels
    grid_.update_voxel(0,0,0, exterior_region_id);

    //Expand the exterior label to the surrounding voxels
    expand_voxel_label({0,0,0}, exterior_region_id);
}
//----------------------------------------------------------------------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------------------------------------------------------------------
//Mark all the voxels that are localized in a lumen 
void automatic_polarizer::mark_luminal_voxels(const std::vector<cell_ptr>& cell_lst){

    //Now we can exploit the fact that the remaining unmarked voxels are necessarily voxels 
    //that are within a lumen. 

    for(unsigned voxel_x_id = 0; voxel_x_id < grid_.get_nb_voxels()[0]; voxel_x_id++){
    for(unsigned voxel_y_id = 0; voxel_y_id < grid_.get_nb_voxels()[1]; voxel_y_id++){
    for(unsigned voxel_z_id = 0; voxel_z_id < grid_.get_nb_voxels()[2]; voxel_z_id++){

        //Get the content of the voxel
        const std::optional<unsigned short> voxel_content = grid_.get_voxel_content(voxel_x_id, voxel_y_id, voxel_z_id);

        //If the voxel is not already marked then we add it to the list of neighbors
        if(!voxel_content.has_value()){
            grid_.update_voxel(voxel_x_id, voxel_y_id, voxel_z_id, lumen_region_id);
        }
    }}}
}
//----------------------------------------------------------------------------------------------------------------------------------------------------


//----------------------------------------------------------------------------------------------------------------------------------------------------
//Given a voxel with a certain label and position in the grid, this method expand the label to the 
//surrounding voxels until it reaches the boundary of the grid or a boundary voxel
void automatic_polarizer::expand_voxel_label(
        const std::array<unsigned,3>& voxel_index,
        const unsigned short label    
    ) noexcept
{
    
    //This vector contains the indices of all the voxels that need to be visited
    //Initially it contains all the neighbors of the voxel given in argument
    std::set<std::array<unsigned,3>> voxel_to_visit_set = get_neighbor_voxels(voxel_index);

    //We loop over all the voxels in the list
    while(voxel_to_visit_set.size() > 0){

        //Pop out the first element
        std::array<unsigned,3> current_voxel_index = *(voxel_to_visit_set.begin());
        voxel_to_visit_set.erase(voxel_to_visit_set.begin());

        //Get the content of the voxel
        const std::optional<unsigned short> voxel_content = grid_.get_voxel_content(current_voxel_index[0], current_voxel_index[1], current_voxel_index[2]);

        //If the voxel has not been marked yet then we update its label and add its neighbors to the list of voxels to visit
        if(!voxel_content.has_value()){

            //Update the voxel with the label given in argument
            grid_.update_voxel(current_voxel_index[0], current_voxel_index[1], current_voxel_index[2], label);

            //Get the neighbors of the current voxel that are not already marked
            //and add them to the list of voxels to visit
            std::set<std::array<unsigned,3>> current_voxel_neighbors = get_neighbor_voxels(current_voxel_index);
            voxel_to_visit_set.insert(current_voxel_neighbors.begin(), current_voxel_neighbors.end());
        }
    }
}
//----------------------------------------------------------------------------------------------------------------------------------------------------


//----------------------------------------------------------------------------------------------------------------------------------------------------
//Given the index of a voxel, this method returns the indices of all the neighboring voxels
//It takes into account the boundary of the grid. 
std::set<std::array<unsigned,3>> automatic_polarizer::get_neighbor_voxels(std::array<unsigned, 3> voxel_index) const noexcept{

    const auto [voxel_x_id, voxel_y_id, voxel_z_id] = voxel_index;

    //Get the maximum index of the voxels in the grid
    const auto [max_voxel_x_id, max_voxel_y_id, max_voxel_z_id] = grid_.get_nb_voxels();

    //The start and end position of the voxels that have to be visited
    const unsigned start_voxel_x_id = (voxel_x_id == 0) ? 0 : voxel_x_id - 1;
    const unsigned start_voxel_y_id = (voxel_y_id == 0) ? 0 : voxel_y_id - 1;
    const unsigned start_voxel_z_id = (voxel_z_id == 0) ? 0 : voxel_z_id - 1;

    const unsigned end_voxel_x_id = (voxel_x_id == max_voxel_x_id - 1) ? voxel_x_id : voxel_x_id + 1;
    const unsigned end_voxel_y_id = (voxel_y_id == max_voxel_y_id - 1) ? voxel_y_id : voxel_y_id + 1;
    const unsigned end_voxel_z_id = (voxel_z_id == max_voxel_z_id - 1) ? voxel_z_id : voxel_z_id + 1;


    //Create a vector with the indices of the 6 surrounding voxels
    //If there are duplicate voxels then they will be removed by the set
    std::set<std::array<unsigned, 3>> neighbor_voxel_set{
        {start_voxel_x_id,  voxel_y_id,         voxel_z_id},
        {end_voxel_x_id,    voxel_y_id,         voxel_z_id},

        {voxel_x_id,        start_voxel_y_id,   voxel_z_id},
        {voxel_x_id,        end_voxel_y_id,     voxel_z_id},

        {voxel_x_id,        voxel_y_id,         start_voxel_z_id},
        {voxel_x_id,        voxel_y_id,         end_voxel_z_id}
    };


    return neighbor_voxel_set;
}
//----------------------------------------------------------------------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------------------------------------------------------------------

//As it name indicates, use the discretize space to polarize all the faces of epithelial cells 
void automatic_polarizer::polarize_faces_based_on_contacts_with_discretized_space(cell_ptr c) const noexcept(false){
 
    //If the cell is an epithelial cell
    if(c->get_cell_type_id() == 0){

        //Loop over the faces of the epithelial cell
        for(face& f: c->face_lst_){

            if(f.is_used()){
                 /*
                region_id = 1: cytoplasm
                region_id = 2: exterior space
                region_id = 3: lumen
                */
                unsigned short region_id = get_region_in_contact_with_face(c, f);
                assert(region_id == cytoplasm_region_id || region_id == exterior_region_id || region_id == lumen_region_id);

                //If the face is in contact with the cytoplasm
                if(region_id == lumen_region_id){
                    f.set_face_type_id(1);
                }
                else if(region_id == cytoplasm_region_id){
                    f.set_face_type_id(1);
                }
                else if(region_id == exterior_region_id){
                    f.set_face_type_id(0);
                }
                else{
                    std::runtime_error("Automatic polarization: The face is in contact with an unknown region");
                }
            }
        }
    }
}


//----------------------------------------------------------------------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------------------------------------------------------------------
//Determine the region in contact with the face by shooting a ray from the face in the direction of its normal
//and returns the region id of the first voxel that is not a boundary voxel along the ray
unsigned short automatic_polarizer::get_region_in_contact_with_face(cell_ptr c, const face& f) const noexcept(false){

    //Get the position of the first node of the face
    const node& n1 = c->get_node_lst()[f.n1_id()];
    const vec3& n1_pos = n1.pos();

    //Get the normal of the face
    const vec3& face_normal = f.get_normal();

    //The origin of the ray is the position of the first point of the face
    const vec3 ray_origin = n1_pos;

    //The end point of the ray is computed by adding the normal of the face to the position of the first point of the face
    const vec3 ray_end = n1_pos + (face_normal * ray_maximum_length_);

    //The code that follows was adapted from the following source:
    //Real time collision detection by Christer Ericson, chapter 7.4.2 (Uniform Grid Intersection Test)
    const auto [x1, y1, z1] = ray_origin.to_array();
    const auto [x2, y2, z2] = ray_end.to_array();

    //Get the index of the first voxel that is intersected by the ray
    auto [i,j,k] = grid_.get_3d_voxel_index(ray_origin);

    // Determine in which primary direction to step
    unsigned di = ((x1 < x2) ? 1 : ((x1 > x2) ? -1 : 0));
    unsigned dj = ((y1 < y2) ? 1 : ((y1 > y2) ? -1 : 0));
    unsigned dk = ((z1 < z2) ? 1 : ((z1 > z2) ? -1 : 0));


    //Determine how far we should advance along the ray to be in the next voxel of the grid which is along the ray
    double tx, ty, tz;

    //Start by the x-direction
    const double min_x = grid_voxel_size_ * std::floor(x1 / grid_voxel_size_);
    const double min_y = grid_voxel_size_ * std::floor(y1 / grid_voxel_size_);
    const double min_z = grid_voxel_size_ * std::floor(z1 / grid_voxel_size_);

    const double max_x = min_x + grid_voxel_size_;
    const double max_y = min_y + grid_voxel_size_;
    const double max_z = min_z + grid_voxel_size_;


    if(almost_equal(face_normal.dx(), 0.)){tx = std::numeric_limits<double>::infinity();}
    else{tx = (face_normal.dx() > 0.) ? (x1 - min_x) / std::abs(x2 - x1) : tx = (max_x - x1) / std::abs(x2 - x1);}

    if(almost_equal(face_normal.dy(), 0.)){ty = std::numeric_limits<double>::infinity();}
    else{ty = (face_normal.dy() > 0.) ? (y1 - min_y) / std::abs(y2 - y1) : ty = (max_y - y1) / std::abs(y2 - y1);}

    if(almost_equal(face_normal.dz(), 0.)){tz = std::numeric_limits<double>::infinity();}
    else{tz = (face_normal.dz() > 0.) ? (z1 - min_z) / std::abs(z2 - z1) : tz = (max_z - z1) / std::abs(z2 - z1);}

    // Determine deltax/deltay, how far (in units of t) one must step along the directed line segment for the horizontal/vertical
    // movement (respectively) to equal the width/height of a cell
    double delta_tx = grid_voxel_size_ / std::abs(x2 - x1);
    double delta_ty = grid_voxel_size_ / std::abs(y2 - y1);
    double delta_tz = grid_voxel_size_ / std::abs(z2 - z1);

    //If the face normal is colinear with a row of boundary buckets the first region which is not a boundary will
    //be the exterior space, the face will be marked as basal, even if it should have it should have been lateral.
    size_t max_iteration = std::max(grid_.get_nb_voxels()[0], std::max(grid_.get_nb_voxels()[1], grid_.get_nb_voxels()[2])) * 3;
    size_t iteration = 0;

   while(iteration < max_iteration) {
        iteration++;

        //Get the label of the bucket
        std::optional<unsigned short> region_id_opt = grid_.get_voxel_content(i,j,k);
        assert(region_id_opt.has_value());
        const unsigned short region_id = region_id_opt.value();

        if (region_id == boundary_region_id){
            if (tx <= ty && tx <= tz) { // tx smallest, step in x
                tx += delta_tx;
                i += di;
            }
            else if (ty <= tx && ty <= tz) { // ty smallest, step in y
                ty += delta_ty;
                j += dj;
            }
            else{ // tz smallest, step in z
                tz += delta_tz;
                k += dk;
            }
            continue;
        }

        //If the voxel is not a boundary voxel, we return its label
        return region_id;
    }

    throw automatic_polarization_exception("The ray did not intersect any voxel in the grid");
}
//---------------------------------------------------------------------------------------------------------------








//----------------------------------------------------------------------------------------------------------------------------------------------------
