#ifndef DEF_USPG_3D
#define DEF_USPG_3D


#include <optional>

#include "uspg_abstract.hpp"

/*
    uspg_3d is a Uniform Space Partitionning Grid where each voxel can contain only 
    one object. 
*/


template <typename T>
class uspg_3d: public uspg_abstract
{

    private: 

        //Each element of this vector corresponds to a voxel of the discretized space
        //The objects contained in a voxel are in turn stored in a std::forward_list
        std::vector<std::optional<T>> voxel_lst_;

        friend class uspg_3d_tester;


    public:

        //Make sure the uspg cannot be default instantiated
        uspg_3d() = default;                           //default constructor
        uspg_3d(const uspg_3d& c) = delete;           //copy constructor
        uspg_3d(uspg_3d&& c) = delete;                //move constructor
        uspg_3d& operator=(const uspg_3d& c) = default;//copy assignment operator
        uspg_3d& operator=(uspg_3d&& c) = default;     //move assignment operator 


//----------------------------------------------------------------------------------------------------------------------------------------------------        
        //Construct the USPG based on its dimensions
        uspg_3d(const double min_x, const double min_y, const double min_z, 
                const double max_x, const double max_y, const double max_z, 
                const double voxel_size, const size_t nb_objects) noexcept(false){

            voxel_size_ = voxel_size;
            update_dimensions(nb_objects, min_x, min_y, min_z, max_x, max_y, max_z);
        }
//----------------------------------------------------------------------------------------------------------------------------------------------------


//----------------------------------------------------------------------------------------------------------------------------------------------------
    //Update the grid dimensions
    void update_dimensions(const size_t nb_objects, const double min_x, const double min_y, const double min_z, 
        const double max_x, const double max_y, const double max_z) noexcept(false) override{

            assert(min_x < max_x);
            assert(min_y < max_y);
            assert(min_z < max_z);

            //Clear the grid
            voxel_lst_.clear();

            //Add a small padding to make sure that the point (max_x, max_y, max_z) is in the grid
            constexpr double delta = std::numeric_limits<double>::epsilon();

            //Recompute the number of voxels in the grid
            nb_voxels_x_ = static_cast<unsigned>(std::ceil((max_x + delta - min_x) / voxel_size_));
            nb_voxels_y_ = static_cast<unsigned>(std::ceil((max_y + delta - min_y) / voxel_size_));
            nb_voxels_z_ = static_cast<unsigned>(std::ceil((max_z + delta - min_z) / voxel_size_));
            const size_t total_nb_voxels = nb_voxels_x_ * nb_voxels_y_ * nb_voxels_z_;

            //The new grid dimensions
            min_x_ = min_x - delta; 
            min_y_ = min_y - delta; 
            min_z_ = min_z - delta; 

            max_x_ = min_x + nb_voxels_x_ * voxel_size_;
            max_y_ = min_y + nb_voxels_y_ * voxel_size_;
            max_z_ = min_z + nb_voxels_z_ * voxel_size_;


            //Resize the grid
            voxel_lst_.resize(total_nb_voxels);
    }
//----------------------------------------------------------------------------------------------------------------------------------------------------


//----------------------------------------------------------------------------------------------------------------------------------------------------
    //Returns the content of the voxel
    const std::optional<T> get_voxel_content(const unsigned voxel_x_id, const unsigned voxel_y_id, const unsigned voxel_z_id) const noexcept{
        //Get the id of the voxel
        size_t voxel_id = get_voxel_index(voxel_x_id, voxel_y_id, voxel_z_id);
        return voxel_lst_[voxel_id];
    }
//----------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------------------------------------------------------------------
    //Place an object in the uspg
    void place_object(const T& object, const double pos_x, const double pos_y, const double pos_z) noexcept {
        //Get the index of the voxel where the object should be inserted
        const size_t voxel_id = get_voxel_index(pos_x, pos_y, pos_z);

        //Insert a copy of the object in the voxel
        voxel_lst_[voxel_id] = object;
    }

    void place_object(T&& object, const double pos_x, const double pos_y, const double pos_z) noexcept {
        //Get the index of the voxel where the object should be inserted
        const size_t voxel_id = get_voxel_index(pos_x, pos_y, pos_z);

        //Insert a copy of the object in the voxel
        voxel_lst_[voxel_id] = object;
    }


    void place_object(const T& object, const vec3& pos) noexcept{
        place_object(object, pos.dx(), pos.dy(), pos.dz());
    }

    void place_object(T&& object, const vec3& pos) noexcept{
        place_object(object, pos.dx(), pos.dy(), pos.dz());
    }
//----------------------------------------------------------------------------------------------------------------------------------------------------


//----------------------------------------------------------------------------------------------------------------------------------------------------
    //This method is similar to place object except that it directly takes the voxel index as input
    void update_voxel(const unsigned voxel_x_index, const unsigned voxel_y_index,const unsigned voxel_z_index,  const T& object) noexcept{
        //Get the id of the voxel
        size_t voxel_id = get_voxel_index(voxel_x_index, voxel_y_index, voxel_z_index);

        //Insert a copy of the object in the voxel
        voxel_lst_[voxel_id] = object;
    }
//----------------------------------------------------------------------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------------------------------------------------------------------
    //Given a position in space returns a list of copies the neighboring objects
    std::forward_list<T> get_neighborhood(const vec3& pos) const noexcept{
        return get_neighborhood(pos.dx(), pos.dy(), pos.dz());
    }
    
    std::forward_list<T> get_neighborhood(const double pos_x, const double pos_y, const double pos_z) const noexcept{
        //Make sure the object is in the grid
        assert(pos_x >= min_x_ && pos_x < max_x_);
        assert(pos_y >= min_y_ && pos_y < max_y_);
        assert(pos_z >= min_z_ && pos_z < max_z_);

        //Get the discretized 3D position of the object in the grid
        const unsigned voxel_x_id = static_cast<unsigned>(std::floor((pos_x -  min_x_) / voxel_size_));
        const unsigned voxel_y_id = static_cast<unsigned>(std::floor((pos_y -  min_y_) / voxel_size_));
        const unsigned voxel_z_id = static_cast<unsigned>(std::floor((pos_z -  min_z_) / voxel_size_));

        return get_neighborhood(voxel_x_id, voxel_y_id, voxel_z_id);
    }



    //Get all the objects in the voxels arround the voxel of the given object
    std::forward_list<T> get_neighborhood(const unsigned object_voxel_x_id, const unsigned object_voxel_y_id, const unsigned object_voxel_z_id) const noexcept{
        //Store the nighboring objects in this vector
        std::forward_list<T> neighboring_objects; 

        //The start and end position of the voxels that have to be visited
        const size_t start_voxel_x_id = object_voxel_x_id == 0 ? 0 : object_voxel_x_id - 1;
        const size_t start_voxel_y_id = object_voxel_y_id == 0 ? 0 : object_voxel_y_id - 1;
        const size_t start_voxel_z_id = object_voxel_z_id == 0 ? 0 : object_voxel_z_id - 1;

        const size_t end_voxel_x_id = object_voxel_x_id == nb_voxels_x_ - 1 ? nb_voxels_x_ : object_voxel_x_id + 2;
        const size_t end_voxel_y_id = object_voxel_y_id == nb_voxels_y_ - 1 ? nb_voxels_y_ : object_voxel_y_id + 2;
        const size_t end_voxel_z_id = object_voxel_z_id == nb_voxels_z_ - 1 ? nb_voxels_z_ : object_voxel_z_id + 2;

        //Use a triple nested for loop to visit all the surrounding voxels
        for(size_t voxel_x_id =  start_voxel_x_id; voxel_x_id < end_voxel_x_id; voxel_x_id++){
        for(size_t voxel_y_id =  start_voxel_y_id; voxel_y_id < end_voxel_y_id; voxel_y_id++){
        for(size_t voxel_z_id =  start_voxel_z_id; voxel_z_id < end_voxel_z_id; voxel_z_id++){

            //Get the id of the current voxel
            const size_t voxel_id = voxel_z_id * nb_voxels_x_ * nb_voxels_y_ + voxel_y_id * nb_voxels_x_ + voxel_x_id;

            //Get the content of the voxel
            std::optional<T> voxel_content = voxel_lst_[voxel_id];

            if(voxel_content) neighboring_objects.push_front(voxel_content.value());
            
        }}}

        return neighboring_objects;
    }
//----------------------------------------------------------------------------------------------------------------------------------------------------


//----------------------------------------------------------------------------------------------------------------------------------------------------
    //Return the total content of the grid
    std::forward_list<T> get_grid_content() const noexcept{

        //Store a copy of all the element contained in the grid
        std::forward_list<T> grid_content; 

        //Use a triple nested for loop to visit all the surrounding voxels
        for(size_t voxel_x_id =  0; voxel_x_id < nb_voxels_x_; voxel_x_id++){
        for(size_t voxel_y_id =  0; voxel_y_id < nb_voxels_y_; voxel_y_id++){
        for(size_t voxel_z_id =  0; voxel_z_id < nb_voxels_z_; voxel_z_id++){

            //Get the id of the current voxel
            const size_t voxel_id = voxel_z_id * nb_voxels_x_ * nb_voxels_y_ + voxel_y_id * nb_voxels_x_ + voxel_x_id;

            //Get the content of the voxel
            std::optional<T> voxel_content = voxel_lst_[voxel_id];

            if(voxel_content) grid_content.push_front(voxel_content.value());
            
        }}}

        return grid_content;
    }
//----------------------------------------------------------------------------------------------------------------------------------------------------


};


#endif