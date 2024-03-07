#ifndef DEF_USPG_ABSTRACT
#define DEF_USPG_ABSTRACT

#include <vector>
#include <forward_list>
#include <string>
#include <iterator>

#include "utils.hpp"
#include "custom_exception.hpp"
#include "vec3.hpp"


/*
    USPG stands for Uniform Space Partitionning Grid. This is a classical optimization technique 
    to find neighborhing objects in a space. To learn more about the theory, please refer to the excellent 
    book: Real time collision detection by Christer Ericson (Chapter 7.1) 
*/

class uspg_abstract
{

    protected: 
        //The origin of the grid
        double min_x_, min_y_, min_z_;
        double max_x_, max_y_, max_z_;

        //The number of voxels in the different directions
        unsigned nb_voxels_x_, nb_voxels_y_, nb_voxels_z_;

        //The side length of the voxels
        double voxel_size_;

        friend class contact_model_face_face;

    public:

//----------------------------------------------------------------------------------------------------------------------------------------------------
    //Make sure the uspg_abstract cannot be default instantiated
    uspg_abstract() = default;                        //default constructor
    uspg_abstract(const uspg_abstract& c) = delete;           //copy constructor
    uspg_abstract(uspg_abstract&& c) = delete;                //move constructor
    uspg_abstract& operator=(const uspg_abstract& c) = default;//copy assignment operator
    uspg_abstract& operator=(uspg_abstract&& c) = default;     //move assignment operator 

//----------------------------------------------------------------------------------------------------------------------------------------------------
    double get_voxel_size() const noexcept {return voxel_size_;}
    std::array<double, 3>   get_min_corner() const noexcept {return {min_x_, min_y_, min_z_};}
    std::array<double, 3>   get_max_corner() const noexcept {return {max_x_, max_y_, max_z_};}
    std::array<unsigned, 3> get_nb_voxels() const noexcept {return {nb_voxels_x_, nb_voxels_y_, nb_voxels_z_};}
//----------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------------------------------------------------------------------
    std::array<unsigned, 3> get_3d_voxel_index(const vec3& pos) const noexcept{
        return get_3d_voxel_index(pos.dx(), pos.dy(), pos.dz());
    }

    std::array<unsigned, 3> get_3d_voxel_index(const double pos_x, const double pos_y, const double pos_z) const noexcept{
        //Make sure the object is in the grid
        assert(pos_x >= min_x_ && pos_x <= max_x_);
        assert(pos_y >= min_y_ && pos_y <= max_y_);
        assert(pos_z >= min_z_ && pos_z <= max_z_);

        //Get the voxel in the grid that contains the point. The first part with almost equal is needed in case the 
        //points is located at the boundary of the grid
        const unsigned voxel_x_id = static_cast<unsigned>(std::floor((pos_x -  min_x_) / voxel_size_));
        const unsigned voxel_y_id = static_cast<unsigned>(std::floor((pos_y -  min_y_) / voxel_size_));
        const unsigned voxel_z_id = static_cast<unsigned>(std::floor((pos_z -  min_z_) / voxel_size_));

        return {voxel_x_id, voxel_y_id, voxel_z_id};
        
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------


//----------------------------------------------------------------------------------------------------------------------------------------------------
    //Given a position in space, indicate in which voxel the object is contained
    size_t get_voxel_index(const vec3& pos) const noexcept{
        const auto [pos_x, pos_y, pos_z] = pos.to_array();
        return get_voxel_index(pos_x, pos_y, pos_z);
    }

    
    size_t get_voxel_index(const double pos_x, const double pos_y, const double pos_z) const noexcept{

        const auto [voxel_x_id, voxel_y_id, voxel_z_id] = get_3d_voxel_index(pos_x, pos_y, pos_z);

        return get_voxel_index(voxel_x_id, voxel_y_id, voxel_z_id);
    }
    
    size_t get_voxel_index(const unsigned voxel_x_id, const unsigned voxel_y_id, const unsigned voxel_z_id) const noexcept{
        assert(voxel_x_id < nb_voxels_x_);
        assert(voxel_y_id < nb_voxels_y_);
        assert(voxel_z_id < nb_voxels_z_);

        //The position of the object in the 1D grid
        const size_t voxel_id = voxel_z_id * nb_voxels_x_ * nb_voxels_y_ + voxel_y_id * nb_voxels_x_ + voxel_x_id;
        return voxel_id;
    }
//----------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------------------------------------------------------------------
    virtual void update_dimensions(const size_t nb_objects, const double min_x, const double min_y, const double min_z, 
    const double max_x, const double max_y, const double max_z) noexcept(false) = 0;
//----------------------------------------------------------------------------------------------------------------------------------------------------
        
//----------------------------------------------------------------------------------------------------------------------------------------------------
    //Make sure the memory usage does not exceed a certain amount of the available
    //memory RAM
    //virtual void check_memory_limit(const size_t nb_objects, const double fraction_of_RAM = 1.) const noexcept(false) = 0; 
//----------------------------------------------------------------------------------------------------------------------------------------------------




};

#endif