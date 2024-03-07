
#ifndef DEF_AUTOMATIC_POLARIZER
#define DEF_AUTOMATIC_POLARIZER


#include <iostream>
#include <vector>
#include <set>
#include <array>
#include <functional>

#include <omp.h>


#include "global_configuration.hpp"

#include "utils.hpp"

#include "custom_structures.hpp"
#include "custom_exception.hpp"

#include "cell.hpp"
#include "face.hpp"
#include "node.hpp"
#include "uspg_3d.hpp"




/*The automatic polarizer manages the face types of the epithelial cells. It first detects if there 
are lumens in the geometry and then mark the faces in contact with the lumen as apical. The faces that are 
in contact with another cell are marked lateral, and the faces in contact with the exterior space are 
marked as basal. 

The automatic polarizer detects luminal, cytoplasmic and the exterior space by discretizing the 
simulation space. Then the cell surfaces are used to create boundaries in this discretized space.
and finally the different regions are detected with the Hoshen–Kopelman algorithm.
*/


//The following lines don't cost us any computation time, these are used to assign 
//to each voxel of the discretized space a label that corresponds to the region it belongs to
constexpr unsigned short boundary_region_id  = 0;
constexpr unsigned short cytoplasm_region_id = 1;
constexpr unsigned short exterior_region_id  = 2;
constexpr unsigned short lumen_region_id     = 3;



class automatic_polarizer{

        protected:

                //The unifom space partitionning grid is used to discretize the simulation space in voxels
                uspg_3d<unsigned short> grid_;
                const double grid_voxel_size_;

                //The maximum length of a ray used to detect with which region a face is in contact
                double ray_maximum_length_;
                


        public:
                automatic_polarizer() = delete;                                         //default constructor
                automatic_polarizer(const automatic_polarizer& v) = delete;           //copy constructor
                automatic_polarizer(automatic_polarizer&& v) = delete;                //move constructor
                automatic_polarizer& operator=(const automatic_polarizer& v) = delete;//copy assignment operator
                automatic_polarizer& operator=(automatic_polarizer&& v) = delete;    //move assignment operator 


                //Constructor
                automatic_polarizer(const double max_edge_length) noexcept: grid_voxel_size_(max_edge_length){};


                inline const uspg_3d<unsigned short>& get_grid() const noexcept {return grid_;};

                //Update the dimensions and the content of the discretization grid
                void polarize_faces(const std::vector<cell_ptr>& cell_lst) noexcept(false);

                //Update the dimensions of the discretization grid
                void update_grid_dimensions(const std::vector<cell_ptr>& cell_lst) noexcept;

                //Mark all the boundary voxels
                void mark_boundary_voxels(const std::vector<cell_ptr>& cell_lst) noexcept;

                //Mark all the voxels that are localized in a cytoplasm
                void mark_cytoplasm_voxels(const std::vector<cell_ptr>& cell_lst) noexcept;

                //Mark all the voxels that are localized in the exterior space 
                void mark_exterior_space_voxels(const std::vector<cell_ptr>& cell_lst);

                //Mark all the voxels that are localized in a lumen 
                void mark_luminal_voxels(const std::vector<cell_ptr>& cell_lst);

                //Given a voxel with a certain label and position in the grid, this method expand the label to the 
                //surrounding voxels until it reaches the boundary of the grid or another voxel with a different label
                void expand_voxel_label(const std::array<unsigned,3>& voxel_index,const unsigned short label) noexcept;
                
                //Given the index of a voxel, this method returns the indices of all the neighboring voxels
                //It takes into account the boundary of the grid
                std::set<std::array<unsigned,3>> get_neighbor_voxels(std::array<unsigned, 3> voxel_index) const noexcept;

                //As it name indicates, use the discretize space to polarize all the faces of epithelial cells 
                void polarize_faces_based_on_contacts_with_discretized_space(cell_ptr c) const noexcept(false);

                //Determine the region in contact with the face by shooting a ray from the face in the direction of its normal
                //and returns the region id of the first voxel that is not a boundary voxel along the ray
                unsigned short get_region_in_contact_with_face(cell_ptr c, const face& f) const noexcept(false);



};

#endif