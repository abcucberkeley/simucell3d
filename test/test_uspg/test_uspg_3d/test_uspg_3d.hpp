#ifndef DEF_TEST_USPG_3D
#define DEF_TEST_USPG_3D

#include <numeric>

#include "uspg_3d.hpp"



class uspg_3d_tester
{   

    public:
        uspg_3d_tester() = default;

        int update_dimensions_test() const;
        int place_object_test() const;
        int get_neighborhood_test() const;
        int get_grid_content_test() const;
        int update_voxel_test() const;



        

};

#endif