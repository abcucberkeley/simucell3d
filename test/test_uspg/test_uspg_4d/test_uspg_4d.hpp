#ifndef DEF_TEST_USPG_4D
#define DEF_TEST_USPG_4D

#include <numeric>

#include "uspg_4d.hpp"



class uspg_4d_tester
{   

    public:
        uspg_4d_tester() = default;

        int update_dimensions_test() const;
        int place_object_test() const;
        int get_neighborhood_test() const;
        int get_grid_content_test() const;

        

};

#endif