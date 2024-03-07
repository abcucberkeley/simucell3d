#include "test_uspg_4d.hpp"

/*

Test the method of the cell class

*/


//---------------------------------------------------------------------------------------------------------
int uspg_4d_tester::update_dimensions_test() const{

    //The dimensions of the grid
    double min_x = 0.; double max_x = 10.;
    double min_y = 0.; double max_y = 9.5;
    double min_z = 0.; double max_z = 9.;
    unsigned nb_objects = 10;
    double voxel_size = 1.; 
    
    //Create a uspg_4d grid with int
    uspg_4d<int> grid(min_x, min_y, min_z, max_x, max_y, max_z, voxel_size, nb_objects);

    //Check the dimensions of the grid
    bool t1 = grid.nb_voxels_x_ == 10;
    bool t2 = grid.nb_voxels_y_ == 10;
    bool t3 = grid.nb_voxels_z_ == 9;
    bool t4 = grid.voxel_size_ == voxel_size;

    //Check the size of the grid voxel lst
    const size_t total_nb_voxels = 10 * 10 * 9;
    bool t5 = grid.voxel_lst_.size() == total_nb_voxels;

    return !(t1 && t2 && t3 && t4 && t5);
}
//---------------------------------------------------------------------------------------------------------



//---------------------------------------------------------------------------------------------------------
int uspg_4d_tester::place_object_test() const{

    //The dimensions of the grid
    double min_x = 0.; double max_x = 3.;
    double min_y = 0.; double max_y = 3.;
    double min_z = 0.; double max_z = 3.;
    unsigned nb_objects = 10;
    double voxel_size = 1.; 

    //Create a grid with 27 voxels
    uspg_4d<int> grid(min_x, min_y, min_z, max_x, max_y, max_z, voxel_size, nb_objects);

    //Create and store some integers. Those have to be stored otherwise their pointers
    //because the grid store pointers and not copies of the objects inserted
    std::vector<int> int_lst(3);
    std::iota(int_lst.begin(), int_lst.end(), 1);

    //Place the 3 integers at 3 different locations in the grid
    grid.place_object(int_lst[0], 0., 0., 0.);
    grid.place_object(int_lst[1], (unsigned) 1, (unsigned) 1, (unsigned) 1);
    grid.place_object(int_lst[2], 2.5, 2.5, 2.5);

    //Make sure the objects are located where expected
    const auto voxel_000_content  = grid.get_voxel_content(0, 0, 0);
    const auto voxel_111_content  = grid.get_voxel_content(1, 1, 1);
    const auto voxel_222_content  = grid.get_voxel_content(2, 2, 2);


    //Make sure the voxel has one object in it
    bool t1 = std::distance(voxel_000_content.begin(), voxel_000_content.end()) == 1;
    bool t2 = std::distance(voxel_111_content.begin(), voxel_111_content.end()) == 1;
    bool t3 = std::distance(voxel_222_content.begin(), voxel_222_content.end()) == 1;

    //Make sure the content is correct
    bool t4 = false, t5 = false, t6 = false;
    if(t1) t4 = (voxel_000_content.front()) == 1;
    if(t2) t5 = (voxel_111_content.front()) == 2;
    if(t3) t6 = (voxel_222_content.front()) == 3;

    return !(t1 && t2 && t3 && t4 && t5 && t6);
}





//---------------------------------------------------------------------------------------------------------



//---------------------------------------------------------------------------------------------------------
int uspg_4d_tester::get_neighborhood_test() const{

    //The dimensions of the grid
    double min_x = 0.; double max_x = 3.;
    double min_y = 0.; double max_y = 3.;
    double min_z = 0.; double max_z = 3.;
    unsigned nb_objects = 10;
    double voxel_size = 1.; 

    //Create a grid with 27 voxels
    uspg_4d<int> grid(min_x, min_y, min_z, max_x, max_y, max_z, voxel_size, nb_objects);

    //Create and store some integers. Those have to be stored otherwise their pointers
    //because the grid store pointers and not copies of the objects inserted
    std::vector<int> int_lst(27);
    std::iota(int_lst.begin(), int_lst.end(), 1);

    int counter = 0;
    //In each voxel insert 2 integers
    for(int z = 0; z < 3 ; z++){
    for(int y = 0; y < 3 ; y++){
    for(int x = 0; x < 3 ; x++){
        
        grid.place_object(int_lst[counter], x, y, z);

        counter++;
    }}}

    //Get all the objects around the voxel 000 and convert the content into a vector
    auto neighborhood_000_lst  = grid.get_neighborhood((unsigned) 0, (unsigned) 0, (unsigned) 0);
    std::vector<int> neighborhood_000_vec;
    for(auto elt: neighborhood_000_lst) neighborhood_000_vec.push_back(elt);
    std::sort(neighborhood_000_vec.begin(), neighborhood_000_vec.end());

    //Make sure the correct objects are in the neighborhood of the vector 000
    std::vector<int> correct_neighborhood_000{1, 2, 4, 5, 10, 11, 13, 14};
    bool t1 = std::equal(correct_neighborhood_000.begin(), correct_neighborhood_000.end(), neighborhood_000_vec.begin());


    //Repeat the same operation with the voxel 111 which is at the center of the grid
    auto neighborhood_111_lst  = grid.get_neighborhood((unsigned) 1, (unsigned) 1, (unsigned) 1);
    std::vector<int> neighborhood_111_vec;
    for(auto elt: neighborhood_111_lst) neighborhood_111_vec.push_back(elt);
    std::sort(neighborhood_111_vec.begin(), neighborhood_111_vec.end());

    //This voxel should be in contact with all the other voxels
    bool t2 = std::equal(neighborhood_111_vec.begin(), neighborhood_111_vec.end(), int_lst.begin());


    return !(t1 && t2);
}
//---------------------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------------------
int uspg_4d_tester::get_grid_content_test() const{
    //The dimensions of the grid
    double min_x = 0.; double max_x = 3.;
    double min_y = 0.; double max_y = 3.;
    double min_z = 0.; double max_z = 3.;
    unsigned nb_objects = 10;
    double voxel_size = 1.; 

    //Create a grid with 27 voxels
    uspg_4d<int> grid(min_x, min_y, min_z, max_x, max_y, max_z, voxel_size, nb_objects);

    //Create and store some integers. Those have to be stored otherwise their pointers
    //because the grid store pointers and not copies of the objects inserted
    std::vector<int> int_lst(27);
    std::iota(int_lst.begin(), int_lst.end(), 1);

    int counter = 0;
    //In each voxel insert 2 integers
    for(int z = 0; z < 3 ; z++){
    for(int y = 0; y < 3 ; y++){
    for(int x = 0; x < 3 ; x++){
        
        grid.place_object(int_lst[counter], x, y, z);

        counter++;
    }}}

    const std::forward_list<int>  grid_content = grid.get_grid_content();

    //Create a copy of the grid content
    std::vector<int> grid_content_copy;
    for(auto obj: grid_content) grid_content_copy.push_back(obj);

    //Sort the grid content
    std::sort(grid_content_copy.begin(), grid_content_copy.end());

    return !(std::equal(int_lst.begin(), int_lst.end(), grid_content_copy.begin()));
}
//---------------------------------------------------------------------------------------------------------




//---------------------------------------------------------------------------------------------------------
// The main function
int main (int argc, char** argv){

    //Check that the command line input is correctly formatted
    assert(argc == 2); 

    //Get the name of the test to run
    std::string test_name = argv[1];
    
    //Run the selected test
    uspg_4d_tester tester;

    if (test_name == "update_dimensions_test")     return tester.update_dimensions_test();
    if (test_name == "get_neighborhood_test")           return tester.get_neighborhood_test();
    if (test_name == "place_object_test")               return tester.place_object_test();
    if (test_name == "get_grid_content_test")           return tester.get_grid_content_test();




    std::cout << "TEST NAME :" << test_name << " DOES NOT EXIST" << std::endl;
    return 1;

}
//---------------------------------------------------------------------------------------------------------
