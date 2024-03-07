#include "statistics_writer.hpp"


/*
The statistics of the cells (volumes, areas, etc.) can either be saved in a csv file or in a string. 
The string can be returned at the end of the simulation. This is useful for passing the simulation 
statistics to the python interpreter if simucell3d is called via the python interface. New cell attributes 
can be added to the CSV file or the final string by simply adding the corresponding function in 
the vector file_data_mapper_lst in the file ./include/io/mesh_data.hpp. 
*/





//-----------------------------------------------------------------------------------------------
//The trivial constructor to write the data in a CSV file
csv_file_statistics_writer::csv_file_statistics_writer(const std::string& output_file_path) 
noexcept(false):  abstract_statistics_writer(), output_file_path_(output_file_path){

    starting_time_point_ = std::chrono::system_clock::now();

    // Create the output filestream object
    std::ofstream file_stream = std::ofstream(output_file_path_);
    if(!file_stream){throw intialization_exception("Unable to write csv file: " + output_file_path_);}

    //The first three columns are reserved for the iteration number, the computation time and the simulation time
    file_stream << "iteration" << sep << "computation_time_(hh::mm:ss)" << sep << "simulation_time" << sep;

    //Loop over the functions defined in the vector file_data_mapper_lst which is located in mesh_data.hpp
    for(auto& mapper : file_data_mapper_lst){
        //Write the header of the csv file
        file_stream << mapper.value_name_ << sep;
    }

    file_stream << "\n";
    file_stream.close();
}
//-----------------------------------------------------------------------------------------------




//-----------------------------------------------------------------------------------------------
void csv_file_statistics_writer::write_data(
    const unsigned iteration, 
    const double simulation_time,  
    const std::vector<cell_ptr>& cell_lst
) noexcept(false){

    //Open the file 
    std::ofstream file_stream = std::ofstream(output_file_path_, std::ios::app);
    if(!file_stream){throw intialization_exception("Unable to open csv file: " + output_file_path_);}

    //Compute the computation time
    auto elapsed_seconds = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now() - starting_time_point_).count();

    const size_t hh =  static_cast<size_t>(elapsed_seconds/3600);
    const size_t mm = (static_cast<size_t>(elapsed_seconds/60))%60;
    const size_t ss =  static_cast<size_t>(elapsed_seconds%60);

    //Stores the computation time in a string
    char buffer[50];
    int n = sprintf(buffer, "%02ld:%02ld:%02ld", hh, mm, ss);
    std::string computation_time_str(buffer);

    //Loop over the cells
    for(cell_ptr c: cell_lst){

        //Start by saving the iteration number
        file_stream << std::to_string(iteration) << sep << computation_time_str << sep << format_number(simulation_time, "%.2e") << sep;

        //Loop over the functions defined in the vector file_data_mapper_lst which is located in mesh_data.hpp
        for(auto& mapper : file_data_mapper_lst){

            //Save the data of the cell
            file_stream << mapper.value_extractor_(c) << sep;
        }

        file_stream << "\n";
    }


    //Close the file 
    file_stream.close();

}
//-----------------------------------------------------------------------------------------------




//-----------------------------------------------------------------------------------------------
//The trivial constructor
string_statistics_writer::string_statistics_writer() 
noexcept(false):  abstract_statistics_writer(){

    starting_time_point_ = std::chrono::system_clock::now();

    //The first column is reserved for the iteration number
    str_stream << "iteration" << sep << "computation_time_(hh::mm:ss)" << sep << "simulation_time" << sep;

    //Loop over the functions defined in the vector file_data_mapper_lst which is located in mesh_data.hpp
    for(auto& mapper : file_data_mapper_lst){
        //Write the header of the csv file
        str_stream << mapper.value_name_ << sep;
    }

    str_stream << "\n";
}
//-----------------------------------------------------------------------------------------------




//-----------------------------------------------------------------------------------------------
void string_statistics_writer::write_data(
    const unsigned iteration, 
    const double simulation_time,  
    const std::vector<cell_ptr>& cell_lst
) noexcept(false){

    //Compute the computation time
    auto elapsed_seconds = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now() - starting_time_point_).count();

    const size_t hh =  static_cast<size_t>(elapsed_seconds/3600);
    const size_t mm = (static_cast<size_t>(elapsed_seconds/60))%60;
    const size_t ss =  static_cast<size_t>(elapsed_seconds%60);

    char buffer[50];
    int n = sprintf(buffer, "%02ld:%02ld:%02ld", hh, mm, ss);
    std::string computation_time_str(buffer);

    //Loop over the cells
    for(cell_ptr c: cell_lst){

        //Start by saving the iteration number
        str_stream << std::to_string(iteration) << sep << computation_time_str << sep << format_number(simulation_time, "%.2e") << sep;

        //Loop over the functions defined in the vector file_data_mapper_lst which is located in mesh_data.hpp
        for(auto& mapper : file_data_mapper_lst){

            //Save the data of the cell
            str_stream << mapper.value_extractor_(c) << sep;
        }
        str_stream << "\n";
    }
}
//-----------------------------------------------------------------------------------------------



