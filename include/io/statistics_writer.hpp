#ifndef DEF_STATISTICS_WRITER
#define DEF_STATISTICS_WRITER

#include <string>
#include <fstream>
#include <iostream>
#include <streambuf>
#include <sstream>
#include <cassert>
#include <vector>
#include <chrono>

#include "custom_structures.hpp"
#include "custom_exception.hpp"
#include "utils.hpp"
#include "mesh_data.hpp"
#include "cell.hpp"



/*
The statistics of the cells (volumes, areas, etc.) can either be saved in a csv file or in a string. 
The string can be returned at the end of the simulation. This is useful for passing the simulation 
statistics to the python interpreter if simucell3d is called via the python interface. New cell attributes 
can be added to the CSV file or the final string by simply adding the corresponding function in 
the vector file_data_mapper_lst in the file ./include/io/mesh_data.hpp. 
*/




//-----------------------------------------------------------------------------------------------
//This is an abstract class that cannot be instantiated
class abstract_statistics_writer{
    protected:
        std::string sep = ",";

        std::chrono::time_point<std::chrono::system_clock> starting_time_point_;

    public:
        abstract_statistics_writer() = default;                                              //default constructor
        abstract_statistics_writer(const abstract_statistics_writer& v) = delete;            //copy constructor
        abstract_statistics_writer(abstract_statistics_writer&& v) = delete;                 //move constructor
        abstract_statistics_writer& operator=(const abstract_statistics_writer& v) = delete; //copy assignment operator
        abstract_statistics_writer& operator=(abstract_statistics_writer&& v) = default;     //move assignment operator 

        virtual void write_data(
            const unsigned iteration, 
            const double simulation_time,  
            const std::vector<cell_ptr>& cell_lst
        ) noexcept(false) = 0;
};
//-----------------------------------------------------------------------------------------------




//-----------------------------------------------------------------------------------------------
class csv_file_statistics_writer: public abstract_statistics_writer{
    private:

        //The path to the csv file where the data will be written
        std::string output_file_path_;

    public:
        csv_file_statistics_writer() = default;                                    //default constructor
        csv_file_statistics_writer(const csv_file_statistics_writer& v) = delete;           //copy constructor
        csv_file_statistics_writer(csv_file_statistics_writer&& v) = delete;                //move constructor
        csv_file_statistics_writer& operator=(const csv_file_statistics_writer& v) = delete;//copy assignment operator
        csv_file_statistics_writer& operator=(csv_file_statistics_writer&& v) = default;     //move assignment operator 

        //The trivial constructor
        csv_file_statistics_writer(const std::string& output_file_path) noexcept(false);


        void write_data(
            const unsigned iteration, 
            const double simulation_time,  
            const std::vector<cell_ptr>& cell_lst
        ) noexcept(false) override;


};
//-----------------------------------------------------------------------------------------------








//-----------------------------------------------------------------------------------------------
class string_statistics_writer: public abstract_statistics_writer{
    private:

        //The data will be written in this string
        std::stringstream str_stream;


    public:                              //default constructor
        string_statistics_writer(const string_statistics_writer& v) = delete;           //copy constructor
        string_statistics_writer(string_statistics_writer&& v) = delete;                //move constructor
        string_statistics_writer& operator=(const string_statistics_writer& v) = delete;//copy assignment operator
        string_statistics_writer& operator=(string_statistics_writer&& v) = default;     //move assignment operator 

        //The trivial constructor
        string_statistics_writer() noexcept(false);

        //Store all the simulation data in the string stream
        void write_data(
            const unsigned iteration, 
            const double simulation_time,  
            const std::vector<cell_ptr>& cell_lst
        ) noexcept(false) override;

        //Return the string containing the data
        std::string get_string() const noexcept(false){return str_stream.str();}
};
//-----------------------------------------------------------------------------------------------


#endif