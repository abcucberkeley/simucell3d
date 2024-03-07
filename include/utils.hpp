#ifndef DEF_UTILS
#define DEF_UTILS


#include <type_traits>
#include <functional>
#include <limits>
#include <algorithm>
#include <exception>
#include <string>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <cassert>
#include <array>
#include <omp.h>

/*
Contains a bunch of useful functions
*/





//Check if 2 numbers are equal given a numeric tolerance
//---------------------------------------------------------------------------------------
template<class T>
typename std::enable_if<!std::numeric_limits<T>::is_integer, bool>::type
inline almost_equal(const T x, const T y, const int ulp =  2)
{
    // the machine epsilon has to be scaled to the magnitude of the values used
    // and multiplied by the desired precision in ULPs (units in the last place)
    return std::fabs(x-y) <= std::numeric_limits<T>::epsilon() * std::fabs(x+y) * ulp
        // unless the result is subnormal
        || std::fabs(x-y) < std::numeric_limits<T>::min();
}
//---------------------------------------------------------------------------------------




//---------------------------------------------------------------------------------------
//Print the content of a container. If the element in the container cannot be printed
//you'll get a huge error during compilation
template<typename T>
inline void print_container(T const& cont) 
{
    for (typename T::const_iterator it = cont.begin(); it != cont.end(); ++it) {
        std::cout << *it << " "; 
    }

    std::cout << std::endl; 
};
//---------------------------------------------------------------------------------------


//---------------------------------------------------------------------------------------
//Sum the element of a container with a custom function. std::partial_sum is terribly limited

//T2 should be a scalar that can be summed
template<typename T1, typename T2>
inline  std::vector<T2> partial_sum_vector(const std::vector<T1>& lst, const std::function<T2(const T1&)>& func){

    //Make sure there are things to add
    assert(lst.size() != 0);

    //Store the sequence of partial sums in this vector
    std::vector<T2> partial_sum_lst{func(lst.front())}; 
    partial_sum_lst.reserve(lst.size());

    //Compute the partial sums
    for (size_t i = 0, j = 1; j < lst.size();  i = j++){
        partial_sum_lst.push_back(partial_sum_lst[i] + func(lst[j]));
    }

    return partial_sum_lst;
};
//---------------------------------------------------------------------------------------




//---------------------------------------------------------------------------------------
//Given a list of indices (to_remove) this function removes the element with the corresponding indices
// in vector. THE INDICES MUST BE SORTED IN ASCENDING ORDER
template<typename T, typename U>
void remove_index(std::vector<T>& vector, std::vector<U>& to_remove){
    auto vector_base = vector.begin();

    typename std::vector<T>::size_type down_by = 0;
    for (auto iter = to_remove.cbegin(); 
                iter < to_remove.cend(); 
                iter++, down_by++){
        typename std::vector<T>::size_type next = iter + 1 == to_remove.cend() ? vector.size() : *(iter + 1);

        std::move(vector_base + *iter + 1, vector_base + next, vector_base + *iter - down_by);
    }
    vector.resize(vector.size() - to_remove.size());
}
//---------------------------------------------------------------------------------------


//---------------------------------------------------------------------------------------
//Flatten a 2D vector
template<typename T>
std::vector<T> flatten(const std::vector<std::vector<T>> &orig)
{   
    std::vector<T> ret;
    for(const auto &v: orig)
        ret.insert(ret.end(), v.begin(), v.end());                                                                                         
    return ret;
}  
//---------------------------------------------------------------------------------------



//---------------------------------------------------------------------------------------
//Format a string of a number into a given format
template<typename T>
inline std::string format_number(T number, std::string fmt){
    char buffer_char[30];
    sprintf(buffer_char, fmt.c_str(), number);
    return std::string(buffer_char);
}
//---------------------------------------------------------------------------------------


//---------------------------------------------------------------------------------------
inline std::string lower_string(std::string& str){
    std::transform(str.begin(), str.end(), str.begin(), [](unsigned char c){ return std::tolower(c); });
    return str;
}
//---------------------------------------------------------------------------------------




//---------------------------------------------------------------------------------------
//This function launches a loop in parallel. If an exception is thrown, the handler stores it 
//and wait for the other threads to finish. Then it rethrows the exception. It prevents the
//program from crashing when an exception is thrown in a parallel loop. Note that the
template<typename T>
inline void parallel_exception_handler(
    const std::vector<T>& vec,             //Contains the vector of objects to iterate over
    const std::function<void(T)>& func     //Note that the function should not modify the size of the vector
) noexcept(false){

    //Create an exception pointer
    std::exception_ptr e_ptr;

    #pragma omp parallel for
    for(size_t i = 0; i < vec.size(); i++){

        try{
            //This might throw an exception
            func(vec[i]);
        }

        //Capture the exception thrown by func
        catch(...){
        
            #pragma omp critical
            {
                //Store the exception
                e_ptr = std::current_exception(); 
            }
        }
    }

    //If an exception was thrown during the execution of the parallel loop, rethrow it
    if (e_ptr) std::rethrow_exception(e_ptr);
}
//---------------------------------------------------------------------------------------



//---------------------------------------------------------------------------------------------------------
inline double cot(const double angle){return 1. / std::tan(angle);}
//---------------------------------------------------------------------------------------------------------








#endif
