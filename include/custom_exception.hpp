#ifndef DEF_CUSTOM_EXCEPTION
#define DEF_CUSTOM_EXCEPTION


#include <exception>
#include <string>




//---------------------------------------------------------------------------------------
//Throw this exception when there is a problem with the automatic polarization
class automatic_polarization_exception: public std::exception
{
    private:
        std::string error_msg_;

public:
    automatic_polarization_exception(const std::string& error_msg) throw() :error_msg_(error_msg) {}

    virtual const char* what() const throw() {return error_msg_.c_str();}
     
    //Destructor is probably not needed
    virtual ~automatic_polarization_exception() throw(){}
};

//---------------------------------------------------------------------------------------






//---------------------------------------------------------------------------------------
//Throw this exception when the solver starts to produce NaNs
class unstable_simulation_exception: public std::exception
{
    private:
        std::string error_msg_;

public:
    unstable_simulation_exception(const std::string& error_msg) throw() :error_msg_(error_msg) {}

    virtual const char* what() const throw() {return error_msg_.c_str();}
     
    //Destructor is probably not needed
    virtual ~unstable_simulation_exception() throw(){}
};

//---------------------------------------------------------------------------------------




//---------------------------------------------------------------------------------------
//Exception related to the parameter xml file reader
class intialization_exception: public std::exception
{
    private:
        std::string error_msg_;

public:
    intialization_exception(const std::string& error_msg) throw() :error_msg_(error_msg) {}

    virtual const char* what() const throw() {return error_msg_.c_str();}
     
    //Destructor is probably not needed
    virtual ~intialization_exception() throw(){}
};

//---------------------------------------------------------------------------------------


//---------------------------------------------------------------------------------------
//Throw this exception if something goes wrong during the division of cell
class division_exception: public std::exception
{
    private:
        std::string error_msg_;

public:
    division_exception(const std::string& error_msg) throw() :error_msg_(error_msg) {}

    virtual const char* what() const throw() {return error_msg_.c_str();}
     
    //Destructor is probably not needed
    virtual ~division_exception() throw(){}
};

//---------------------------------------------------------------------------------------




//---------------------------------------------------------------------------------------
//Exception related to the parameter xml file reader
class parameter_reader_exception: public std::exception
{
    private:
        std::string error_msg_;

public:
    parameter_reader_exception(const std::string& error_msg) throw() :error_msg_(error_msg) {}

    virtual const char* what() const throw() {return error_msg_.c_str();}
     
    //Destructor is probably not needed
    virtual ~parameter_reader_exception() throw(){}
};

//---------------------------------------------------------------------------------------


//---------------------------------------------------------------------------------------
//Exception related to the ball pivoting algorithm
class bpa_exception: public std::exception
{
    private:
        std::string error_msg_;

public:
    bpa_exception(const std::string& error_msg) throw() :error_msg_(error_msg) {}

    virtual const char* what() const throw() {return error_msg_.c_str();}
     
    //Destructor is probably not needed
    virtual ~bpa_exception() throw(){}
};
//---------------------------------------------------------------------------------------




//---------------------------------------------------------------------------------------
//When the mesh of a cell is not manifold (holes in the mesh or edges connected to more than 2 face)
class mesh_integrity_exception: public std::exception
{
    private:
        std::string error_msg_;

public:
    mesh_integrity_exception(const std::string& error_msg) throw() :error_msg_(error_msg) {}

    virtual const char* what() const throw() {return error_msg_.c_str();}
     
    //Destructor is probably not needed
    virtual ~mesh_integrity_exception() throw(){}
};
//---------------------------------------------------------------------------------------


//---------------------------------------------------------------------------------------
//Exception related to a problem in the input geometry
class initial_triangulation_exception: public std::exception
{
    private:
        std::string error_msg_;

public:
    initial_triangulation_exception(const std::string& error_msg) throw() :error_msg_(error_msg) {}

    virtual const char* what() const throw() {return error_msg_.c_str();}
     
    //Destructor is probably not needed
    virtual ~initial_triangulation_exception() throw(){}
};
//---------------------------------------------------------------------------------------



//---------------------------------------------------------------------------------------
//Exception indicating a bad formatting of the input mesh file
class mesh_reader_exception: public std::exception
{
    private:
        std::string error_msg_;

public:
    mesh_reader_exception(const std::string& error_msg) throw() :error_msg_(error_msg) {}

    virtual const char* what() const throw() {return error_msg_.c_str();}
     
    //Destructor is probably not needed
    virtual ~mesh_reader_exception() throw(){}
};
//---------------------------------------------------------------------------------------


//---------------------------------------------------------------------------------------
//Exception related to a problem during the writing of a mesh file
class mesh_writer_exception: public std::exception
{
    private:
        std::string error_msg_;

public:
    mesh_writer_exception(const std::string& error_msg) throw() :error_msg_(error_msg) {}

    virtual const char* what() const throw() {return error_msg_.c_str();}
     
    //Destructor is probably not needed
    virtual ~mesh_writer_exception() throw(){}
};
//---------------------------------------------------------------------------------------




//---------------------------------------------------------------------------------------
//Exception related to a problem during the writing of a mesh file
class uspg_exception: public std::exception
{
    private:
        std::string error_msg_;

public:
    uspg_exception(const std::string& error_msg) throw() :error_msg_(error_msg) {}

    virtual const char* what() const throw() {return error_msg_.c_str();}
     
    //Destructor is probably not needed
    virtual ~uspg_exception() throw(){}
};
//---------------------------------------------------------------------------------------








#endif