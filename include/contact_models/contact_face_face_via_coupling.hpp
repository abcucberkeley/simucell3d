#ifndef DEF_CONTACT_FACE_FACE_VIA_COUPLING 
#define DEF_CONTACT_FACE_FACE_VIA_COUPLING


#include "global_configuration.hpp"

#define _USE_MATH_DEFINES

#include "contact_model_abstract.hpp"
#include <math.h>


class contact_face_face_via_coupling: public contact_model_abstract{

    private:
        static constexpr double max_dot_product_adhesion_  = std::cos(90 * M_PI / 180.0); 
        static constexpr double max_dot_product_repulsion_ = std::cos(90 * M_PI / 180.0); 


    public:
        contact_face_face_via_coupling() = default;                                          //default constructor
        contact_face_face_via_coupling(const contact_face_face_via_coupling& v) = delete;           //copy constructor
        contact_face_face_via_coupling(contact_face_face_via_coupling&& v) = delete;                //move constructor
        contact_face_face_via_coupling& operator=(const contact_face_face_via_coupling& v) = delete;//copy assignment operator
        contact_face_face_via_coupling& operator=(contact_face_face_via_coupling&& v) = default;     //move assignment operator 

        //The trivial constructor of the contact model
        contact_face_face_via_coupling(const global_simulation_parameters& sim_parameters) noexcept(false);
            
        void run(const std::vector<cell_ptr>& cell_lst) noexcept override;

        //Find the faces that are within a distance below the contact cutoff and apply adhesive or repulsive forces 
        void resolve_all_contacts(const std::vector<cell_ptr>& cell_lst) noexcept;


        //Prevent the surfaces of the cells to interpenetrate by applying repulsive forces on the surfaces.
        //These forces are only applied if the 2 cells are not epithelial
        void resolve_contact(cell_ptr c1, cell_ptr c2, node& n, face* f) const noexcept;


};
#endif
