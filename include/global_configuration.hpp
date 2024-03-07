#ifndef DEF_GLOBAL_CONFIGURATION
#define DEF_GLOBAL_CONFIGURATION


//This file contains different compilation macros that can be used to activate or deactivate the compilation of 
//some parts of the code. For instance these tags can be used to reduce the attributes of the 
//mesh to the minimum and thus optimize performance. They can also be used to choose the model of polarization
//or the way the equations of motions are solved.


//If this tag is set to true, the face will stored their contact energies (adhesion and repulsion energies), this adds 8 more bytes to the face class
#define FACE_STORE_CONTACT_ENERGY false


//Write the normals of the cell surfaces at the level of the faces and vertices. This adds a lot more data to the mesh files
#define WRITE_NORMALS_IN_OUTPUT_MESH_FILE true


//The faces can be polarize by different techniques. The following index can be used to choose the polarization mode:
//  - 0 : No polarization
//  - 1 : Polarization based on the contacts of the cell
//  - 2 : Polarization based on the discretization of the simulation space
#define POLARIZATION_MODE_INDEX 1

/*
The technique used to compute the contact forces can be chosen via the following index:
    - 0 : The contact forces are computed between pairs of faces and vertices by creating springs between the surfaces
    - 1 : The contact forces are computed mechanically linking pairs of adjacent nodes
    - 2 : The contact forces are computed mechanically linking pairs of adjacent faces

*/
#define CONTACT_MODEL_INDEX 1


/*The dynamic model and the time integration scheme used to run the simulation can be chosen via the following index:
    - 0 : Full equations of motions solved with semi-implicit euler  (1st order integration scheme)
    - 1 : Overdamped equations of motions solved with forward euler  (1st order integration scheme)
*/ 
#define DYNAMIC_MODEL_INDEX 0



#endif