import vtk 
from sys import argv
import numpy as np 
from os import path, _exit
import re

"""
    In some meshes, each face is defined as a vtkCell. This is the case by default in stl, ply and obj meshes. 
    This script regroups all the faces belonging to the same cell and form one single vtkCell from them.
    The produced meesh can then be used to run simulations.
"""




def get_polydata_from_file(filename):
    """
    Read a mesh file and return a vtkPolyData object

    Parameters:
    -----------

    filename: str
        Path to the mesh file

    Returns:
    --------

    polydata: vtkPolyData
        The mesh as a vtkPolyData object

    """

    extension = filename.split(".")[-1]

    if extension == "vtk":
        reader = vtk.vtkUnstructuredGridReader()
    elif extension == "vtp":
        reader = vtk.vtkXMLPolyDataReader()
    elif extension == "stl":
        reader = vtk.vtkSTLReader()
    elif extension == "ply":
        reader = vtk.vtkPLYReader()
    elif extension == "obj":
        reader = vtk.vtkOBJReader()
    elif extension == "vtu":
        reader = vtk.vtkXMLUnstructuredGridReader()
    else: 
        raise ValueError("File extension not recognized")

    reader.SetFileName(filename)
    reader.Update()

    #The output might be at this stage an unstructured grid or a polydata
    reader_output =  reader.GetOutput()

    #We make sure to convert the mesh to polydata
    if isinstance(reader_output, vtk.vtkUnstructuredGrid):
        geom_filter = vtk.vtkGeometryFilter()
        geom_filter.SetInputData(reader_output)
        geom_filter.Update()
        polydata = geom_filter.GetOutput()
        return polydata
    elif isinstance(reader_output, vtk.vtkPolyData):
        return reader_output
    else:
        raise ValueError("The mesh could not be converted to polydata")




def get_cell_ids(polydata):
    """
    Get the cell to which each face and node belongs

    Parameters:
    -----------

    polydata: vtkPolyData
        The mesh as a vtkPolyData object

    Returns:
    --------

    face_cell_id: np.array()
        The id cell to which each face belongs
    
    point_cell_id: np.array()
        The id cell to which each node belongs
    
    """



    #Apply the connectivity filter to the mesh
    connectivity_filter = vtk.vtkConnectivityFilter()
    connectivity_filter.SetInputData(polydata)
    connectivity_filter.SetExtractionModeToAllRegions()
    connectivity_filter.ColorRegionsOn()
    connectivity_filter.Update()
    polydata = connectivity_filter.GetOutput()

    #Get the cell to which each face belongs and same for the nodes
    face_cell_id = np.array(polydata.GetCellData().GetArray("RegionId"))
    point_cell_id = np.array(polydata.GetPointData().GetArray("RegionId"))
    return face_cell_id, point_cell_id





def separate_cell_points_and_faces(polydata, face_cell_id, point_cell_id):
    """
    This method separates in dictionaries the points and faces of each cell. 
    Then this dictionaries can be used to recreate the individual cells

    Parameters:
    -----------

    polydata: vtkPolyData
        The mesh as a vtkPolyData object

    face_cell_id: np.array()
        The id cell to which each face belongs

    point_cell_id: np.array()
        The id cell to which each node belongs

    Returns:
    --------
    face_dic: dict
        A dictionary where the keys are the cell ids and the values are the faces belonging to that cell

    point_dic: dict
        A dictionary where the keys are the cell ids and the values are the points belonging to that cell

    """
    
    #Get the number of cells
    num_cells =  np.unique(face_cell_id).shape[0]  

    #Create a dictionnary where for each cell will be stored in lists all the triangles that
    #belong to that cell, the key is the cell id
    face_dic = {k: [] for k in range(num_cells)}

    #Same than the dictionnary above, here we store all the point IDs of the points belonging 
    #to a given cell
    point_dic = {k: [] for k in range(num_cells)}

    #Loop over the points Ids and store in which cell they belong
    for pointId in range(polydata.GetNumberOfPoints()):
        point_dic[point_cell_id[pointId]].append(pointId)


    #Loop over all the triangles and store them in the cell dic based on their cellIDs
    for faceId in range(polydata.GetNumberOfCells()):

        #Extract the points of the face
        face_points = [polydata.GetCell(faceId).GetPointId(i) for i in range(polydata.GetCell(faceId).GetNumberOfPoints())]
        face_dic[face_cell_id[faceId]].append(face_points)

    return face_dic, point_dic



def write_unstructured_grid(polydata, face_dic, point_dic, output_file_path):
    """
    This method directly writes a file at the legacy unstructured grid vtk format. 

    Parameters:
    -----------

    polydata: vtkPolyData
        The mesh as a vtkPolyData object

    face_dic: dict
        A dictionary where the keys are the cell ids and the values are the faces belonging to that cell

    point_dic: dict
        A dictionary where the keys are the cell ids and the values are the points belonging to that cell

    output_file_path:
        The path where the unstructured grid will be written


    Returns:
    --------

    None
    """

    #Create the file
    f = open(output_file_path, "w")

    #Write the header 
    f.write("# vtk DataFile Version 4.2\nvtk output\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS {} float\n".format(polydata.GetNumberOfPoints()))

    #Loop over the points of the polydata
    vtk_points = polydata.GetPoints()
    
    #Convert the point coordinates to a list
    point_lst = [list(vtk_points.GetPoint(i)) for i in range(vtk_points.GetNumberOfPoints())]

    #Write the points in the list
    for i in range(point_lst.__len__()):
        for point_coord in point_lst[i]:
            f.write("{:.5E} ".format(point_coord))
        if (i + 1) % 3 == 0:  f.write("\n")



    #Compute the number of integer needed to save each cell
    nb_integers = face_dic.__len__() * 2 
    for cell_id in face_dic.keys():
        nb_integers += sum([len(face_dic[cell_id][i]) + 1 for i in range(face_dic[cell_id].__len__())])

    #Write the cell data
    f.write("\nCELLS {} {}\n".format(face_dic.__len__(), nb_integers))

    #Loop over the cells
    for cell_id in face_dic.keys():

        #Get the number of faces of the cell
        nb_faces = face_dic[cell_id].__len__()

        #Get the number of integer needed to save the cell
        nb_integers_cell = sum([len(face_dic[cell_id][i]) + 1 for i in range(face_dic[cell_id].__len__())]) + 1

        f.write("{} {} ".format(nb_integers_cell, nb_faces))

        #Loop over the faces of the cell
        for face in face_dic[cell_id]:

            #Write the number of points of the face
            f.write("{} ".format(face.__len__()))

            #Write the points of the face
            for point in face:
                f.write("{} ".format(point))

        #Write a new line
        f.write("\n")
 

    #Indicate the type of cell
    f.write("\nCELL_TYPES {}\n".format(face_dic.__len__()))
    for cell_id in face_dic.keys():
        f.write("42\n")


if __name__ == "__main__":

    assert argv.__len__() == 3, "Wrong command line argument. Example of correct usage\n"+\
    "python convert_to_cell_mesh.py /path/to/mesh/file /path/to/output/file"

    #The absoltue path the simulation folder
    path_to_input_mesh_file = argv[1].strip()
    path_to_output_mesh_file = argv[2].strip()

    assert path.exists(path_to_input_mesh_file), "The path to the input mesh file does not exist."

    #Read the mesh
    polydata = get_polydata_from_file(path_to_input_mesh_file)

    #Get the cell to which each face belongs
    face_cell_id, point_cell_id = get_cell_ids(polydata)

    #Separate the faces and points of each cell
    face_dic, point_dic = separate_cell_points_and_faces(polydata, face_cell_id, point_cell_id)


    #Get the cells unstructured grids
    write_unstructured_grid(polydata, face_dic, point_dic, path_to_output_mesh_file)