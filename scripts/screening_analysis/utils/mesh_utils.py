import os
import trimesh as tm
from trimesh.voxel import creation
import vtk
from typing import List, Union, Dict

"""
Utility function for mesh loading, extraction, conversion and preparation.
"""

#---------------------------------------------------------------------------------------------------------------------------------------------------------
def convert_vtk_cell_to_trimesh(
        vtk_cell: vtk.vtkCell
) -> tm.Trimesh:
    """
    Convert a vtk cell to a trimesh cell

    Parameters:
    -----------
    vtk_cell: vtk.vtkCell
        The vtk cell to convert

    Returns:
    --------
    trimesh_cell: tm.Trimesh
        The trimesh cell
    """

    #Get the point id of all the points of the cell
    global_point_id_lst = [vtk_cell.GetPointId(i) for i in range(vtk_cell.GetNumberOfPoints())]

    #Get the coordinates of all the points of the cell
    point_coordinate_lst = [vtk_cell.GetPoints().GetPoint(i) for i in range(vtk_cell.GetNumberOfPoints())]

    #Create a map between the global point id and the local point ids
    global_to_local_point_id_map = {global_point_id_lst[i] : i for i in range(len(global_point_id_lst))}

    faces_lst = []

    #Loop over the faces of the cell
    for i in range(vtk_cell.GetNumberOfFaces()):
            
            #Get the face
            face = vtk_cell.GetFace(i)
    
            #Get the number of points of the face
            num_points = face.GetNumberOfPoints()
    
            #Create a list to store the points of the face
            face_points_lst = []
    
            #Loop over the points of the face
            for j in range(num_points):
    
                #Get the global id of the point
                global_point_id = face.GetPointId(j)

                #Convert the global id to the local id
                local_point_id = global_to_local_point_id_map[global_point_id]
    
                #Append the point to the list
                face_points_lst.append(local_point_id)
    
            #Append the face to the list
            faces_lst.append(face_points_lst)


    #Create the trimesh cell
    trimesh_cell = tm.Trimesh(point_coordinate_lst, faces_lst)

    return trimesh_cell
#---------------------------------------------------------------------------------------------------------------------------------------------------------



#---------------------------------------------------------------------------------------------------------------------------------------------------------
def get_cell_geometries_from_file(
        path_to_input_mesh: str,
        cached_meshes: Union[None, Dict[int, List[tm.Trimesh]]]
    ) -> List[tm.Trimesh]:
    """
    Load the initial mesh file and return each cell as a trimesh object in a list

    Parameters:
    -----------
        path_to_input_mesh: (str)
            The path to the input mesh file

        cached_meshes: (Union[None, Dict[int, List[tm.Trimesh]]])
            A dictionary in which keys are indexes associated to already loaded meshes and
            values are lists of cell meshes in trimesh format.
            If `None`, meshes lists are not cached.

    Returns:
    --------
        cell_lst: (List[tm.Trimesh])
            The list of the cells meshes in trimesh format
    """

    #get number of file
    file_name = os.path.basename(path_to_input_mesh)
    file_id = file_name.split("_")[-1][:-4]

    if isinstance(cached_meshes, dict) and file_id in cached_meshes:
        return cached_meshes[file_id]
    
    reader = vtk.vtkUnstructuredGridReader()
    reader.SetFileName(path_to_input_mesh)
    reader.Update()
    ugrid = reader.GetOutput()

    cell_mesh_lst =  []
    for i in range(ugrid.GetNumberOfCells()):
        vtk_cell = ugrid.GetCell(i)
        trimesh_cell_mesh = convert_vtk_cell_to_trimesh(vtk_cell)
        cell_mesh_lst.append(trimesh_cell_mesh)
    
    if isinstance(cached_meshes, dict):
        cached_meshes[file_id] = cell_mesh_lst

    return cell_mesh_lst
#---------------------------------------------------------------------------------------------------------------------------------------------------------
