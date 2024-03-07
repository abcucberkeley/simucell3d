import numpy as np 
import vtk 
from os import path, _exit, getcwd, listdir 
from scipy import optimize


"""
This script calculates the angle between the cells of a cell doublet.
"""


#------------------------------------------------------------------------------------------------------------------
def read_mesh(path_to_mesh:str):
    """
    Reads the mesh from the path given as parameter.

    Parameters
    ----------
    path_to_mesh : str
        Path to the mesh.


    Returns
    -------
    ugrid : vtkUnstructuredGrid
    """

    reader = vtk.vtkUnstructuredGridReader()
    reader.SetFileName(path_to_mesh)
    reader.Update()
    ugrid = reader.GetOutput()
    return ugrid
#------------------------------------------------------------------------------------------------------------------



#------------------------------------------------------------------------------------------------------------------
def get_cell_apical_face_lst(cell_id:int, mesh:vtk.vtkUnstructuredGrid()):
    """
    
    Get all the apical faces belonging to a certain cell
    
    Parameters:
    -----------
    cell_id : int
        Id of the cell.

    mesh : vtk.vtkUnstructuredGrid
        Mesh containing the cells.

    Returns:
    --------
    apical_face_lst : list
        List containing the apical faces belonging to the cell.

        
    """

    
    #Use a threshold filter to get all the faces of the cell
    threshold_filter_1 = vtk.vtkThreshold()
    threshold_filter_1.SetInputData(mesh)
    threshold_filter_1.SetInputArrayToProcess(0, 0, 0, vtk.vtkDataObject.FIELD_ASSOCIATION_CELLS, "face_cell_id")
    threshold_filter_1.ThresholdBetween(cell_id, cell_id)
    threshold_filter_1.Update()


    #Use a second thresold to get the apical faces of the cell
    threshold_filter_2 = vtk.vtkThreshold()
    threshold_filter_2.SetInputData(threshold_filter_1.GetOutput())
    threshold_filter_2.SetInputArrayToProcess(0, 0, 0, vtk.vtkDataObject.FIELD_ASSOCIATION_CELLS, "face_type_id")
    threshold_filter_2.ThresholdBetween(0, 0)
    threshold_filter_2.Update()

    return threshold_filter_2.GetOutput()
#------------------------------------------------------------------------------------------------------------------



#------------------------------------------------------------------------------------------------------------------
def slice_faces(mesh : vtk.vtkUnstructuredGrid()):
    """Slice the faces along the z axis and return the slice.
    
    Parameters:
    -----------
    mesh : vtk.vtkUnstructuredGrid
        Mesh containing the faces.

    Returns:
    --------
    slice_mesh : vtk.vtkPolyData
        Slice of the faces.
    """

    #Create the vtk plane that will cut the cells
    vtk_plane = vtk.vtkPlane()
    vtk_plane.SetOrigin(mesh.GetCenter())
    vtk_plane.SetNormal([0, 1, 0])

    slice_filter = vtk.vtkPlaneCutter()
    slice_filter.SetInputData(mesh)
    slice_filter.SetPlane(vtk_plane)
    slice_filter.Update()

    return slice_filter.GetOutput()
#------------------------------------------------------------------------------------------------------------------




#------------------------------------------------------------------------------------------------------------------

def calc_distance_between_centroid_and_points(centroid_ar, points_ar):
    """
    Compute the distances of the points in points_lst to the centroid (centroid_x, centroid_y)
    
    Parameters:
    -----------
    centroid_ar : np.ndarray shape (2,)
        Centroid of the points.

    points_lst : np.ndarray shape (n, 2)
        List of points.

    Returns:
    --------
    distances_lst : np.ndarray shape (n,)
        List of distances.
    """
    return np.linalg.norm(points_ar - centroid_ar, axis=1)



def loss_function(centroid_ar, points_ar):
    """
    Compute the distances of the points in points_lst to the centroid (centroid_x, centroid_y)
    and then substract the mean of the distances.
    
    Parameters:
    -----------
    centroid_ar : np.ndarray shape (2,)
        Centroid of the points.

    points_lst : np.ndarray shape (n, 2)
        List of points.

    Returns:
    --------
    distances_lst : np.ndarray shape (n,)
        List of updated distances.
    """
    distance_ar = calc_distance_between_centroid_and_points(centroid_ar, points_ar)
    return distance_ar - distance_ar.mean()




def jacobian_function(centroid_ar, points_ar):
    """
    This method computes the jacobian of the loss function.
    The axis corresponding to derivatives must be coherent with the col_deriv option of leastsq"""
    xc, yc     = centroid_ar.tolist()
    df2b_dc    = np.empty((2, points_ar.shape[0]))

    Ri = calc_distance_between_centroid_and_points(centroid_ar, points_ar)
    df2b_dc[0] = (xc - points_ar[:,0])/Ri  # dR/dxc
    df2b_dc[1] = (yc - points_ar[:,1])/Ri  # dR/dyc
    df2b_dc    = df2b_dc - df2b_dc.mean(axis=1)[:, np.newaxis]

    return df2b_dc





def fit_circle_to_points(point_ar):
    """
    Fit a circle to the points in points_lst.

    Parameters:
    -----------
    points_lst : np.ndarray shape (n, 2)
        List of points.

    Returns:
    --------
    centroid_ar : np.ndarray shape (2,)
        Centroid of the points.

    radius : float
        Radius of the circle.
    """
    
    #Initial guess for the centroid and radius
    centroid_ar = point_ar.mean(axis=0)

    #Compute the centroid and radius that minimize the loss function
    new_centroid_ar, output = optimize.leastsq(loss_function, centroid_ar, args=(point_ar), Dfun=jacobian_function, col_deriv=True)

    #Compute the radius
    radius = calc_distance_between_centroid_and_points(new_centroid_ar, point_ar).mean()
    return new_centroid_ar, radius

#------------------------------------------------------------------------------------------------------------------




#------------------------------------------------------------------------------------------------------------------
def save_circle(centroid_ar, radius, filename:str):
    """
    Save the circle defined by the centroid and radius to a file.

    Parameters:
    -----------
    centroid_ar : np.ndarray shape (2,)
        Centroid of the points.

    radius : float
        Radius of the circle.

    filename : str
        Name of the file to save the circle to.
    """
    
    #Create the circle
    circle = vtk.vtkRegularPolygonSource()
    circle.SetCenter(centroid_ar[0], centroid_ar[1], centroid_ar[2])
    circle.SetRadius(radius)
    circle.SetNormal(0, 1, 0)
    circle.SetNumberOfSides(50)
    circle.Update()

    #Save the circle to a file
    vtk_writer = vtk.vtkXMLPolyDataWriter()
    vtk_writer.SetFileName(filename)
    vtk_writer.SetInputData(circle.GetOutput())
    vtk_writer.Write()

#------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------
def calc_angle_between_cells(path_to_mesh:str, index:int):
    """
    Calculate the angle between the cells of a cell doublet.

    Parameters:
    -----------
    path_to_mesh : str
        Path to the mesh containing the cell doublet geometry.
    
    Returns:
    --------
    theta : float
        Angle between the cells of the cell doublet.
    """

    #Read the mesh
    ugrid_mesh = read_mesh(path_to_mesh)

    #Get the apical faces of the cells
    cell_0_apical_faces_mesh = get_cell_apical_face_lst(0, ugrid_mesh)
    cell_1_apical_faces_mesh = get_cell_apical_face_lst(1, ugrid_mesh)

    #Slice the apical faces along the y-axis
    cell_0_apical_faces_slice_mesh = slice_faces(cell_0_apical_faces_mesh)
    cell_1_apical_faces_slice_mesh = slice_faces(cell_1_apical_faces_mesh)

    #Save the slices to a file
    vtk_writer_0 = vtk.vtkXMLPolyDataWriter()
    vtk_writer_0.SetFileName(f"cell_0_apical_faces_slice_{index}.vtp")
    vtk_writer_0.SetInputData(cell_0_apical_faces_slice_mesh)
    vtk_writer_0.Write()

    vtk_writer_1 = vtk.vtkXMLPolyDataWriter()
    vtk_writer_1.SetFileName(f"cell_1_apical_faces_slice_{index}.vtp")
    vtk_writer_1.SetInputData(cell_1_apical_faces_slice_mesh)
    vtk_writer_1.Write()

    #Extract the points coordinates from the slices
    cell_0_apical_faces_slice_points = cell_0_apical_faces_slice_mesh.GetPoints()
    cell_1_apical_faces_slice_points = cell_1_apical_faces_slice_mesh.GetPoints()

    #Get the points coordinates as numpy arrays
    points_cell_1_ar = np.array([cell_0_apical_faces_slice_points.GetPoint(i) for i in range(cell_0_apical_faces_slice_points.GetNumberOfPoints())])
    points_cell_2_ar = np.array([cell_1_apical_faces_slice_points.GetPoint(i) for i in range(cell_1_apical_faces_slice_points.GetNumberOfPoints())])

    #Delete the frist and last 5 points
    points_cell_1_ar = np.delete(points_cell_1_ar, [0, 1, 2, 3, 4, -1, -2, -3, -4, -5], axis=0)
    points_cell_2_ar = np.delete(points_cell_2_ar, [0, 1, 2, 3, 4, -1, -2, -3, -4, -5], axis=0)

    centroid_3d_cell_1_ar = np.mean(points_cell_1_ar, axis=0)
    centroid_3d_cell_2_ar = np.mean(points_cell_2_ar, axis=0)

    #Remove the y coordinate of the points
    points_2d_cell_1_ar = np.delete(points_cell_1_ar, 1, axis=1)
    points_2d_cell_2_ar = np.delete(points_cell_2_ar, 1, axis=1)

    #Fit circles to the points
    centroid_2d_cell_1_ar, radius_1 = fit_circle_to_points(points_2d_cell_1_ar)
    centroid_2d_cell_2_ar, radius_2 = fit_circle_to_points(points_2d_cell_2_ar)

    #Reconstitude the 3d centroid
    centroid_3d_cell_1_ar = np.array([centroid_2d_cell_1_ar[0], centroid_3d_cell_1_ar[1], centroid_2d_cell_1_ar[1]])
    centroid_3d_cell_2_ar = np.array([centroid_2d_cell_2_ar[0], centroid_3d_cell_2_ar[1], centroid_2d_cell_2_ar[1]])

    save_circle(centroid_3d_cell_1_ar, radius_1, filename=f"cell_1_circle_{index}.vtp")
    save_circle(centroid_3d_cell_2_ar, radius_2, filename=f"cell_2_circle_{index}.vtp")


    d = np.linalg.norm(centroid_3d_cell_1_ar - centroid_3d_cell_2_ar)
    if np.abs(radius_1 - radius_2) < d and d < radius_1 + radius_2:    
        l = np.sqrt((np.power(radius_1+radius_2, 2) - d*d)*(d*d - np.power(radius_1-radius_2, 2))) / d
        theta = np.arcsin(l / (2. * radius_1)) + np.arcsin(l / (2. * radius_2))
        return theta / 2.   
    return float("nan")






if __name__ == "__main__":

    #Get the path to the directory where this script is stored 
    path_to_script_folder = path.dirname(path.realpath(__file__))

    #Get the path to the simucell3d python library
    path_to_build = path.abspath(path.join(path_to_script_folder, "..", "..", "..", "..", "build"))

    #Path to the directory where the meshes are stored
    path_to_mesh_folder = path.join(path_to_build, "cell_doublet_angle_output", "simulation_1", "face_data")

    #List the files in this folder
    mesh_filename_lst = listdir(path_to_mesh_folder)

    #Sort the files
    mesh_filename_lst = sorted(mesh_filename_lst, key=lambda x: int(x.split("_")[1].split(".")[0]))

    #Get the last file
    path_to_mesh = path.join(path_to_mesh_folder, mesh_filename_lst[-1])

    theta = calc_angle_between_cells(path_to_mesh, 0)
    print("\nangle between cells:", theta, "\n")





    #Path to the mesh containing the cell doublet geometry
    #path_to_mesh = path.join(getcwd(), "test_geometry.vtk")



