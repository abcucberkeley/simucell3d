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
    
    Get all the latreal faces belonging to a certain cell
    
    Parameters:
    -----------
    cell_id : int
        Id of the cell.

    mesh : vtk.vtkUnstructuredGrid
        Mesh containing the cells.

    Returns:
    --------
    apical_face_lst : list
        List containing the latreal faces belonging to the cell.

        
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
def get_intersection_points_between_circles(centroid_1_ar, radius_1, centroid_2_ar, radius_2):
    """
    Return the 2 intersection points between two circles.

    Parameters:
    -----------

    centroid_1_ar : np.ndarray shape (3,)
        Centroid of the first circle.

    radius_1 : float
        Radius of the first circle.

    centroid_2_ar : np.ndarray shape (3,)
        Centroid of the second circle.

    radius_2 : float
        Radius of the second circle.

    Returns:
    --------
    intersection_points_lst : list
        List containing the 2 intersection points. in 3D space.
    """

    #Remove the y coordinate of the centroids
    centroid_1_2d_ar = np.delete(centroid_1_ar.copy(), 1)
    centroid_2_2d_ar = np.delete(centroid_2_ar.copy(), 1)

    d = np.linalg.norm(np.array(centroid_1_2d_ar) - np.array(centroid_2_2d_ar))
    
    # Check for no intersection or one circle being completely within the other
    if d > radius_1 + radius_2 or d < abs(radius_1 - radius_2):
        return None
    
    # Calculate the intersection points if the circles intersect
    a = (radius_1**2 - radius_2**2 + d**2) / (2 * d)
    h = np.sqrt(radius_1**2 - a**2)
    
    x2 = centroid_1_2d_ar[0] + (a * (centroid_2_2d_ar[0] - centroid_1_2d_ar[0])) / d
    y2 = centroid_1_2d_ar[1] + (a * (centroid_2_2d_ar[1] - centroid_1_2d_ar[1])) / d
    
    intersection1 = [x2 + h * (centroid_2_2d_ar[1] - centroid_1_2d_ar[1]) / d, y2 - h * (centroid_2_2d_ar[0] - centroid_1_2d_ar[0]) / d]
    intersection2 = [x2 - h * (centroid_2_2d_ar[1] - centroid_1_2d_ar[1]) / d, y2 + h * (centroid_2_2d_ar[0] - centroid_1_2d_ar[0]) / d]
    
    #Add the y coordinate to the intersection points
    intersection_point_1 = np.array([intersection1[0], centroid_1_ar[1], intersection1[1]])
    intersection_point_2 = np.array([intersection2[0], centroid_1_ar[1], intersection2[1]])

    return [intersection_point_1, intersection_point_2]

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
    cell_2_apical_faces_mesh = get_cell_apical_face_lst(2, ugrid_mesh)


    #Slice the apical faces along the y-axis
    cell_0_apical_faces_slice_mesh = slice_faces(cell_0_apical_faces_mesh)
    cell_1_apical_faces_slice_mesh = slice_faces(cell_1_apical_faces_mesh)
    cell_2_apical_faces_slice_mesh = slice_faces(cell_2_apical_faces_mesh)


    #Save the slices to a file
    vtk_writer_0 = vtk.vtkXMLPolyDataWriter()
    vtk_writer_0.SetFileName(f"cell_0_apical_faces_slice_{index}.vtp")
    vtk_writer_0.SetInputData(cell_0_apical_faces_slice_mesh)
    vtk_writer_0.Write()

    vtk_writer_1 = vtk.vtkXMLPolyDataWriter()
    vtk_writer_1.SetFileName(f"cell_1_apical_faces_slice_{index}.vtp")
    vtk_writer_1.SetInputData(cell_1_apical_faces_slice_mesh)
    vtk_writer_1.Write()

    vtk_writer_2 = vtk.vtkXMLPolyDataWriter()
    vtk_writer_2.SetFileName(f"cell_2_apical_faces_slice_{index}.vtp")
    vtk_writer_2.SetInputData(cell_2_apical_faces_slice_mesh)
    vtk_writer_2.Write()
  
    #Extract the points coordinates from the slices
    cell_0_apical_faces_slice_points = cell_0_apical_faces_slice_mesh.GetPoints()
    cell_1_apical_faces_slice_points = cell_1_apical_faces_slice_mesh.GetPoints()
    cell_2_apical_faces_slice_points = cell_2_apical_faces_slice_mesh.GetPoints()


    #Get the points coordinates as numpy arrays
    points_cell_1_ar = np.array([cell_0_apical_faces_slice_points.GetPoint(i) for i in range(cell_0_apical_faces_slice_points.GetNumberOfPoints())])
    points_cell_2_ar = np.array([cell_1_apical_faces_slice_points.GetPoint(i) for i in range(cell_1_apical_faces_slice_points.GetNumberOfPoints())])
    points_cell_3_ar = np.array([cell_2_apical_faces_slice_points.GetPoint(i) for i in range(cell_2_apical_faces_slice_points.GetNumberOfPoints())])


    #Delete the frist and last 5 points
    points_cell_1_ar = np.delete(points_cell_1_ar, [0, 1, 2, 3, 4, -1, -2, -3, -4, -5], axis=0)
    points_cell_2_ar = np.delete(points_cell_2_ar, [0, 1, 2, 3, 4, -1, -2, -3, -4, -5], axis=0)
    points_cell_3_ar = np.delete(points_cell_3_ar, [0, 1, 2, 3, 4, -1, -2, -3, -4, -5], axis=0)


    centroid_3d_cell_1_ar = np.mean(points_cell_1_ar, axis=0)
    centroid_3d_cell_2_ar = np.mean(points_cell_2_ar, axis=0)
    centroid_3d_cell_3_ar = np.mean(points_cell_2_ar, axis=0)

    #Remove the y coordinate of the points
    points_2d_cell_1_ar = np.delete(points_cell_1_ar, 1, axis=1)
    points_2d_cell_2_ar = np.delete(points_cell_2_ar, 1, axis=1)
    points_2d_cell_3_ar = np.delete(points_cell_3_ar, 1, axis=1)


    #Fit circles to the points
    centroid_2d_cell_1_ar, radius_1 = fit_circle_to_points(points_2d_cell_1_ar)
    centroid_2d_cell_2_ar, radius_2 = fit_circle_to_points(points_2d_cell_2_ar)
    centroid_2d_cell_3_ar, radius_3 = fit_circle_to_points(points_2d_cell_3_ar)

    #Reconstitude the 3d centroid
    centroid_3d_cell_1_ar = np.array([centroid_2d_cell_1_ar[0], centroid_3d_cell_1_ar[1], centroid_2d_cell_1_ar[1]])
    centroid_3d_cell_2_ar = np.array([centroid_2d_cell_2_ar[0], centroid_3d_cell_2_ar[1], centroid_2d_cell_2_ar[1]])
    centroid_3d_cell_3_ar = np.array([centroid_2d_cell_3_ar[0], centroid_3d_cell_2_ar[1], centroid_2d_cell_3_ar[1]])

    save_circle(centroid_3d_cell_1_ar, radius_1, filename=f"cell_1_circle_{index}.vtp")
    save_circle(centroid_3d_cell_2_ar, radius_2, filename=f"cell_2_circle_{index}.vtp")
    save_circle(centroid_3d_cell_3_ar, radius_3, filename=f"cell_3_circle_{index}.vtp")

    #Get the intersection points between the circles
    cell_1_and_2_intersection_points = get_intersection_points_between_circles(centroid_3d_cell_1_ar, radius_1, centroid_3d_cell_2_ar, radius_2)
    cell_2_and_3_intersection_points = get_intersection_points_between_circles(centroid_3d_cell_2_ar, radius_2, centroid_3d_cell_3_ar, radius_3)
    cell_3_and_1_intersection_points = get_intersection_points_between_circles(centroid_3d_cell_3_ar, radius_3, centroid_3d_cell_1_ar, radius_1)

    #Save the intersection points to a vtk mesh file 
    cell_1_and_2_intersection_mesh = vtk.vtkPolyData()
    cell_2_and_3_intersection_mesh = vtk.vtkPolyData()
    cell_3_and_1_intersection_mesh = vtk.vtkPolyData()

    cell_1_and_2_intersection_mesh.SetPoints(vtk.vtkPoints())
    cell_2_and_3_intersection_mesh.SetPoints(vtk.vtkPoints())
    cell_3_and_1_intersection_mesh.SetPoints(vtk.vtkPoints())

    for intersection_point in cell_1_and_2_intersection_points:
        cell_1_and_2_intersection_mesh.GetPoints().InsertNextPoint(intersection_point[0], intersection_point[1], intersection_point[2])

    for intersection_point in cell_2_and_3_intersection_points:
        cell_2_and_3_intersection_mesh.GetPoints().InsertNextPoint(intersection_point[0], intersection_point[1], intersection_point[2])

    for intersection_point in cell_3_and_1_intersection_points:
        cell_3_and_1_intersection_mesh.GetPoints().InsertNextPoint(intersection_point[0], intersection_point[1], intersection_point[2])

    vtk_writer_1 = vtk.vtkXMLPolyDataWriter()
    vtk_writer_1.SetFileName(f"cell_1_and_2_intersection_points_{index}.vtp")
    vtk_writer_1.SetInputData(cell_1_and_2_intersection_mesh)
    vtk_writer_1.Write()

    vtk_writer_2 = vtk.vtkXMLPolyDataWriter()
    vtk_writer_2.SetFileName(f"cell_2_and_3_intersection_points_{index}.vtp")
    vtk_writer_2.SetInputData(cell_2_and_3_intersection_mesh)
    vtk_writer_2.Write()

    vtk_writer_3 = vtk.vtkXMLPolyDataWriter()
    vtk_writer_3.SetFileName(f"cell_3_and_1_intersection_points_{index}.vtp")
    vtk_writer_3.SetInputData(cell_3_and_1_intersection_mesh)
    vtk_writer_3.Write()

    #Remove the y coordinate of the intersection points
    cell_1_and_2_intersection_points_2d_ar = np.delete(cell_1_and_2_intersection_points.copy(), 1, axis=1)
    cell_2_and_3_intersection_points_2d_ar = np.delete(cell_2_and_3_intersection_points.copy(), 1, axis=1)
    cell_3_and_1_intersection_points_2d_ar = np.delete(cell_3_and_1_intersection_points.copy(), 1, axis=1)

    #Create the vectors between the intersection points 
    vector_1_2 = cell_1_and_2_intersection_points_2d_ar[1] - cell_1_and_2_intersection_points_2d_ar[0]
    vector_2_3 = cell_2_and_3_intersection_points_2d_ar[1] - cell_2_and_3_intersection_points_2d_ar[0]
    vector_3_1 = cell_3_and_1_intersection_points_2d_ar[1] - cell_3_and_1_intersection_points_2d_ar[0]

    #Save the 3 vectors to a mesh file, the origin of the vectors is the first intersection point
    vector_1_2_mesh = vtk.vtkPolyData()
    vector_2_3_mesh = vtk.vtkPolyData()
    vector_3_1_mesh = vtk.vtkPolyData()

    vector_1_2_mesh.SetPoints(vtk.vtkPoints())
    vector_2_3_mesh.SetPoints(vtk.vtkPoints())
    vector_3_1_mesh.SetPoints(vtk.vtkPoints())

    vector_1_2_mesh.GetPoints().InsertNextPoint(cell_1_and_2_intersection_points_2d_ar[0,0], cell_1_and_2_intersection_points[0][1], cell_1_and_2_intersection_points_2d_ar[0,1])
    vector_2_3_mesh.GetPoints().InsertNextPoint(cell_2_and_3_intersection_points_2d_ar[0,0], cell_2_and_3_intersection_points[0][1], cell_2_and_3_intersection_points_2d_ar[0,1])
    vector_3_1_mesh.GetPoints().InsertNextPoint(cell_3_and_1_intersection_points_2d_ar[0,0], cell_3_and_1_intersection_points[0][1], cell_3_and_1_intersection_points_2d_ar[0,1])

    # Create a vtk.vtkDoubleArray to store the vector
    vector_array_1_2 = vtk.vtkDoubleArray()
    vector_array_2_3 = vtk.vtkDoubleArray()
    vector_array_3_1 = vtk.vtkDoubleArray()

    vector_array_1_2.SetNumberOfComponents(3)  # Assuming a 3D vector
    vector_array_2_3.SetNumberOfComponents(3)  # Assuming a 3D vector
    vector_array_3_1.SetNumberOfComponents(3)  # Assuming a 3D vector

    # Set the vector components in the array
    vector_array_1_2.InsertNextTuple3(vector_1_2[0], 0, vector_1_2[1])
    vector_array_2_3.InsertNextTuple3(vector_2_3[0], 0, vector_2_3[1])
    vector_array_3_1.InsertNextTuple3(vector_3_1[0], 0, vector_3_1[1])

    # Add the vector array to the point data
    vector_1_2_mesh.GetPointData().SetVectors(vector_array_1_2)
    vector_2_3_mesh.GetPointData().SetVectors(vector_array_2_3)
    vector_3_1_mesh.GetPointData().SetVectors(vector_array_3_1)

    #Write the mesh
    vtk_writer_4 = vtk.vtkXMLPolyDataWriter()
    vtk_writer_4.SetFileName(f"vector_1_2_{index}.vtp")
    vtk_writer_4.SetInputData(vector_1_2_mesh)
    vtk_writer_4.Write()

    vtk_writer_5 = vtk.vtkXMLPolyDataWriter()
    vtk_writer_5.SetFileName(f"vector_2_3_{index}.vtp")
    vtk_writer_5.SetInputData(vector_2_3_mesh)
    vtk_writer_5.Write()

    vtk_writer_6 = vtk.vtkXMLPolyDataWriter()
    vtk_writer_6.SetFileName(f"vector_3_1_{index}.vtp")
    vtk_writer_6.SetInputData(vector_3_1_mesh)
    vtk_writer_6.Write()

    #We calculate the cosinus of the angles between the vectors
    cos_theta_ij = np.dot(vector_3_1, vector_2_3) / (np.linalg.norm(vector_3_1) * np.linalg.norm(vector_2_3))
    cos_theta_ik = np.dot(vector_1_2, vector_2_3) / (np.linalg.norm(vector_1_2) * np.linalg.norm(vector_2_3))

    #Calculate the estimated ratio of tension
    estimated_ratio = - (cos_theta_ij + cos_theta_ik)
    print(f"Estimated ratio of tension: {estimated_ratio}")


if __name__ == "__main__":


    #Path to the directory where the meshes are stored
    #path_to_mesh_folder = path.join(getcwd(), "simulation_results")


    calc_angle_between_cells(r"/home/steve/SimuCell3D_v2/scripts/python_bindings/cell_triplet_test_geom.vtk", 0)


    #Path to the mesh containing the cell doublet geometry
    #path_to_mesh = path.join(getcwd(), "test_geometry.vtk")



