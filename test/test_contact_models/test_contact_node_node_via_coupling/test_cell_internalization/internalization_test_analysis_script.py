import numpy as np 
import vtk 
from os import path, _exit, getcwd, listdir 
from scipy import optimize




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
def get_cell_lateral_face_lst(cell_id:int, mesh:vtk.vtkUnstructuredGrid()):
    """
    
    Get all the lateral faces belonging to a certain cell
    
    Parameters:
    -----------
    cell_id : int
        Id of the cell.

    mesh : vtk.vtkUnstructuredGrid
        Mesh containing the cells.

    Returns:
    --------
    lateral_face_lst : list
        List containing the lateral faces belonging to the cell.

        
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
    threshold_filter_2.ThresholdBetween(1, 1)
    threshold_filter_2.Update()

    return threshold_filter_2.GetOutput()
#------------------------------------------------------------------------------------------------------------------


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
#------------------------------------------------------------------------------------------------------------------





#------------------------------------------------------------------------------------------------------------------
def slice_faces(
        division_plane_center:np.array, 
        division_plane_normal: np.array, 
        mesh :  vtk.vtkUnstructuredGrid()
    ):
    """Slice the faces along the y axis and return the slice.
    
    Parameters:
    -----------
    division_plane_center : np.ndarray
        Center of the plane that will cut the cells.

    division_plane_normal : np.ndarray
        Normal of the plane that will cut the cells.

    mesh : vtk.vtkUnstructuredGrid
        Mesh containing the faces.

    Returns:
    --------
    slice_mesh : vtk.vtkPolyData
        Slice of the faces.
    """

    #Create the vtk plane that will cut the cells
    vtk_plane = vtk.vtkPlane()
    vtk_plane.SetOrigin(division_plane_center)
    vtk_plane.SetNormal(division_plane_normal)

    slice_filter = vtk.vtkPlaneCutter()
    slice_filter.SetInputData(mesh)
    slice_filter.SetPlane(vtk_plane)
    slice_filter.Update()

    return slice_filter.GetOutput()
#------------------------------------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------------------------------------
def fit_circle(mesh_center: np.ndarray, mesh: vtk.vtkUnstructuredGrid()):
    """
    Fit a circle to the points given as parameter.

    Parameters:
    -----------

    mesh_center : np.ndarray shape (3,)
        Center of the mesh.

    mesh : vtk.vtkUnstructuredGrid
        Mesh containing the points.

    Returns:
    --------

    circle_centroid : np.ndarray 3d
        Centroid of the circle.

    circle_radius : float
        Radius of the circle.
    """


    #Extract the points from the mesh
    points_cell_ar = np.array([mesh.GetPoint(i) for i in range(mesh.GetNumberOfPoints())])

    #Remove the y coordinates from the points
    points_cell_ar = np.delete(points_cell_ar, 1, 1)

    #Compute a first guess for the centroid of the circle
    centroid_1st_guess_ar = np.mean(points_cell_ar, axis=0)  

    #Fit a circle to the points
    fitted_centroid, optimization_output = optimize.leastsq(loss_function, centroid_1st_guess_ar, args=(points_cell_ar), Dfun=jacobian_function, col_deriv=True)
  
    #Compute the radius of the fitted circle
    radius = calc_distance_between_centroid_and_points(fitted_centroid, points_cell_ar).mean()

    #Add back the y dimension to the centroid 
    fitted_centroid = np.insert(fitted_centroid, 1, mesh_center[1])

    return fitted_centroid, radius
#------------------------------------------------------------------------------------------------------------------



#------------------------------------------------------------------------------------------------------------------
def save_circle(circle_center , circle_radius, output_path):
    """
    Save the circle in a vtk file. The circle is drawn in the xz plane.

    Parameters:
    -----------
    circle_center : np.ndarray shape (3,)
        Center of the circle.

    circle_radius : float
        Radius of the circle.

    output_path : str
        Path to the output file.
    
    Returns:
    --------

    None
    """

    #Create the circle
    circle = vtk.vtkRegularPolygonSource()
    circle.SetCenter(circle_center[0], circle_center[1], circle_center[2])
    circle.SetRadius(circle_radius)
    circle.SetNormal(0, 1, 0)
    circle.SetNumberOfSides(100)
    circle.Update()

    #Save the circle to a file
    vtk_writer = vtk.vtkXMLPolyDataWriter()
    vtk_writer.SetFileName(output_path)
    vtk_writer.SetInputData(circle.GetOutput())
    vtk_writer.Write()


#------------------------------------------------------------------------------------------------------------------

def find_circle_intersection(circle_center_1, circle_center_2, r1, r2 ):


    # Calculate the distance between the centers of the circles
    d = np.linalg.norm(circle_center_1 - circle_center_2)

    # Check if the circles are separate or one is contained within the other
    if d > r1 + r2 or d < abs(r1 - r2):
        # Circles are separate or one is contained within the other, no intersection
        return []

    # Calculate the intersection points
    a = (r1**2 - r2**2 + d**2) / (2 * d)
    h = np.sqrt(r1**2 - a**2)

    x1, y1, z1 = circle_center_1
    x2, y2, z2 = circle_center_2

    # Calculate the coordinates of the intersection points
    x3 = x1 + a * (x2 - x1) / d
    z3 = z1 + a * (z2 - z1) / d

    intersection1 = np.array([x3 + h * (z2 - z1) / d, y1, z3 - h * (x2 - x1) / d])
    intersection2 = np.array([x3 - h * (z2 - z1) / d, y1, z3 + h * (x2 - x1) / d])

    return [intersection1, intersection2]
#------------------------------------------------------------------------------------------------------------------



if __name__ == "__main__":

    #Get the path to the directory where this script is stored 
    path_to_script_folder = path.dirname(path.realpath(__file__))

    #Get the path to the mesh file
    path_to_mesh_file = path.join(path_to_script_folder, "test_geom_2_face_data.vtk")

    #Read the mesh
    ugrid = read_mesh(path_to_mesh_file)

    #Get the center of the mesh
    mesh_center = np.array(ugrid.GetCenter())

    #Get the apical faces of the  2 cells
    cell_0_apical_faces = get_cell_apical_face_lst(0, ugrid)
    cell_1_apical_faces = get_cell_apical_face_lst(1, ugrid)

    #Get teh lateral faces of the 2 cells
    cell_lateral_faces = get_cell_lateral_face_lst(0, ugrid)

    #Slice the faces along the y axis
    cell_0_apical_faces_slice = slice_faces(mesh_center, np.array([0.,1.,0.]), cell_0_apical_faces)
    cell_1_apical_faces_slice = slice_faces(mesh_center, np.array([0.,1.,0.]), cell_1_apical_faces)
    cell_lateral_faces_slice  = slice_faces(mesh_center, np.array([0.,1.,0.]), cell_lateral_faces)

    #Fit a circle to the 3 slices
    #Fit circles to the points
    cell_0_apical_faces_slice_centroid, cell_0_apical_faces_slice_radius = fit_circle(mesh_center, cell_0_apical_faces_slice)
    cell_1_apical_faces_slice_centroid, cell_1_apical_faces_slice_radius = fit_circle(mesh_center, cell_1_apical_faces_slice)
    cell_lateral_faces_slice_centroid,  cell_lateral_faces_slice_radius  = fit_circle(mesh_center, cell_lateral_faces_slice)

    #Save the 3 fitted circles
    save_circle(cell_0_apical_faces_slice_centroid, cell_0_apical_faces_slice_radius, path.join(path_to_script_folder, "cell_0_apical_faces_slice.vtp"))    
    save_circle(cell_1_apical_faces_slice_centroid, cell_1_apical_faces_slice_radius, path.join(path_to_script_folder, "cell_1_apical_faces_slice.vtp"))
    save_circle(cell_lateral_faces_slice_centroid,  cell_lateral_faces_slice_radius,  path.join(path_to_script_folder, "cell_lateral_faces_slice.vtp"))

    #Compute the intersection points between the circle cell_0_apical_faces_slice and the circle cell_1_apical_faces_slice
    intersection_points = find_circle_intersection(cell_0_apical_faces_slice_centroid, cell_1_apical_faces_slice_centroid, cell_0_apical_faces_slice_radius, cell_1_apical_faces_slice_radius)

    #Only the second intersection point is important
    intersection_point = intersection_points[1]

    center_1 = cell_0_apical_faces_slice_centroid
    center_2 = cell_1_apical_faces_slice_centroid
    radius_1 = cell_0_apical_faces_slice_radius
    radius_2 = cell_1_apical_faces_slice_radius

    d = np.linalg.norm(center_1  - center_2)
    if np.abs(radius_1 - radius_2) < d and d < radius_1 + radius_2:    
        l = np.sqrt((np.power(radius_1+radius_2, 2) - d*d)*(d*d - np.power(radius_1-radius_2, 2))) / d
        theta = np.arcsin(l / (2. * radius_1)) + np.arcsin(l / (2. * radius_2))
        print("theta", theta / 2.)


    #Print the angles between the center of the circles and the intersection points 
    R_c = cell_lateral_faces_slice_radius
    R_1 = cell_0_apical_faces_slice_radius
    R_2 = cell_1_apical_faces_slice_radius

    #The 3 raddii are located on the same line, project the intersection point on this line
    projected_intersection_point = np.array([intersection_point[0], cell_0_apical_faces_slice_centroid[1], cell_0_apical_faces_slice_centroid[2]])

    #Compute the angles
    psi_1 = np.arccos(np.dot(projected_intersection_point - cell_0_apical_faces_slice_centroid,  intersection_point - cell_0_apical_faces_slice_centroid) / (np.linalg.norm(projected_intersection_point - cell_0_apical_faces_slice_centroid) * np.linalg.norm(intersection_point - cell_0_apical_faces_slice_centroid)))
    psi_2 = np.arccos(np.dot(projected_intersection_point - cell_1_apical_faces_slice_centroid,  intersection_point - cell_1_apical_faces_slice_centroid) / (np.linalg.norm(projected_intersection_point - cell_1_apical_faces_slice_centroid) * np.linalg.norm(intersection_point - cell_1_apical_faces_slice_centroid)))
    psi_c = np.arccos(np.dot(projected_intersection_point - cell_lateral_faces_slice_centroid,   intersection_point - cell_lateral_faces_slice_centroid)  / (np.linalg.norm(projected_intersection_point - cell_lateral_faces_slice_centroid)  * np.linalg.norm(intersection_point - cell_lateral_faces_slice_centroid)))
    psi_c = np.pi - psi_c

    print("psi_1 = ", psi_1)
    print("psi_2 = ", psi_2)
    print("psi_c = ", psi_c)

    print("R1 = ", R_1)
    print("R2 = ", R_2)
    print("Rc = ", R_c)

