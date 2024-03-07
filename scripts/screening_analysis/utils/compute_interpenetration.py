import os
import sys
import numpy as np
import trimesh as tm
from sys import argv
from trimesh.collision import CollisionManager
from tqdm import tqdm
from itertools import product
from typing import Optional, Tuple, Union
from filelock import FileLock

sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))
from mesh_utils import get_cell_geometries_from_file


#----------------------------------------------------------------------------------------
def generate_points_on_sphere(
        center: np.ndarray[float], 
        radius: float, 
        num_points: int
    ) -> np.ndarray[float]:
    '''
    Given center coordinates and radius, generate a set of points on the surface of a sphere.

    Parameters:
    -----------
        center: (np.ndarray[float])
            An array of shape (3, ) defining the coordinates of the center of the sphere.

        radius: (float)
            The length of the radius of the sphere to generate.

        num_points: (int)
            The number of points to generate on the surface of the sphere.
            (Must be the square of an integer).

    Returns:
    --------
        coords: (np.ndarray[float]) 
            An array of shape (num_points, 3) storing the generated coordinates on the 
            sphere surface.
    '''
    # Generate grid of values on the ranges of theta and phi
    len_grid = int(np.sqrt(num_points))
    theta_grid = np.linspace(0, 2 * np.pi, len_grid)
    phi_grid = np.linspace(0, np.pi, len_grid)
    theta = [th for th, _ in product(theta_grid, phi_grid)]
    phi = [ph for _, ph in product(theta_grid, phi_grid)]

    # Convert spherical coordinates to Cartesian coordinates
    x = center[0] + radius * np.sin(phi) * np.cos(theta)
    y = center[1] + radius * np.sin(phi) * np.sin(theta)
    z = center[2] + radius * np.cos(phi)

    coords = np.column_stack((x, y, z))
    
    return coords
#----------------------------------------------------------------------------------------



#----------------------------------------------------------------------------------------
def are_interpenetrating(
        mesh_1: tm.Trimesh,
        mesh_2: tm.Trimesh,
        interp_radius: float
) -> bool:
    '''
    Check if two meshes are interpenetrating.

    Parameters:
    -----------

    mesh_1: (tm.Trimesh)
        A cell meshes in the `trimesh` format.
    
    mesh_2: (tm.Trimesh)
        Another cell meshes in the `trimesh` format.

    interp_radius: (float)
        The radius of the sphere used to check interpenetration.

    Returns:
    --------

    `True`, if the 2 meshes interpenetrate each other.
    '''

    # Instantiate a collision manager
    manager = CollisionManager()

    # Add meshes to the manager
    manager.add_object(name='1', mesh=mesh_1)
    manager.add_object(name='2', mesh=mesh_2)

    # Get colliding idxs
    collision, contacts = manager.in_collision_internal(return_data=True)

    # Get intersection points and their barycenter
    if collision:
        contact_pts = np.asarray([contact.point for contact in contacts])
        collision_barycenter = np.mean(contact_pts, axis=0)
        collision_sphere = generate_points_on_sphere(
            collision_barycenter, radius=interp_radius, num_points=25
        )
        return np.all(mesh_1.contains(collision_sphere)) and np.all(mesh_2.contains(collision_sphere))
    else:
        return False
#----------------------------------------------------------------------------------------



#----------------------------------------------------------------------------------------
def check_interpenetration_in_geometry(
        path_to_geometry_file: str,
        avg_edge_len: Optional[float] = None,
        remove_shell: Optional[bool] = True,
) -> Union[Tuple[int, int], float]:
    '''
    Given the path to a .vtk file containing a list of cell meshes, return a list
    of indexes of colliding cells.

    Parameters:
    -----------

    path_to_geometry_file: (str)
        The path to a .vtk file containing a list of cell meshes.
    
    remove_shell: (Optional[bool], default=True)
        If `True` remove the cell mesh corresponding to the shell
        (i.e., the largest mesh).

    avg_edge_len: (Optional[float], default=None)
        The average edge length over all the cell meshes in the geometry.
        If `None`, `avg_edge_len` is computed from scratch.
    
    Returns:
    --------

    (Union[Tuple[int, int], float])
        If interpenetrating cells are found, a tuple of indexes of the first 2 interpenetrating cells is returned.
        Otherwise the current value of avg_edge_len is returned.
    '''
    
    cell_mesh_lst = get_cell_geometries_from_file(path_to_geometry_file, None)

    # Remove ECM shell, since not interesting to check interpenetration
    if remove_shell:
        # remove largest mesh (shell)
        num_vertices = [len(mesh.vertices) for mesh in cell_mesh_lst]
        try:
            largest_mesh_idx, _ = max(enumerate(num_vertices), key=lambda x: x[1])
        except:
            print(f'Number of cells: {len(cell_mesh_lst)}')
            print(f'Number of vertices per mesh: {num_vertices}')
        del cell_mesh_lst[largest_mesh_idx]

    # If necessary, compute the average edge length across the cell meshes in the current geometry
    if not avg_edge_len:
        edge_lens = []
        for mesh in cell_mesh_lst:
            edges = mesh.edges_unique
            vertices = mesh.vertices
            edges_length = list(map(
                lambda x: np.linalg.norm(vertices[x[0]] - vertices[x[1]]),
                edges
            ))
            edge_lens.append(np.mean(edges_length))
        avg_edge_len = np.mean(edge_lens)

    # Check interpenetration for each pair of meshes
    for i in tqdm(range(len(cell_mesh_lst))):
        for j in range(i+1, len(cell_mesh_lst)):
            if are_interpenetrating(cell_mesh_lst[i], cell_mesh_lst[j], 2*avg_edge_len):
                return i, j

    return avg_edge_len
#----------------------------------------------------------------------------------------



#----------------------------------------------------------------------------------------
def check_simulation_run_interpenetration(
        path_to_mesh_dir: str,
        check_every: Optional[int] = 4,
        show_info: Optional[bool] = True
) -> bool:
    '''
    Given a path to the directory that stores .vtk files for the different
    simulation iterations (e.g., ./simulation_outputs/simulation_{id}/cell_data),
    returns `True` if at any iteration, two meshes interpenetrate each other. 

    NOTE:
    The function stops checking once reached an iteration showing interpenetration.
    
    Parameters:
    -----------

    path_to_mesh_dir: (str)
        The path to the directory that stores .vtk files for the different
        simulation iterations (e.g., ./simulation_outputs/simulation_{id}/cell_data)
    
    check_every: (Optional[int], default=3)
        Specify every how many iterations interpenetration is checked.
    
    show_info: (Optional[true], default=True)
        If `True`, print the info about when interpenetration is detected.

    Returns:
    --------

    (bool)
        `True` is at least two cells interpenetrate at any iteration. 
    '''
    # Note: it is more likely to  have interpenetration at the last iterations.
    # Therefore, we start our search from there
    num_files = len(os.listdir(path_to_mesh_dir))

    count = 1
    average_edge_length = None # needed to check interpenetration
    interp_idxs = None
    for i in range(num_files, 0, -check_every):
        mesh_file = f'result_{i}.vtk'
        print('-----------------------------------------------------------')
        print(f'Analyzing iteration {count}/{(int(num_files//check_every))} -> {mesh_file} ...')
        count += 1

        iter_id = mesh_file.split('.')[0].split('_')[1]

        res = check_interpenetration_in_geometry(
            path_to_geometry_file=os.path.join(path_to_mesh_dir, mesh_file),
            avg_edge_len=average_edge_length,
            remove_shell=True
        )    

        # If `res` is a tuple of int -> interp_idxs
        # If `res` is a float -> average_edge_length
        if isinstance(res, float):
            average_edge_length = res
        elif isinstance(res, tuple):
            interp_idxs = res 

        if interp_idxs:
            if show_info:
                print(f'Found interpenetration between cells: {interp_idxs} at iteration {iter_id}.')
            return True
        
    if show_info:
        print('No interpenetration between cells was found.')
    return False
#----------------------------------------------------------------------------------------



if __name__ == '__main__':

    #Check that the name of the screen was given in command line input
    assert argv.__len__() == 3 , "Wrong command line arguments:"+\
        "python3 check_mesh_collision.py path/to/simulation_outputs/sim_{id}/cell_data/ path/to/output/file.txt"

    #Get the arguments
    path_to_cell_meshes = argv[1].strip()
    assert os.path.exists(path_to_cell_meshes), "The given cell_data folder does not exist: " + path_to_cell_meshes
    output_file_path = argv[2].strip()
    assert os.path.isfile(output_file_path), "The given output file does not exist: " + output_file_path

    print(f'Analyzing {path_to_cell_meshes}')

    # Run the algorithm
    res = check_simulation_run_interpenetration(
        path_to_mesh_dir=path_to_cell_meshes,
        check_every=10, 
        show_info=True
    )

    # Acquire the file lock before writing to the output file
    sim_id = os.path.basename(path_to_cell_meshes.replace('/cell_data', '')).split('_')[1]
    with FileLock(output_file_path + ".lock"):
        with open(output_file_path, "a") as file:
            # Write the result to the output file
            file.write(f'{sim_id} {str(res)} \n')
