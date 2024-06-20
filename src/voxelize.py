"""
.. _voxelize_surface_mesh_example:

Voxelize a Surface Mesh
~~~~~~~~~~~~~~~~~~~~~~~

Create a voxel model (like legos) of a closed surface or volumetric mesh.

"""


#    docker run -it -v ${pwd}/simulation_results:/app/simulation_results  simucell3d_docker_img:latest ../parameters_default_dynamic_dan.xml



import numpy as np
from scipy.interpolate import LinearNDInterpolator
import pyvista as pv
import pyvistaqt as pvqt
import fast_simplification
import matplotlib.pyplot as plt
import matplotlib.collections as collections

from tqdm import tqdm
from tifffile import imwrite
from pathlib import Path
import os
import shutil
from line_profiler_pycharm import profile

def round_to_nearest(value, base):
    value = np.asarray(value)

    value /= base
    value = np.round(value)
    value *= base
    return value


def slice_to_array(slc, normal, origin, name, ni=50, nj=50):
    """Converts a PolyData slice to a 2D NumPy array.

    It is crucial to have the true normal and origin of
    the slicing plane

    Parameters
    ----------
    slc : PolyData
        The slice to convert.
    normal : tuple(float)
        the normal of the original slice

        ** only works if normal is along z! Bug in code where z = np.array([0]).  Need to flatten to 2D after rotation I think, not before.

    origin : tuple(float)
        the origin of the original slice
    name : str
        The scalar array to fetch from the slice
    ni : int
        The resolution of the array in the i-direction
    nj : int
        The resolution of the array in the j-direction

    """
    # Make structured grid
    x = np.linspace(slc.bounds[0], slc.bounds[1], ni)
    y = np.linspace(slc.bounds[2], slc.bounds[3], nj)
    z = np.array([0])
    plane = pv.StructuredGrid(*np.meshgrid(x, y, z))

    vx = np.array([0., 0., 1.])
    direction = normal / np.linalg.norm(normal)
    if np.array_equal(vx, direction):
        pass
    else:
        # rotate and translate grid to be ontop of the slice
        vx -= vx.dot(direction) * direction
        vx /= np.linalg.norm(vx)
        vy = np.cross(direction, vx)
        rmtx = np.array([vx, vy, direction])
        plane.points = plane.points.dot(rmtx)

    plane.points -= plane.center
    plane.points += origin

    # resample the data
    sampled = plane.sample(slc, tolerance=slc.length * 0.5)
    # Fill bad data
    sampled[name][~sampled["vtkValidPointMask"].view(bool)] = np.nan

    # plot the 2D array
    array = sampled[name].reshape(sampled.dimensions[0:2])
    return array

@profile
def voxelize_volume_with_bounds(mesh, density, check_surface=True):
    """Voxelize mesh to create a RectilinearGrid voxel volume.

    Creates a voxel volume that encloses the input mesh and discretizes the cells
    within the volume that intersect or are contained within the input mesh.
    ``InsideMesh``, an array in ``cell_data``, is ``1`` for cells inside and ``0`` outside.

    Parameters
    ----------
    mesh : pyvista.DataSet
        Mesh to voxelize.

    bounds :

    check_surface : bool, default: True
        Specify whether to check the surface for closure. If on, then the
        algorithm first checks to see if the surface is closed and
        manifold. If the surface is not closed and manifold, a runtime
        error is raised.

    Returns
    -------
    pyvista.RectilinearGrid
        RectilinearGrid as voxelized volume with discretized cells.

    See Also
    --------
    pyvista.voxelize
    pyvista.DataSetFilters.select_enclosed_points

    Examples
    --------
    """

    # check and pre-process input mesh
    # unique_cell_ids = np.unique(mesh['face_cell_id'])


    # need to use a plane to intersect the mesh, then extact those points and voxelize that.

    surface = mesh.extract_geometry()
    # surface = mesh.extract_points(mesh['face_cell_id'] == unique_cell_ids[0]).extract_geometry()  # filter preserves topology
    if not surface.faces.size:
        # we have a point cloud or an empty mesh
        raise ValueError('Input mesh must have faces for voxelization.')
    if not surface.is_all_triangles:
        # reduce chance for artifacts, see gh-1743
        surface.triangulate(inplace=True)


    if density is None:
        density = mesh.length / 100
    if isinstance(density, (int, float, np.number)):
        density_x, density_y, density_z = [density] * 3
    elif isinstance(density, np.ndarray):
        density_x, density_y, density_z = density
    else:
        raise TypeError(f'Invalid density {density!r}, expected number or array-like.')


    x_min, x_max, y_min, y_max, z_min, z_max = mesh.bounds
    x = np.arange(x_min, x_max, density_x)
    y = np.arange(y_min, y_max, density_y)
    z = np.arange(z_min, z_max, density_z)

    # Create a RectilinearGrid
    voi = pv.RectilinearGrid(x, y, z)

    # get part of the voi within the mesh's bounding surface.
    selection = voi.select_enclosed_points(surface, tolerance=0.0, check_surface=False)  # takes the longest
    mask_vol = selection.point_data['SelectedPoints'].view(np.bool_)

    # Get voxels that fall within input mesh boundaries
    # 1. Take the voi points that are enclosed by the surface (extract_points).
    # 2. Get what Cell IDs these associate with  ["vtkOrginialCellIds"]
    # 3. Just take the unique cell ids (if a cell encloses two or more surface points, just return its ID once)
    # 4. Set a cell value ("InsideMesh") for the voi according to these

    voi_points_with_mesh = voi.extract_points(np.argwhere(mask_vol))

    cell_ids, idx = np.unique(voi_points_with_mesh["vtkOriginalCellIds"], return_index=True)
    # face_cell_ids = voi_points_with_mesh[idx]['face_cell_id']

    # Create new element of grid where all cells _within_ mesh boundary are
    # given new name 'MeshCells' and a discrete value of 1
    voi['InsideMesh'] = np.zeros(voi.n_cells)
    voi['InsideMesh'][cell_ids] = 1

    return voi



input_dir =  Path.cwd().joinpath(r"simulation_results/dynamic_simulation/face_data")

output_dir = Path(os.path.join(input_dir.parent, "voxelized"))
shutil.rmtree(output_dir, ignore_errors=True)
output_dir.mkdir(parents=True, exist_ok=True)


meshfiles = sorted(Path(input_dir).iterdir(), key=os.path.getmtime, reverse=True)

for meshfile in tqdm(meshfiles, unit=" files", position=0, leave=True):
    if meshfile.suffix.endswith(".vtk"):

        print(f"Reading {meshfile}")
        mesh = pv.read(meshfile)

        # print(mesh.array_names)
        mesh.set_active_scalars('face_cell_id')
        surf = mesh.extract_surface()


        CameraPosition = (-0.00028617924571107807, 1.2141999718551233e-05, 1.2238500858074985e-05)
        CameraFocalPoint = (1.2123850410716841e-05, 1.2141999718551233e-05, 1.2238500858074985e-05)
        CameraViewUp = (0.0, 0.0, 1.0)
        cpos = [CameraPosition, CameraFocalPoint, CameraViewUp]
        # mesh.plot(scalars='face_cell_id', cpos=cpos, opacity=0.99, specular=0.3)

        volume = np.zeros((64,1024,1024), dtype=np.uint16)


        # mesh_coords = mesh.cell_centers().points
        # mesh_values = mesh.get_array('face_cell_id')
        # vol_x = np.linspace(0, 1, volume.shape[2])
        # vol_y = np.linspace(0, 1, volume.shape[1])
        # vol_z = np.linspace(0, 1, volume.shape[0])
        #
        # v_xi, v_yi, v_zi = np.meshgrid(vol_x, vol_y, vol_z, indexing='xy')
        #
        # # interpolate
        # interp = LinearNDInterpolator(mesh_coords, mesh_values)
        #
        # volume = interp(vol_x, vol_y, vol_z)
        # np.griddata(mesh_coords, mesh_values, (v_xi, v_yi, v_zi), method='linear')

        center = np.array([0,0,0])
        slice_normal = np.array([0, 0, 1])   # stupid (x,y,z)
        point_a = center + slice_normal * 1e-4
        point_b = center - slice_normal * 1e-4
        myline = pv.Line(point_a, point_b, resolution=volume.shape[0]-1)

        # p = pvqt.BackgroundPlotter()
        # pv.global_theme.allow_empty_mesh = True
        # p.show_bounds(grid=True, location='back')
        # p.title = f'{meshfile.stem}'

        name = 'face_cell_id'
        cm = plt.colormaps['tab20c']
        plt.style.use('dark_background')

        if True:
            for idx, line_point in enumerate(myline.points):

                slice = surf.slice(normal=slice_normal, origin=line_point)
                # p.add_mesh(slice)

                unique_cell_ids = np.unique(slice[name]).astype(int)
                fig, ax = plt.subplots(figsize=(1.024, 1.024), dpi=1000)
                for color_idx, unique_cell_id in enumerate(unique_cell_ids):
                    single = slice.threshold([unique_cell_id, unique_cell_id], scalars=name).extract_surface()


                    starts = single.lines[1:][::3]
                    ends = single.lines[2:][::3]
                    x_starts = single.points[starts][:, 0]
                    y_starts = single.points[starts][:, 1]
                    x_ends = single.points[ends][:, 0]
                    y_ends = single.points[ends][:, 1]

                    x = np.vstack([x_starts, x_ends])
                    y = np.vstack([y_starts, y_ends])

                    segments = np.array([x, y]).transpose()


                    # Create a LineCollection object
                    Blue = unique_cell_id & 255
                    Green = (unique_cell_id >> 8) & 255
                    Red = (unique_cell_id >> 16) & 255

                    color = (Red / 255, Green / 255, Blue / 255)
                    lc = collections.LineCollection(segments, colors=color, antialiased=False, linewidths=.1, facecolors=color)
                    ax.add_collection(lc)

                ax.set_xlim([mesh.bounds[0], mesh.bounds[1]])
                ax.set_ylim([mesh.bounds[2], mesh.bounds[3]])
                ax.axis('off')
                # plt.show()

                fig.canvas.draw()

                # Convert the canvas to a raw RGB buffer then back to the u16 we
                buf = fig.canvas.renderer.buffer_rgba()
                ncols, nrows = fig.canvas.get_width_height()
                image = np.frombuffer(buf, dtype=np.uint8).reshape(nrows, ncols, 4).astype(int)
                r = image[:, :, 0]
                g = image[:, :, 1]
                b = image[:, :, 2]

                slab = r * 256**2 + g * 256 + b     # convert RGB back to face_cell_id

                volume[idx, :, :] = slab    # assign to the slice in the volume
                output_file = Path(f'{output_dir}/{meshfile.stem}.tif')
                imwrite(output_file, volume.astype(np.uint16))
                plt.close()
