"""
.. _voxelize_surface_mesh_example:

Voxelize a Surface Mesh
~~~~~~~~~~~~~~~~~~~~~~~

Create a voxel model (like legos) of a closed surface or volumetric mesh.

"""


#    docker run -it -v ${pwd}/simulation_results:/app/simulation_results  simucell3d_docker_img:latest ../parameters_default_dynamic_dan.xml



import numpy as np
import pyvista as pv
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

    # get part of the mesh within the mesh's bounding surface.
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

        mesh = pv.read(meshfile)
        # print(mesh.array_names)
        mesh.set_active_scalars('face_cell_id')
        unique_cell_ids = np.unique(mesh['face_cell_id'])

        ###############################################################################
        CameraPosition = (-0.00028617924571107807, 1.2141999718551233e-05, 1.2238500858074985e-05)
        CameraFocalPoint = (1.2123850410716841e-05, 1.2141999718551233e-05, 1.2238500858074985e-05)
        CameraViewUp = (0.0, 0.0, 1.0)

        cpos = [CameraPosition, CameraFocalPoint, CameraViewUp]
        # mesh.plot(scalars='face_cell_id', cpos=cpos, opacity=0.99, specular=0.3)


        ###############################################################################
        # Create a voxel model of the bounding surface

        density = 10e-6

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

        # volume = np.zeros((len(x),len(y),len(z)), dtype=int)
        volume = np.zeros((256,256,256), dtype=int)

        if False:

            for i in tqdm(range(mesh.n_cells), desc=f"{meshfile.stem}", unit=' mesh_cells', position=0):
                points = mesh.get_cell(i).points
                # points = mesh.extract_all_edges().cell_centers().points
                for point in points:
                    cell_location = point
                    # cell_location = round_to_nearest(point, 1e-7)
                    # cell_location = round_to_nearest(mesh.get_cell(i).center, 1e-7)
                    cell_location = np.flip(cell_location, axis=0)  # z y x

                    volume_index = np.round((cell_location - np.array([z_min, y_min, x_min])) / density).astype(int)
                    volume_index = np.clip(volume_index, 0, np.array(volume.shape)-1)

                    face_id = mesh['face_cell_id'][i].astype(int)
                    volume[tuple(volume_index)] = face_id


                # print(f"{volume_index}   id={face_id}")

            output_file = Path(f'{output_dir}/{meshfile.stem}.tif')
            imwrite(output_file, volume)
            # print(f'Done. {volume.shape=}  shape={len(x), len(y), len(z)}, unique_cells={unique_cell_ids}, output file={output_file}')



        total_mesh = mesh.threshold([min(unique_cell_ids), max(unique_cell_ids)], scalars='face_cell_id')
        voxels = voxelize_volume_with_bounds(total_mesh, density=density)
        dimensions = np.array(voxels.dimensions) - 1
        dimensions = np.flip(dimensions, axis=0)
        imwrite(f'total.tif', voxels['InsideMesh'].reshape(dimensions))
        print(f'Done. {dimensions=}')

        # for unique_cell_id in unique_cell_ids:
        #     single = mesh.threshold([unique_cell_id, unique_cell_id], scalars='face_cell_id')
        #     single_biocell = pv.voxelize_volume(single, density=density)
        #     dimensions = np.array(single_biocell.dimensions) - 1
        #     dimensions = np.flip(dimensions, axis=0)
        #     single_cell = single_biocell['InsideMesh'].reshape(dimensions)
        #     imwrite(f'junk{unique_cell_id}.tif', single_cell)
        #     print(f'{unique_cell_id=} {dimensions=}')



        # print('Visualization.  Press "q" to quit and continue.')
        # p = pv.Plotter()
        # p.add_mesh(voxels, color=True, show_edges=False, opacity=0.5)
        # p.add_mesh(mesh, color="lightblue", opacity=0.5)
        # p.show(cpos=cpos)
        break
