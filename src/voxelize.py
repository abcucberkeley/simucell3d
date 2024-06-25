"""
Voxelize a Simulation
~~~~~~~~~~~~~~~~~~~~~~~

Create a voxel model (like legos) of a closed surface or volumetric mesh.

1. First you need to build Docker Image:

docker build -t simucell3d_docker_img .


2. Next, you need to generate simulation results (these get saved in a mesh file format: .vtk)

docker run -it -v ${pwd}/simulation_results:/app/simulation_results -v ${pwd}/data:/data -v ${pwd}:/app/params simucell3d_docker_img:latest /app/params/parameters_default_dynamic_dan.xml


3.  Then you can voxelize the results


Output will be two files:
    1. outline (basically the bubbles' boundaries)
    2. labels (should be the filled interior of each bubble)

"""


import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt
import matplotlib.collections as collections

from tqdm import tqdm
from tifffile import imwrite
from pathlib import Path
import os
import shutil
from line_profiler_pycharm import profile

from timeit import default_timer as timer
from datetime import timedelta


def round_to_nearest(value, base):
    value = np.asarray(value)

    value /= base
    value = np.round(value)
    value *= base
    return value


@profile
def main(shape: tuple = (256, 1024, 1024),
         input_dir: Path = Path.cwd().joinpath(r"simulation_results/dynamic_simulation/face_data"),
         extent: float = 1e-4,  # Z boundary distance from zero in meters
         ):

    volume_outlines = np.zeros(shape=shape, dtype=np.uint16)
    volume_labels = np.zeros_like(volume_outlines)

    output_dir_outlines = Path(os.path.join(input_dir.parent, "voxelized_outlines"))
    output_dir_labels = Path(os.path.join(input_dir.parent, "voxelized_labels"))

    shutil.rmtree(output_dir_outlines, ignore_errors=True)
    shutil.rmtree(output_dir_labels, ignore_errors=True)
    output_dir_outlines.mkdir(parents=True, exist_ok=True)
    output_dir_labels.mkdir(parents=True, exist_ok=True)

    print(f"Output files: {output_dir_outlines}  and  {output_dir_labels}")

    meshfiles = sorted(Path(input_dir).iterdir(), key=os.path.getmtime, reverse=True)
    # meshfiles = meshfiles[1401::100]

    bounds = None

    for meshfile in tqdm(meshfiles, unit=" files", position=0, leave=True):
        if meshfile.suffix.endswith(".vtk"):

            print(f"Reading {meshfile}")
            mesh = pv.read(meshfile)

            # meshsize = np.array([mesh.bounds[0]-mesh.bounds[1],
            #                      mesh.bounds[2]-mesh.bounds[3],
            #                      mesh.bounds[4]-mesh.bounds[5]])
            #
            # bounds_okay = all(np.abs(meshsize) > 1e-4)
            # print(f"Bounds ok: {bounds_okay}  {meshfile.name}")
            # if bounds_okay:
            #     continue
            # else:
            #     break

            if bounds is None:
                bounds = mesh.bounds

            # print(mesh.array_names)
            mesh.set_active_scalars('face_cell_id')
            surf = mesh.extract_surface()


            CameraPosition = (-0.00028617924571107807, 1.2141999718551233e-05, 1.2238500858074985e-05)
            CameraFocalPoint = (1.212385041071684e-05, 1.2141999718551233e-05, 1.2238500858074985e-05)
            CameraViewUp = (0.0, 0.0, 1.0)
            cpos = [CameraPosition, CameraFocalPoint, CameraViewUp]
            # mesh.plot(scalars='face_cell_id', cpos=cpos, opacity=0.99, specular=0.3)

            center = np.array([0,0,0])
            slice_normal = np.array([0, 0, 1])   # stupid (x,y,z)
            point_a = center + slice_normal * extent
            point_b = center - slice_normal * extent
            z_slice_locations = pv.Line(point_a, point_b, resolution=volume_outlines.shape[0]-1)

            # p = pvqt.BackgroundPlotter()
            # pv.global_theme.allow_empty_mesh = True
            # p.show_bounds(grid=True, location='back')
            # p.title = f'{meshfile.stem}'

            name = 'face_cell_id'
            plt.style.use('dark_background')

            dpi = 1000
            figsize = (volume_outlines.shape[2] / dpi, volume_outlines.shape[1] / dpi)
            fig_outlines, ax_outlines = plt.subplots(num=10, clear=True, figsize=figsize, dpi=dpi)
            fig_labels, ax_labels = plt.subplots(num=11, clear=True, figsize=figsize, dpi=dpi)

            for z_slice_idx, line_point in enumerate(z_slice_locations.points):
                slice = surf.slice(normal=slice_normal, origin=line_point)
                unique_cell_ids = np.unique(slice[name]).astype(int)

                ax_labels.cla()
                ax_outlines.cla()


                # for c in ax_outlines.collections:  # possibly better to use: for c in plt.lines (see comment)
                #     if c.get_gid() == id:
                #         c.remove()
                #
                # for c in ax_labels.collections:  # possibly better to use: for c in plt.lines (see comment)
                #     if c.get_gid() == id:
                #         c.remove()

                for color_idx, unique_cell_id in enumerate(unique_cell_ids):
                    single_outline = (slice
                                      .threshold([unique_cell_id, unique_cell_id], scalars=name)
                                      .extract_surface())

                    starts = single_outline.lines[1:][::3]
                    ends = single_outline.lines[2:][::3]
                    x_starts = single_outline.points[starts][:, 0]
                    y_starts = single_outline.points[starts][:, 1]
                    x_ends = single_outline.points[ends][:, 0]
                    y_ends = single_outline.points[ends][:, 1]

                    x = np.vstack([x_starts, x_ends])
                    y = np.vstack([y_starts, y_ends])

                    segments = np.array([x, y]).transpose()

                    Blue = unique_cell_id & 255
                    Green = (unique_cell_id >> 8) & 255
                    Red = (unique_cell_id >> 16) & 255
                    color = (Red / 255, Green / 255, Blue / 255)

                    # Create a LineCollection object
                    lc = collections.LineCollection(segments, colors=color, antialiased=False, linewidths=.1)
                    ax_outlines.add_collection(lc)

                    destinations = []
                    for a, b in zip(x_ends, y_ends):
                        idxs = np.where(x_starts == a)
                        if len(idxs[0]) == 0:
                            index = np.argmin(np.abs(x_starts-a))
                            destinations.append(index)
                        else:
                            done = False
                            for idx in idxs[0]:
                                if y_starts[idx] == b and not done:
                                    destinations.append(idx)
                                    done = True
                            if not done:
                                index = np.argmin(np.abs(y_starts - b))
                                destinations.append(index)

                    ordered_destinations = []
                    i = 0
                    for destination in destinations:
                        ordered_destinations.append(i)
                        i = destinations[i]

                    x_patch_verts = x_starts[ordered_destinations]
                    y_patch_verts = y_starts[ordered_destinations]

                    patches = np.array([x_patch_verts, y_patch_verts]).transpose()


                    Blue = unique_cell_id & 255
                    Green = (unique_cell_id >> 8) & 255
                    Red = (unique_cell_id >> 16) & 255
                    color = (Red / 255, Green / 255, Blue / 255)

                    # Create a PolyCollection object (a list of patches)
                    pc = collections.PolyCollection(verts=[patches], facecolors=color, antialiased=False, closed=True)
                    ax_labels.add_collection(pc)

                ax_outlines.set_xlim([bounds[0], bounds[1]])
                ax_outlines.set_ylim([bounds[2], bounds[3]])
                ax_outlines.axis('off')

                fig_outlines.canvas.draw()

                # Convert the canvas to a raw RGB buffer, next go back to the u16 we wanted.
                buf = fig_outlines.canvas.renderer.buffer_rgba()
                ncols, nrows = fig_outlines.canvas.get_width_height()
                image = np.frombuffer(buf, dtype=np.uint8).reshape(nrows, ncols, 4).astype(int)
                r = image[:, :, 0]
                g = image[:, :, 1]
                b = image[:, :, 2]

                slab = r * 256**2 + g * 256 + b     # convert RGB back to face_cell_id
                volume_outlines[z_slice_idx, :, :] = slab.astype(np.uint16)    # assign to the slice in the volume



                ax_labels.set_xlim([bounds[0], bounds[1]])
                ax_labels.set_ylim([bounds[2], bounds[3]])
                ax_labels.axis('off')

                fig_labels.canvas.draw()

                # Convert the canvas to a raw RGB buffer, next go back to the u16 we wanted.
                buf = fig_labels.canvas.renderer.buffer_rgba()
                ncols, nrows = fig_labels.canvas.get_width_height()
                image = np.frombuffer(buf, dtype=np.uint8).reshape(nrows, ncols, 4).astype(int)
                r = image[:, :, 0]
                g = image[:, :, 1]
                b = image[:, :, 2]

                slab = r * 256**2 + g * 256 + b     # convert RGB back to face_cell_id
                volume_labels[z_slice_idx, :, :] = slab.astype(np.uint16)    # assign to the slice in the volume

                # plt.close('all')

        output_file_outlines = Path(f'{output_dir_outlines}/{meshfile.stem}_outlines.tif')
        imwrite(output_file_outlines, volume_outlines.astype(np.uint16), compression='deflate', dtype=np.uint16)
        output_file_labels = Path(f'{output_dir_labels}/{meshfile.stem}_labels.tif')
        imwrite(output_file_labels, volume_labels.astype(np.uint16), compression='deflate', dtype=np.uint16)

if __name__ == '__main__':
    start = timer()

    main()

    end = timer()
    print(timedelta(seconds=end - start))
