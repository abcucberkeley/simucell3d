import numpy as np
import pyvista as pv

from tqdm import tqdm
from pathlib import Path
from line_profiler_pycharm import profile
from scipy.spatial import distance_matrix
from vtkmodules.vtkIOLegacy import vtkUnstructuredGridWriter

from timeit import default_timer as timer
from datetime import timedelta

def round_to_nearest(value, base):
    value = np.asarray(value)

    value /= base
    value = np.round(value)
    value *= base
    return value

@profile
def uniform_mesh(
        input_mesh_path: Path = None,
        extent: float = 51e-6,  # boundary distance from zero in meters
        ):
    mesh = pv.read(input_mesh_path)
    mesh.translate(np.array(mesh.center) * -1, inplace=True)
    print(f'{mesh.array_names}')

    bounds = mesh.bounds
    size = (bounds[0]-bounds[1],
            bounds[2]-bounds[3],
            bounds[4]-bounds[5])
    size = max(np.abs(np.array(size)))

    estimated_n = round_to_nearest(extent**3 / size**3, 1).astype(int)

    print(f'extent = {extent}, input_mesh = {input_mesh_path.name}, input_mesh size = {size},\nroughly {estimated_n} points')
    rng = np.random.default_rng()

    n = estimated_n
    low = -extent/2 - size/2
    high = extent + size/2
    grid = np.array(rng.uniform(low=low, high=high, size=3))
    grid = np.expand_dims(grid, axis=0)
    for i in range(n):
        nucleation_point = np.array(rng.uniform(low=low, high=high, size=3))
        distances = distance_matrix(np.expand_dims(nucleation_point, axis=0), grid)
        if np.all(distances > size):
            grid = np.vstack((grid, nucleation_point))
    print(f"  found {grid.shape[0]} points that don't intersect")

    total = pv.MultiBlock()
    for i, position in tqdm(enumerate(grid), unit=' cells', desc='Creating meshes at each nucleation point', total=grid.shape[0]):
        next_mesh = mesh.translate(position, inplace=False, transform_all_input_vectors=True)
        next_mesh['Cell_id'] = next_mesh['Cell_id'] + i
        total.append(next_mesh)
    total = total.combine()
    print(f"Final grid:\n {total}, {total.n_cells=}, {total.number_of_cells=}")

    output_file = Path(f"{input_mesh_path.with_suffix('')}_meshed.vtk")
    print(f"Saving...")
    writer = vtkUnstructuredGridWriter()
    writer.SetFileVersion(42)
    writer.SetFileName(output_file)
    writer.SetInputData(total)
    writer.SetFileTypeToASCII()  # needs this for Simucell3d to load it
    writer.Write()
    # input.save(output_file, binary=False)  # Garbage.  Won't get read by
    print(f"Saved: {output_file.resolve()}")

    ## Check if file was written properly
    # total = pv.read(output_file)
    # print(f"Saved grid:\n {total}, {total.n_cells=}, {total.number_of_cells=}")


if __name__ == '__main__':
    start = timer()

    uniform_mesh(Path("../data/input_meshes/sphere.vtk"))

    end = timer()
    print(timedelta(seconds=end - start))
