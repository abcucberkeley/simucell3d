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
from scipy.spatial import distance_matrix

def round_to_nearest(value, base):
    value = np.asarray(value)

    value /= base
    value = np.round(value)
    value *= base
    return value

@profile
def uniform_mesh(n,
         input_mesh_path: Path = None,
         extent: float = 51e-6,  # boundary distance from zero in meters
         ):
    input = pv.read(input_mesh_path)
    input.translate(np.array(input.center) * -1, inplace=True)
    print(f'{input.array_names}')

    bounds = input.bounds
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

    mesh = input
    total = pv.MultiBlock()
    for i, position in tqdm(enumerate(grid), unit=' cells', desc='Creating meshes at each nucleation point', total=grid.shape[0]):
        if total is None:
            total = mesh.copy(deep=True)
        else:
            next = mesh.translate(position, inplace=False, transform_all_input_vectors=False)
            next['Cell_id'] = next['Cell_id'] + i
            total.append(next)
    total = total.combine()
    print(f"Final grid:\n {total}")
    print(f"Cell_id's: {total['Cell_id'][:20]}...")


    output_file = Path(f"{input_mesh_path.with_suffix('')}_meshed.vtk")
    print(f"Saving...")
    total.save(output_file, binary=False)
    print(f"Saved: {output_file.resolve()}")

if __name__ == '__main__':
    uniform_mesh(100, Path("../data/input_meshes/sphere.vtk"))