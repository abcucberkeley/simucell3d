import pandas as pd
import numpy as np
from itertools import product
from typing import Dict

"""
Generate some fake dataframes to test the dashboard plotting function.
"""

np.random.seed(1234)


# Set the screened parameters
surface_tension = np.geomspace(3e-4, 2.5e-3, 10)
repulsion_strength = np.array([1e8, 5e8, 1e9])
time_step = np.array([5e-9, 1e-8, 5e-8, 1e-7])
parameters = np.asarray(
    [
        [ts, rep, ten] 
        for ts, rep, ten in product(time_step, repulsion_strength, surface_tension) 
])
num_records = parameters.shape[0]


# Simulate the other features
sim_id = np.arange(1, num_records + 1)

interpenetration = np.random.choice([0, 1], size=num_records, p=[0.95, 0.05]).astype(bool)

sphericity_mean = np.random.normal(10, 20, num_records)
sphericity_std = np.random.normal(20, 5, num_records)
sphericity = [[mean, std] for mean, std in zip(sphericity_mean, sphericity_std)]

vol_loss_mean = np.random.normal(2, 5, num_records)
vol_loss_std = np.random.normal(5, 1, num_records)
vol_loss = [[mean, std] for mean, std in zip(vol_loss_mean, vol_loss_std)]

num_iters = np.asarray([
    np.random.randint(0, n, 1)[0]
    for n in range(10, 10 + 10 * num_records, 10)
])
has_crashed = np.random.choice([0, 1], size=num_records, p=[0.95, 0.05]).astype(bool)
num_iterations = [[nit, cr] for nit, cr in zip(num_iters, has_crashed)]

IoU = np.random.uniform(0, 0.45, num_records)
IoU_deriv = np.random.normal(0.025, 0.05, num_records)


# Create the dataframe
column_names = [
    'sim_id', 'surface_tension', 'sphericity', 'interpenetration',
    'volume_loss', 'num_iterations', 'IoU', 'IoU_derivative'
]
data_dict = {
    'sim_id': sim_id,
    'time_step': parameters[:, 0],
    'repulsion_strength': parameters[:, 1],
    'surface_tension': parameters[:, 2], 
    'delta_sphericity': sphericity,
    'interpenetration': interpenetration,
    'volume_loss': vol_loss,
    'num_iterations': num_iterations,
    'IoU': IoU,
    'IoU_derivative': IoU_deriv
}

test_df = pd.DataFrame(data=data_dict, index=None)
test_df.to_csv("./misc/synthetic_data_3_parameters.csv")