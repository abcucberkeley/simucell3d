# Parameter screening

Several scripts have already been written to launch parameter screenings on clusters managed by Slurmp like the Euler cluster.  These scripts are located in the `scripts/parameter_screening` folder. You will have to adapt the last section of the script `launch_parameter_screening.py` if your cluster has a different job management system than slurmp.

<br> 
<br> 

## Installation

Load the required modules
```
env2lmod
module load gcc/9.3.0 openblas/0.3.20 python/3.11.2 cmake/3.25.0
```

Create a virtual environment
```
cd path/to/SimuCell3D/scripts/parameter_screening
python -m venv --system-site-packages venv_screening
```

Activate the virtual env
```
source ./venv_screening/bin/activate
```

Install pandas 
```
OPENBLAS=$OPENBLAS_ROOT/lib/libopenblas.so pip install --ignore-installed --no-deps pandas==1.4.4
```

</br>
</br>

## Launching a parameter screening
The parameter set of each simulation should be stored in the rows of a csv file like `scripts/parameter_screening/parameter_table.csv`: 

sim_id | input_mesh | time_step | damping | epi_bulk_modulus | ... 
---    | ---        | ---        | ---    |---               | ---
1 | .\data\cube | 1e-7 | 1e6 | 2500 | ... |  
2 | .\data\cube | 1e-7 | 1e6 |5000 | ... |  
3 | .\data\cube | 1e-7 | 1e6 |7500 | ... |  

Each row corresponds to a simulation. 

</br>

Start by loading the required modules
```
env2lmod
module load gcc/9.3.0 openblas/0.3.20 python/3.11.2 cmake/3.25.0
```

Activate the virtual env
```
cd path/to/SimuCell3D/scripts/parameter_screening
source ./venv_screening/bin/activate
```



Call the script that automatically launches the screening:
```
python3 launch_parameter_screening.py [path_to_parameter_table.csv] [name_of_the_screen] [path_to_output_folder]
```

Example:
```
python3 launch_parameter_screening.py parameter_table.csv dummy_screen /cluster/scratch/username
```

Please, make sure the section at the end of `launch_parameter_screening.py` is adapted to your need. By default the script generates jobs with a running time of 120h and 28GB of RAM which might be way more than what you need. 


```
sbatch -n 1 --cpus-per-task=4  --time=120:00:00  --mem-per-cpu=7168
```

<br>
<br>

## Collecting the results of a parameter screening

The results of the parameter screening are stored in the specified output folder. The results of each simulation are stored in a folder named after the simulation id. 
<br>
The script `collect_screening_data.py` can be used to create a summary of each simulation. It will create a folder named `summary` in the output folder and store the last time point of each simulation in it. This script can be called from the folder:

```
path/to/SimuCell3D/scripts/parameter_screening
```

In the following way:
```
python3 collect_screening_data.py [path_to_output_folder]
```
