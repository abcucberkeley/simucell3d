# First time installation 

Load the required modules
```
env2lmod
module load gcc/9.3.0 openblas/0.3.20 python/3.11.2 cmake/3.25.0
```

Create a virtual environment
```
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


# Classic usage 
Load the required modules
```
env2lmod
module load gcc/9.3.0 openblas/0.3.20 python/3.11.2 cmake/3.25.0
```

Activate the virtual env
```
source ./venv_screening/bin/activate
```

