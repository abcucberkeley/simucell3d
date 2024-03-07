This folder contains script with various simulation routines that test the contact model: `contact_node_node_via_coupling`. 
This contact model creates mechanical link between the surfaces of adjeacent cells.

</br>

## Installation


```
cd path/to/simucell/build
```

```
cmake -DENABLE_PYTHON_BINDINGS=TRUE -DPYTHON_EXECUTABLE=$(which python3) -DCMAKE_BUILD_TYPE=Release .. && make -j6
```

```
cd path/to/simucell/test/test_contact_models/test_contact_node_node_via_coupling
```

```
python3 -m venv venv && source venv/bin/activate
```

```
pip3 install -r requirements.txt
```

</br>

## Running the tests

The scripts can then be run with the following commands:
```
python3 test_cell_doublet_angle/test_cell_doublet_angle.py
```

</br>

## Debugging the test in vscode
The tests can be debugged in vscode by first compiling the code with the following commands:
``` 
cmake --no-warn-unused-cli \
-DCMAKE_BUILD_TYPE:STRING=Debug \
-DENABLE_PYTHON_BINDINGS=TRUE \
-DPYTHON_EXECUTABLE=$(which python3) \
-DCMAKE_EXPORT_COMPILE_COMMANDS:BOOL=TRUE \
-DCMAKE_C_COMPILER:FILEPATH=/usr/bin/gcc \
-DCMAKE_CXX_COMPILER:FILEPATH=/usr/bin/g++ \
-S/home/srunser/SimuCell3D_v2 \
-B/home/srunser/SimuCell3D_v2/build \
-G "Unix Makefiles"
```
Make sure that the paths to the compiler and debugger are correct. Then call ```make``` with the following command:
```
cmake --build /home/srunser/SimuCell3D_v2/build --config Debug --target all -j 10 
```

Finally, create the following file: `path/to/simucell/.vscode/launch.json` with this content: 

```
{
    "version": "0.2.0",
    "configurations": [
        {
            "name": "python_binding_debug",
            "type": "cppdbg",
            "request": "launch",
            "program": "${workspaceFolder}/test/test_contact_models/test_contact_node_node_via_coupling/venv/bin/python3",
            "stopAtEntry": false,
            "MIMode": "gdb",
            "miDebuggerPath": "/usr/bin/gdb",
            "cwd": "${workspaceFolder}/test/test_contact_models/test_contact_node_node_via_coupling",
            "args": [
                "${file}"
            ],
            "launchCompleteCommand": "exec-run",
            "setupCommands": [
                {
                    "description": "Enable pretty-printing for gdb",
                    "text": "-enable-pretty-printing",
                    "ignoreFailures": true
                }
            ]
        },
    ]
}
```
Open the script that needs to be debugged and press F5. 



