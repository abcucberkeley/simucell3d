The python scripts contained in this folder can be debugged by running the following commands in the terminal:

```
cmake --no-warn-unused-cli -DCMAKE_BUILD_TYPE:STRING=Debug -DENABLE_PYTHON_BINDINGS=TRUE -DPYTHON_EXECUTABLE=$(which python3) -DCMAKE_EXPORT_COMPILE_COMMANDS:BOOL=TRUE -DCMAKE_C_COMPILER:FILEPATH=/usr/bin/gcc -DCMAKE_CXX_COMPILER:FILEPATH=/usr/bin/g++ -S/home/srunser/SimuCell3D_v2 -B/home/srunser/SimuCell3D_v2/build -G "Unix Makefiles"
```
```
cmake --build /home/srunser/SimuCell3D_v2/build --config Debug --target all -j 10 
```

Then, the following `launch.json` file should be created in the folder `/path/to/simucell/.vscode`

```
{
    "version": "0.2.0",
    "configurations": [
        {
            "name": "python_binding_debug",
            "type": "cppdbg",
            "request": "launch",
            "program": "${workspaceFolder}/scripts/python_bindings/bo_env/bin/python3",
            "stopAtEntry": false,
            "MIMode": "gdb",
            "miDebuggerPath": "/usr/bin/gdb",
            "cwd": "${workspaceFolder}/scripts/python_bindings",
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
The python script should be opened in the editor and the debug process can be started by pressing `F5`.

