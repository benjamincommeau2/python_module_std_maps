# Unfinished project. Currently under development. Next step is to modify StdMapPython.cpp to create its own type of python variable as std::map if that is possible.
# Creating a python module that allows the user to store numpy arrays in C++ standard library maps.
# Currently works for default library path installation for python3.6 and numpy1.17. Changing the header library calls for these versions may be required in the C++ python extender script (e.i. StdMapPythonC.cpp and StdMapPython.h) may be required.
```console
foo@bar:~$ bash install_StdMapPythonC_in_python.sh
```
# Compile Python C++ extender.
```console
foo@bar:~$ python3 test_StdMapPythonC.py
```
# Run run python3 test script. Calls multiple python scripts called dtype_StdMapPython.cpp passing different numpy value types and tests if the try catch block in the StdMapPython.cpp catches any errors and displays them in the terminal.
