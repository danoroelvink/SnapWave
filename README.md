# SnapWave
 Fast, implicit, unstructured-grid short wave solver
 
## Compiling and Installation

Executables are installed under `./snapwave/lnx64/bin` or `.\snapwave\x64\<configuration>\bin`.  

### Windows
Compiling is supported under VS with the Intel ifx compiler. Set the ONEAPI_ROOT environment variable before building `snapwave.sln`. The necessary dll's are collected as a post-build step based on the value of this variable.

### Linux
A `makefile` is provided with the repository. Compiling was tested with `ifort`, `ifx` and `gfortran`. Make sure your system netcdf libraries are compatible with the compiler used.
 
For example, to build under Ubuntu:
``` bash
sudo add-apt-repository universe
sudo apt update
sudo apt install libnetcdff-dev gfortran gcc
make clean
make
```
 
