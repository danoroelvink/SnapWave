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
 
To automatically remove module and object files, run with 
``` bash
make clean
make STRIP_BUILD=1
```
## Log files

SnapWave writes its log to the run directory (modelled on XBeach):

| File | Contents |
|---|---|
| `snapwave_log.txt` | Start/end time, every input parameter (defaults marked `(default)`), grid and boundary summaries, one line per time step with solver convergence, total run time and number of warnings |
| `snapwave_warning.txt` | Only created when there are warnings, e.g. unknown or misspelled keywords in `snapwave.inp`, lines without `=`, boundary series not covering the simulation period, solver not converged |
| `snapwave_error.txt` | Only created when SnapWave stops on an input error, with the reason, e.g. a missing input file, an empty `bndfile`, or boundary files whose number of columns does not match the number of boundary points |

The input parameters (including defaults) are also stored in the map and his files, as attributes of the scalar variable `inp`; the attribute `parameters_with_default_values` lists the keywords that were not given in `snapwave.inp`.

In the code, use `writelog(destination, message)` from `snapwave_log` instead of `write(*,*)`. The destination is any combination of `l` (log file), `s` (screen), `w` (warning file) and `e` (error file). Use `writewarning(message)` for warnings and `halt_program(message)` to stop on an error. `num2str()` converts integers and reals to strings, for example `call writelog('ls', 'nodes: '//num2str(no_nodes))`.
