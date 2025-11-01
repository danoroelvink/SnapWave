@echo off
setlocal enabledelayedexpansion
set t0=%time%
del /f INPUT PRINT
del /f *.nc

set "extension=.swn"
set OMP_NUM_THREADS_SWAN=20

for %%f in ("*%extension%") do (
    set "filename=%%~nf"
    echo Running executable with argument: !filename!
    copy /a /y  !filename!!extension! INPUT
)
swan 
rem call D:\OSS\Delft3d_mordev_superbranch\install_all\bin\swan.bat
