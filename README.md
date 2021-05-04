## NuSolarFitter
Early version of the code to fit solar neutrinos with JUNO
Based on the vector-implementation in BxRandom

## Installation

### Requirements
gcc >= 9.3.0
cmake >= 3.4.1
ROOT >= 6.22

Because the software is written with the c++14 standard, the compilers must be recent enough to reflect this. Additionally, ROOT must be installed against this standard.

### @brief
Because the CMakeLists operates on environment variables, the only thing to do before installing NuSolarFitter is to source the correct compilers and "thisroot.sh" (which is probably in your .bashrc anyway). As long as "gcc", "g++", and "ROOTSYS" exist, the software should compile.

## Usage
The program can be compiled using RebuildAndInstall_sh.sh, which automatically builds/install the programm according to the CMakeLists.txt while using the correct compilers. If you are calling the script from outside of the root folder, you can set the SIM_DIR environment variable to where you would like the build/ and install/ folders containing binaries etc. to go.
Alternatively, one can run the program interactively using root and .L
