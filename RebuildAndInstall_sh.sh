#!/bin/sh
################################################################################
# Script to start the CMake-driven rebuild process for the entire 
# project. The project will be rebuilt in a directory 'build'; any preexisting 
# directory 'build' will be deleted. 
################################################################################

# define partent directory for the simulation build and install
SIM_DIR=/home/alexandregoettel/work/PhD/Juno/JunoFit/skeleton

# get script name and call directory
 ME="$(basename "$(readlink -f "$BASH_SOURCE")")"
 ME="[$ME]"
 FROM=$PWD
 
 echo "$ME Starting to build project..."

# get name of project's base directory by taking the parent directory of the
# parent directory of this script (assuming it is in ".../<BASE_DIR>/scripts/)
 BASE_DIR="$(dirname "$(dirname "$(readlink -f "$BASH_SOURCE")")")"
 echo "$ME Project base directory is '$BASE_DIR'".

# replace existing 'build' directory in the project's base directory
 echo "$ME Going to directory where to build the project."
 cd $SIM_DIR
 echo "$ME Removing existing 'build' directory."
 rm -rf build
 echo "$ME Removing existing 'install' directory."
 rm -rf install
 echo "$ME Creating new 'build' directory."
 mkdir build
 echo "$ME Creating new 'install' directory."
 mkdir install
 echo "$ME Going into 'build' directory."
 cd build
 
# build with CMAKE; 
 echo "$ME Calling CMake to build."
 cmake -DCMAKE_INSTALL_PREFIX=$SIM_DIR/install \
       -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ \
       $SIM_DIR
 cmake --build . -j 4

# install with CMAKE
echo "$ME Calling CMake to install."
  cmake --build . --target install
   
# go back to call directory
 echo "$ME Going back to call directory."
 cd $FROM
