#!/bin/bash

# Set up working folder
projectFolder="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
echo "The current folder is $projectFolder"

mkdir libraries

# VTK 
echo "Building VTK"

cd $projectFolder/libraries
git clone --recursive -b v9.4.1 https://gitlab.kitware.com/vtk/vtk.git vtk
cd vtk

mkdir build install
cd build

cmake -DCMAKE_INSTALL_PREFIX="$projectFolder/libraries/vtk/install" -DCMAKE_BUILD_TYPE="Release" ..
make -j 4
make install
    

## CGAL
echo "CGAL"

cd $projectFolder/libraries
git clone --recursive -b v6.0.1 https://github.com/CGAL/cgal cgal
cd ./cgal

mkdir build install
cd build

cmake -DCMAKE_INSTALL_PREFIX="$projectFolder/libraries/cgal/install" -DCMAKE_BUILD_TYPE="Release" ..
make -j 4
make install

## TTK
echo "TTK"

cd $projectFolder/libraries
git clone --recursive -b 1.3.0 https://github.com/topology-tool-kit/ttk ttk
cd ./ttk

mkdir build install
cd build

cmake -DCMAKE_INSTALL_PREFIX="$projectFolder/libraries/ttk/install" -DCMAKE_BUILD_TYPE="Release" -DTTK_BUILD_PARAVIEW_PLUGINS="Off" -DCMAKE_PREFIX_PATH="$projectFolder/libraries/VTK-9.4.1/install" ..
make -j 4
make install





