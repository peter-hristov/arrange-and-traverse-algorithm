#!/bin/bash

# Optionally read how many cores to use
CORES=${1:-$(nproc)}

# Set up working folder
projectFolder="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
echo "The current folder is $projectFolder"

# Set up where are going to clones and compile all dependencies
mkdir libraries

## Boost
echo "Building Boost"

cd $projectFolder/libraries
wget https://archives.boost.io/release/1.88.0/source/boost_1_88_0.tar.gz

tar -xvf boost_1_88_0.tar.gz
cd boost_1_88_0/
./bootstrap.sh --prefix="$projectFolder/libraries/install/boost_1_88_0"
./b2 install -j$CORES

## VTK 
echo "Building VTK"

cd $projectFolder/libraries
git clone --recursive -b v9.4.1 https://gitlab.kitware.com/vtk/vtk.git vtk
cd vtk

mkdir build install
cd build

cmake -DCMAKE_INSTALL_PREFIX="$projectFolder/libraries/install/vtk" -DCMAKE_BUILD_TYPE="Release" ..
make -j$CORES
make install
    

### CGAL
echo "CGAL"

cd $projectFolder/libraries
git clone --recursive -b v6.0.1 https://github.com/CGAL/cgal cgal
cd ./cgal

mkdir build install
cd build

cmake -DCMAKE_PREFIX_PATH="$projectFolder/libraries/install/boost_1_88_0" -DCMAKE_INSTALL_PREFIX="$projectFolder/libraries/install/cgal" -DCMAKE_BUILD_TYPE="Release" ..
make -j$CORES
make install

### TTK
echo "TTK"

cd $projectFolder/libraries
git clone --recursive -b 1.3.0 https://github.com/topology-tool-kit/ttk ttk
cd ./ttk

mkdir build install
cd build

cmake -DCMAKE_INSTALL_PREFIX="$projectFolder/libraries/install/ttk" -DCMAKE_BUILD_TYPE="Release" -DTTK_BUILD_PARAVIEW_PLUGINS="Off" -DCMAKE_PREFIX_PATH="$projectFolder/libraries/install/vtk;$projectFolder/libraries/install/boost_1_88_0" ..
make -j$CORES
make install

### RS Explorer
cd $projectFolder/build
cmake -DCMAKE_PREFIX_PATH="$projectFolder/libraries/install/cgal;$projectFolder/libraries/install/vtk;$projectFolder/libraries/install/ttk" -DCMAKE_EXPORT_COMPILE_COMMANDS=On -DCMAKE_BUILD_TYPE=Release ..
make

cd $projectFolder






