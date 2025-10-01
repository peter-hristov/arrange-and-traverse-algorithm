#!/bin/bash

# Build dependencies
projectFolder="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
echo "The current folder is $projectFolder"

rm -rf $projectFolder/librarires
mkdir $projectFolder/librarires

# VTK 
echo "Building VTK"

cd $projectFolder/librarires
git clone --recursive https://gitlab.kitware.com/vtk/vtk.git VTK-9.4.1
cd ./VTK-9.4.1
git checkout v9.4.1 

mkdir build install
cd build

cmake -DCMAKE_INSTALL_PREFIX="$projectFolder/librarires/VTK-9.4.1/install" -DCMAKE_BUILD_TYPE="Release" ..
make -j 4
make install
    

## CGAL
echo "CGAL"

cd $projectFolder/librarires
git clone --recursive https://github.com/CGAL/cgal cgal
cd ./cgal
git checkout v6.0.1 

mkdir build install
cd build

cmake -DCMAKE_INSTALL_PREFIX="$projectFolder/librarires/cgal/install" -DCMAKE_BUILD_TYPE="Release" ..
make -j 4
make install

# Build application
echo "Reeb Space Project"

cd $projectFolder
cd build

cmake -DCMAKE_PREFIX_PATH="$projectFolder/librarires/cgal/install;$projectFolder/librarires/VTK-9.4.1/install" -DCMAKE_BUILD_TYPE="Release" ..
make -j 4


# Try this on HPC Tetralight
#cmake -DCMAKE_PREFIX_PATH="/proj/reeb-space-storage/users/x_pethr/librarires/cgal-6.0.1/install;/proj/reeb-space-storage/users/x_pethr/librarires/vtk-9.4.1/install;/proj/reeb-space-storage/users/x_pethr/librarires/boost_1_88_0/install" -DCMAKE_BUILD_TYPE="Release" ..
