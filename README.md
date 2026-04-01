# Introduction

Welcome to the RS explorer project folder. The following build instructions have been tested on Ubuntu 22.04 and Ubuntu 24.04. For other Linux distrbutions almost the same steps should work. This project has been build on top the github source code for the arrange and traverse algorithm[1].


This application has the following dependencies:
    VTK  v9.4.1
    TTK  v1.3.0
    CGAL v6.0.1

# Building 
You could install all the dependencies on your own, or use the build script we provide. Our script clones vtk, ttk and cgal into the ./libraries folder and then compiles and intall them in the folder ./libraries/install. 
The ttk lbraris linked to the vtk install in ./libraries/install.
To run the build script with 4 cores:

``` 
bash build.sh 4
```

For using <n> cores run:
```
bash build.sh <n>

```

Our build script also build the RS Visualiser application. If you have simillar version of the dependencies you can build the RS visualiser youself by with cmake:

```
cmake -DCMAKE_PREFIX_PATH="<path_to_cgal_install>;<path_to_vtk_install>;<path_to_ttk_install>" -DCMAKE_BUILD_TYPE=Release ..
```

# Running

To run our application we have provided a number of tests datasets. The simples one is three-sheet-toy.vtu. It has been used in previous Reeb space papers [1,2]. Run with:

```
./build/rsX -f ./data/three-sheet-toy.vtu
```

To explore other features of our application run:

```
./build/rsX -h
```

We have also provided two of the datasets we have used in the paper, torus and ethanediol, as well as their downsampled version, which are faster to compute. Refer to Table A1 in the paper for computation times and memory usage.

You have the option to save a compute Reeb space, so that you can load it on a rerun, which is much faster than recomputing. For example.

Save a Reeb space:
```
./build/rsX -f ./data/torus/downsample-id-1.vtu -s ./data/torus/downsample-id-1.rs
```

Load a Reeb space
```
./build/rsX -f ./data/torus/downsample-id-1.vtu -l ./data/torus/downsample-id-1.rs
```

Note that the .rs file extension is just a matter of convention.

# Controls

## Domain view
Left click and drag             - rotate camera
Right click                     - select a colour from the segmented fiber surface and select the sheet that corresponds to it

## Range view
Left click                      - compute labeled fiber, draw the mouse to compute multiple (enable Trace fiber from the bottom controls to leave a trace)
Right click                     - add a point to the fiber surface control polygon
Backspace                       - remove the last point added to the fiber surface control polygon
Shift + left click              - select all sheets at that point in the range

## Bottom panels
The information shows the selected sheets.
The visibility panel allows the user to show/hide the rendered geometry as well as set opacity.
The clear panel clear rendered geometry.
The compute panels computes things like fiber surfaces (FS), sheet-features, fiber surface from a fiber trace and the option to trace fibers continuously.
The sheet selection panel allows the selection of sheets by ID or by top <n>.


Happy Reeb space exploring!

1. https://github.com/peter-hristov/arrange-and-traverse-algorithm
2. Hristov, P., Sakurai, D., Carr, H., Hotz, I. and Masood, T.B., 2025, August. Arrange and Traverse Algorithm for Computation of Reeb Spaces of Piecewise Linear Maps. In Computer Graphics Forum (Vol. 44, No. 5, p. e70206).
3. Tierny, J. and Carr, H., 2016. Jacobi fiber surfaces for bivariate Reeb space computation. IEEE Transactions on Visualization and Computer Graphics, 23(1), pp.960-969.

