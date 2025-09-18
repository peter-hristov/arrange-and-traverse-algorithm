# Optimisation 13.09.2025
Unit Test the Final reeb space

Make an array isPseudoSingular
Read VTU and directly downsample myself
Optimise memory usage by my tetMesh class

Done - Debug nested faces
Done - Optimise the red-blue computation.
    Done - Email Sebastien 
    Done - Parallelise
Done - Computing the unique components in a preimage graph (for the correspondence graph)
    Done - track the current components without O(fiber size)
    Done - You only need them once per face.
    Done - Use IDs in CGAL directly (get rid of all the maps).
Done - Use id's instead of maps with half-edge and faces.

# Refactor

Done - IO
Done - TetMesh
Done - Arrangement
Done - ReebSpace
Done - Fiber computation
Done - Move the hashes outside of Disjoint set
Done - PlotWidet (Simplify the 2D visualiser)
Done - TracerVisualiserWidget


Doen - (Later) Update the timer for use a hashmap
Done - (Later) Put the Reeb space computation in a different library
Done - (Nope) Is this called twice? - ReebSpace::determineCorrespondence(tetMesh, arrangement, curr); in the BFS?
Done - (Not for now) Switch to using size_t for indices 
    Not as vital for now, it won't be soon that we can compute a Reeb space for a tet with over 32bit vertices.


# Optomise memory:
done - Remove the polygons and colours from PlotWidget
done - remove  - data->arrangementFiberComponents.resize(data->arrangementFacesIdices.size(), -1);

Use direct pointers in uion find? Don't think it could work though



# Tasks
Done - Compute the arrangement of a number of line segments
Done - Traverse the half-edge arrangement data struture (use ccb and twin)
Done - How to see where the segment came from? Use with_history and originating curve
Done - Set up IDs for the faces (use a map for now, store ID field in data later)
Done - Get a running BFS
Done - Lower and upper star for each edge (loop all tets, maybe we don't need a big datastruture)
Done - Compute fibers in each face
Done - Compute preimage graph for each face
Done - Compute connectivity of fibers for each face (via disjoint set)
Done - Compute fiber fiber components assignemtns (make sure differnt sheets show the same)
Done - Compute the Reeb space
Done - Colour the sheets in the range
Done - Flexible fibers

Use better startegy to find fiber components.
Use better maps between various handles and indices and so on.



Later
- Can you tell the algo which segments have the same point? Does that help?
- How do we use rational numbers for precision? Robustness.






















### TODO

Resurrect volume projection.

Write paper.






### DONE
Use OpenMP.
Fix projected triangles when changin u,v fields
Cache histogram and scatterplot.
Make sure you have a consistent logic about deleting FSCP points and you don't do too much work.
Toggles for FS/IS visibility and Union/Intersection.
Add Carthesian grid and 0 lines.
Update object selection.





### Backlog

Do the u and v scalar fields.
Crop & downsize domain.
Recompute min/max values for every field
Change isovalueSpinBox on isofield change



Deal with raw pointers in Data.h






### Key New Features

# Most Important for Now
Once you click an object it is locked and shown

Click objects show only them on the plot accross different timesteps

Selecting one object (via clicking) 
Single object evolution.

Export dat


# Tracking Over Time
Get triangles as output.

# Eventually
Exporting Features

# Bugs
Opacity does not work well.

Isotime



### Crasy Ideas
Abstract Fiber surface (distance functions over other things)






### DONE

# 1st Iteration 14.01.2019 - 18.01.2019
    Compute Isosurfaces
    Produce 2D Histogram
    Compute Fiber Sufaces
    Add interactivity

# 2nd Iteration 28.01.2019 - 01.02.2019
    Clockwise Convex FSCP
    Arcball Camera
    Colourful Scatterplot
    Read NC data files
    Get only big components
    Colourful Scatterplot with mean and standard deviation
    Put up axis and units
    Scale Computation

# 3rd Iteration 04.02.2019 - 08.02.2019
    Added Downsampling for Isosurfaces and Fiber Surfaces.
    Add volumes and exporting to fiber surfaces and isosurfaces.
    Scale computing volumes on the large dataset.
    Add an isovalue histogram
    Filter components by size
    Add slider for changing connected component minmum size
    Add static clipping for dimensions via command line parameters. 
    Make volumes work with downsampling.
    Make it work for any Convex FSCP

# 4rd Iteration 18.02.2019 - 26.02.2019
    Worked on the continuous scatterplot


# 5rd Iteration 04.03.2019 - 08.03.2019
    Reset FS when no points are there
    Add a checkbox to toggel FS/IS
    Add IS flip and plot
    Fix recompute on IS number box.
    Add a usage string
    Added interpolation to FS computation - none, nearest and bilinear
    Add an option to change the scale of the scatterplot
    Add projecting lines instead of points
    Fix Range clipping on the histogram

# 6th Iteration
    Finally Finished Working on the Continuous Scatterplot

### TODO

# Vital

    Figure out how to distribute application
    Refactor scatterplot rendering

    Contour Tree / Merge Tree and easy location of components

    Interface
        Add isovalue labels



    Refactor Code
    Add Documentation
    Call from python scripts using pybind.

# Optional
    Improve interface
    Any FSCP
    Verify correctness using a plot of f and f'
    Click on 3D plot to select objects (rays?)
    Add clang tools to static analyse and format the code.
    Save only the biggest objects.

# Scrapped
    Implement Central Differencing and add a function/gradient plot.
    Add volume size to mask file
    Axis on the 3D plot.
    Add dynamic clipping for dimensions via the qt interface.
    Simplify by removing points or flattening the distribution for the isosurface scalar field.


### NOTES

# Optimisation
    Compute Scatterplot Faster -> to do with weird mem access
    Volumes Take 10s to compute, figure out why. Is it memory access again?
    I got over a second speedup for just not processing unused cubes;
    There is a MAJOR problem with distance field computation. Over 12s on full data. I cannot get it down significatly more with this method.

    Quite alot of time is spent on just reading from the 2D arrays in the PlotWidget for example in getting the FS distance values and computing the histogram

# Things which are slow and I cannot easily solve:

    Scatterplot                     - memory access is slow
    Volumes                         - not sure, probably related to memory access
    Distance Field                  - voronoi diagram GPU z-buffer algorithm
    Isosurface                      - (flexible and octrees and paralle)

# Notes on Optimisation
    Skipping empty triangles gave a 0.9s speedup.
    Went from 0.9s to 2.6s when we go from 0 triangles to 1,587,894 triangles.
