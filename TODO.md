# TODO

[ ] Current
    [x] Save fibers and fiber surfaces.
    [x] Draw fibers and fiber surfaces at the same time.
    [x] Run on TQ
    [x] Histogram of sheet area/volume
    [x] Read and visualise .vtp atoms
    [x] Test correctness
    [x] Fix bug where you click inside a sheet
    [x] Add perturbation to the mesh
    [x] Add fiber point trace
    [x] Save/load reeb space
    [x] Why is it taking so much memory for example for ET-3?
    [x] Potential bug at isabel-4-cropped (crashed)



    [ ] (Optional) Map each regular vertex of the mesh to a sheet (sheet area estimate).
    [ ] (Optional) Optimize extracting only the singular fiber components (not all fiber components of singular fibers)
    [ ] (Optional) Refactor code to make it more usable 

    [x] Only save fiber triangle seeds to save memory
        Evaluating improvements on ET-3 (just RS stage Mb)
        Without seeds           - 807
        Seed fiber graphs       - 7607
        Seed triangles          - 838.41



[x] Faster Interactive Fiber Surfaces

    [x] Using TTK's fiber surface
        [x] Make sure the coordinates match in the range.
        [x] Read their fiber surface.
        [x] Remesh along singular fibers.
        [x] Extract regions based on the segmentation
            [x] Map non-singular vertices of FS triangles to the range
            [x] Extract singular fibers
            [x] Identify components and sheets
        [x] Find the sheet of each region.
        [x] Import TTK.
        [x] Use TTK's fiber surface algorithm.
        [x] Compute connected components
        [x] Start using CGAL's mesh operations and polygon soup operations.
        [x] Why is he CGAL mesh broken? 2 boundary edges for 2 spheres?

        [x] Optimisation
            [x] Speed up the remeshing with in-place Euler operations
            [x] Speed up remeshing with one pass over all isovalues
            [x] Optimise the fiber saturation (we don't need the full flexible fiber, just follow a few triangles)
            [x] Optimising using one single line lookup
            [x] Optimise modpoint computation

        [x] Bugs
            [x] Mini triangles (with area bellow 1e-14), add an epsilong to colouring the vertices, robustness issues.

    [ ] (Optional) Reimplement cropped marching tets.
    [ ] (Optional) Exact surface remeshing
    [ ] (Optional) Reeb graph of the restriction?

[x] Flexible fibers
    [x] Implement flexible fibers via a visible vertex
    [x] Implement flexible fibers via a line. 

    [x] Optimise the flexible fiber computation
        [x] You don't need a full search, it will be at one of the edges.
        [x] You don't need to do this for regular intersections, just used the +- triangles directly.
        [x] Test if this actually works

[ ] Stiched Fiber Surfaces
    [x] Implement the unaffected and regular cases
    [x] Implement the definite case
    [ ] Implement the indefinite case

    [x] Reduce number of triangles
        [x] Speed up the computation based on when triangles are in/out?
        [x] Can we look at this per tet/per triangle of the FS?

[ ] Sheet polygon extraction and shrinking
    [x] Extract the boundary half-edges of sheet polygons (with holes)
    [x] Extract the boundary points of the sheet polygons (with holes)
    [x] Shrink the polygons with  "2D Straight Skeleton and Polygon Offsetting"
    [ ] Test the shrinking, how do the holes behave?


[ ] Interactive UI
    [ ] Click in the range and a histogram of the sheets and some statistics about them.
    [ ] Click in the range and pull the feature fiber surfaces.



Plus of speeding up SFS - exactness, our method


Indeed the remeshing along the singular fiber doesn't need to be so exact, even if it's not i'll be a little off, maybe not even visible.
Then get the regions and try out a few fibers points.
Because we have the tet id, one fiber has one component per tet, so maybe we don't need the edge id, that could be fine.
I could also find the singular point and remesh around it manually.
This could actually work.
Also you could do the Delaneu? Or too complex actually...
