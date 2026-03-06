# TODO

[ ] Other ideas
    [x] Save fibers and fiber surfaces.
    [x] Draw fibers and fiber surfaces at the same time.

    [ ] For each FS, histogram of triangles and trea per sheet.
    [ ] List of intersected sheets and their range area.
    [ ] Map each regular vertex of the mesh to a sheet (sheet area estimate).

[ ] Faster Interactive Fiber Surfaces
    [^] (Optional) Reimplement cropped marching tets.

    [ ] Using TTK's fiber surface
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

        [ ] Optimisation
            [x] Speed up the remeshing with in-place Euler operations
            [x] Speed up remeshing with one pass over all isovalues
            [x] Optimise the fiber saturation (we don't need the full flexible fiber, just follow a few triangles)
--->        [ ] Optimising using one single line lookup
            [ ] Optimise modpoint computation?

        [x] Bugs
            [x] Mini triangles (with area bellow 1e-14), add an epsilong to colouring the vertices, robustness issues.

        [ ] Refactor code to make it more usable 

        [ ] (Optional) Reeb graph of the restriction?

[ ] Stiched Fiber Surfaces
    [x] Implement the unaffected and regular cases
    [x] Implement the definite case
    [ ] Implement the indefinite case

    [x] Reduce number of triangles
        [x] Speed up the computation based on when triangles are in/out?
        [x] Can we look at this per tet/per triangle of the FS?


[ ] Flexible fibers
    [x] Implement flexible fibers via a visible vertex
    [x] Implement flexible fibers via a line. 

    [ ] Optimise the flexible fiber computation
        [x] You don't need a full search, it will be at one of the edges.
        [x] You don't need to do this for regular intersections, just used the +- triangles directly.
        [ ] Test if this actually works

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
