# TODO



[ ] Faster Interactive Fibr Surfaces
    [ ] (Ideal case) Reimplement cropped marching tets.

    [ ] Using TTK's fiber surface
        [x] Make sure the coordinates match in the range.
        [x] Read their fiber surface.
        [x] Remesh along singular fibers.
        [ ] Extract regions based on the segmentation
            [ ] Map non-singular vertices of FS triangles to the range
            [ ] Extract singular fibers
            [ ] Identify components and sheets
        [ ] Find the sheet of each region

        [ ] Reeb graph of the restriction?

[ ] Stiched Fiber Surfaces
    [x] Implement the unaffected and regular cases
    [x] Implement the definite case
    [ ] Implement the indefinite case
    [ ] Speed up the computation based on when triangles are in/out?

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
