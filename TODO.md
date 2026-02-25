# TODO



[ ] Faster Interactive Fibr Surfaces
    [ ] (Ideal case) Reimplement cropped marching tets.

    [ ] Using TTK's fiber surface
        [x] Make sure the coordinates match in the range.
        [x] Read their fiber surface.
        [ ] Extract singular fibers
        [ ] Remesh along singular fibers.
        [ ] Extract regions based on the segmentation
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
