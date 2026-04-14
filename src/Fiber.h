#pragma once

#include <vector>

#include "./FiberPoint.h"
#include "./TetMesh.h"
#include "./Arrangement.h"
#include "./ReebSpace.h"
#include "./ReebSpace2.h"
#include "./SurfaceMesh.h"

namespace fiber
{
    std::vector<FiberPoint> computeLabeledFiber(TetMesh &, Arrangement &, ReebSpace2 &, std::array<double, 2>, const std::set<int> &);
    std::vector<FiberPoint> growSeedSet(const TetMesh &tetMesh, Arrangement &arrangement, ReebSpace2 &reebSpace, const std::array<double, 2> &fiberPoint, const std::vector<std::pair<int, int>> &fiberSeeds);

    SurfaceMesh computeSegmentedFiberSurface(TetMesh &, Arrangement &, ReebSpace2 &, const std::vector<std::array<double, 2>> &, const std::set<int> &selectedSheets = {});


    // ChatGPT generate (supposedly based on Tableu)
    //
    const std::vector<std::array<float, 3>> fiberColours = {
        {0.121f, 0.466f, 0.705f}, // blue
        {1.000f, 0.498f, 0.054f}, // orange
        {0.173f, 0.627f, 0.173f}, // green
        {0.839f, 0.153f, 0.157f}, // red
        {0.580f, 0.404f, 0.741f}, // purple
        {0.549f, 0.337f, 0.294f}, // brown
        {0.890f, 0.467f, 0.761f}, // pink
        //{0.121f, 0.466f, 0.705f}, // blue
        {0.498f, 0.498f, 0.498f}, // gray
        {0.737f, 0.741f, 0.133f}, // olive
        {0.090f, 0.745f, 0.811f}, // cyan
        {0.682f, 0.780f, 0.909f}, // light blue
        {1.000f, 0.733f, 0.471f}, // light orange
        {0.596f, 0.875f, 0.541f}, // light green
        {1.000f, 0.596f, 0.588f}, // light red
        {0.773f, 0.690f, 0.835f}, // light purple
        {0.769f, 0.612f, 0.580f}, // light brown
        {0.969f, 0.714f, 0.824f}, // light pink
        {0.780f, 0.780f, 0.780f}, // light gray
        {0.859f, 0.859f, 0.553f}, // light olive
        {0.620f, 0.855f, 0.898f}  // light cyan
    };



};
