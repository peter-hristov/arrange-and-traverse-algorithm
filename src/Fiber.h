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
    std::vector<FiberPoint> computeFiberSAT(TetMesh &, Arrangement &, ReebSpace2 &, std::array<double, 2>, const std::set<int> &);
    std::vector<FiberPoint> computeFiberFromTriangleSeed(const TetMesh &tetMesh, Arrangement &arrangement, ReebSpace2 &reebSpace, const std::array<double, 2> &fiberPoint, const std::vector<std::pair<int, int>> &fiberSeeds);

    SurfaceMesh computeFiberSurfaceSingularSegment(TetMesh &, Arrangement &, ReebSpace2 &, const std::vector<std::array<double, 2>> &, const std::set<int> &selectedSheets = {});


    // Helper functions



    // Colour map for Reeb space sheets and fibers
    //const std::vector<std::array<float, 3>> fiberColours = {
    //{1.0f, 0.0f, 0.0f},    // Vivid Red
    //{0.0f, 1.0f, 0.0f},    // Bright Green
    //{0.0f, 0.0f, 1.0f},    // Pure Blue
    //{1.0f, 1.0f, 0.0f},    // Bright Yellow
    //{0.0f, 1.0f, 1.0f},    // Cyan
    //{1.0f, 0.0f, 1.0f},    // Magenta
    //{0.58f, 0.0f, 1.0f},   // Deep Purple
    //{0.0f, 0.45f, 0.7f},   // Ocean Blue
    //{1.0f, 0.5f, 0.0f},    // Orange
    //{0.0f, 0.6f, 0.5f},    // Teal
    //{1.0f, 0.84f, 0.0f},   // Gold
    //{0.85f, 0.4f, 0.55f},  // Mauve
    //{0.4f, 0.8f, 1.0f},    // Sky Blue
    //{0.2f, 0.8f, 0.2f},    // Leaf Green
    //{0.9f, 0.3f, 0.3f},    // Coral Red
    //{0.6f, 0.6f, 0.0f}     // Olive
    //};
    //const std::vector<std::array<float, 3>> fiberColours = {
    //{1.0f, 0.0f, 0.0f},    // Vivid Red
    //{0.0f, 1.0f, 0.0f},    // Bright Green
    //{0.0f, 0.0f, 1.0f},    // Pure Blue
    //{1.0f, 1.0f, 0.0f},    // Bright Yellow
    //{0.0f, 1.0f, 1.0f},    // Cyan
    //{1.0f, 0.0f, 1.0f},    // Magenta
    //{0.58f, 0.0f, 1.0f},   // Deep Purple
    //{0.0f, 0.45f, 0.7f},   // Ocean Blue
    //{1.0f, 0.5f, 0.0f},    // Orange
    //{0.0f, 0.6f, 0.5f},    // Teal
    //{1.0f, 0.84f, 0.0f},   // Gold
    //{0.85f, 0.4f, 0.55f},  // Mauve
    //{0.4f, 0.8f, 1.0f},    // Sky Blue
    //{0.2f, 0.8f, 0.2f},    // Leaf Green
    //{0.9f, 0.3f, 0.3f},    // Coral Red
    //{0.6f, 0.6f, 0.0f},    // Olive

    //// Additional 20 colors
    //{0.7f, 0.2f, 0.0f},    // Burnt Orange
    //{0.5f, 0.0f, 0.5f},    // Dark Magenta
    //{0.3f, 0.3f, 0.3f},    // Dark Gray
    //{0.8f, 0.6f, 0.4f},    // Tan
    //{0.6f, 0.2f, 0.2f},    // Brick Red
    //{0.0f, 0.5f, 0.0f},    // Forest Green
    //{0.0f, 0.3f, 0.5f},    // Deep Teal
    //{0.5f, 0.7f, 0.9f},    // Powder Blue
    //{0.9f, 0.7f, 0.8f},    // Pink
    //{0.4f, 0.0f, 0.2f},    // Plum
    //{0.3f, 0.6f, 0.1f},    // Grass Green
    //{0.0f, 0.7f, 0.7f},    // Aqua
    //{0.9f, 0.5f, 0.2f},    // Salmon
    //{0.7f, 0.5f, 0.9f},    // Lavender
    //{0.2f, 0.2f, 0.6f},    // Indigo
    //{0.6f, 0.4f, 0.0f},    // Bronze
    //{0.3f, 0.5f, 0.7f},    // Steel Blue
    //{0.8f, 0.3f, 0.6f},    // Fuchsia
    //{0.2f, 0.7f, 0.3f},    // Spring Green
    //{0.7f, 0.2f, 0.4f}     // Raspberry
    //};

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
