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
    std::vector<FiberPoint> computeFiberFromFiberGraph(const TetMesh &tetMesh, Arrangement &arrangement, ReebSpace2 &reebSpace, const std::array<double, 2> &fiberPoint);

    std::vector<FiberPoint> computeFiber(const TetMesh &tetMesh, Arrangement &arrangement, ReebSpace &reebSpace, const std::array<double, 2> &fiberPoint, const int reebSheetIdOnly);


    std::vector<FiberPoint> computeFiberSAT(TetMesh &, Arrangement &, ReebSpace2 &, std::array<double, 2>, const std::set<int> &);
    std::vector<FiberPoint> computeFiberFromTriangleSeed(const TetMesh &tetMesh, Arrangement &arrangement, ReebSpace2 &reebSpace, const std::array<double, 2> &fiberPoint, const std::vector<std::pair<int, int>> &fiberSeeds);

    std::vector<FiberPoint> processFiberGraph(const TetMesh &tetMesh, Arrangement &arrangement, ReebSpace2 &reebSpace, const std::array<double, 2> &fiberPoint, FiberGraph &pg, const std::set<int> activeSheets);
    std::vector<FiberPoint> processFiberGraph2(const TetMesh &tetMesh, Arrangement &arrangement, ReebSpace2 &reebSpace, const std::array<double, 2> &fiberPoint, FiberGraph &pg, const std::set<int> activeSheets);

    std::vector<FiberPoint> computeFiberSurface(TetMesh &, Arrangement &, ReebSpace2 &, const std::vector<std::array<double, 2>> & , int);

    SurfaceMesh computeFiberSurfaceSingularSegment(TetMesh &, Arrangement &, ReebSpace2 &, const std::vector<std::array<double, 2>> & , int);

    std::vector<FiberPoint> computeFiberPointsFromSurfaceMesh(SurfaceMesh &, ReebSpace2 &, const std::set<int> &);


    std::vector<int> extractPath(const int &start, const std::unordered_map<int, std::vector<int>> &fgAdj, std::vector<bool> &visited);
    std::vector<int> extractCycle(const int &start, const std::unordered_map<int, std::vector<int>> &fgAdj, std::vector<bool> &visited);
    std::pair<std::map<int, std::vector<int>>, std::map<int, std::vector<int>>> buildFiberGraphPathsAndCycles(const TetMesh &tetMesh, ReebSpace2 &reebSpace, FiberGraph fg);



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
        {0.894f, 0.102f, 0.110f}, // red
        {0.216f, 0.494f, 0.722f}, // blue
        {0.596f, 0.306f, 0.639f}, // purple
        {1.000f, 0.498f, 0.000f}, // orange
        {1.000f, 1.000f, 0.200f}, // yellow
        {0.651f, 0.337f, 0.157f}, // brown
        {0.969f, 0.506f, 0.749f}, // pink
        {0.600f, 0.600f, 0.600f}, // gray
        {0.000f, 0.745f, 0.933f}, // cyan
        {0.941f, 0.894f, 0.259f}, // gold
        {0.745f, 0.157f, 0.819f}, // magenta
        {0.255f, 0.412f, 0.882f}, // royal blue
        {0.180f, 0.800f, 0.443f}, // lime green
        {0.902f, 0.647f, 0.294f}, // tangerine
        {0.510f, 0.216f, 0.384f}, // dark purple
        {0.376f, 0.682f, 0.776f}, // teal
        {0.961f, 0.737f, 0.369f}, // light gold
        {0.851f, 0.420f, 0.090f}, // burnt orange
        {0.722f, 0.451f, 0.200f},  // sienna
        {0.302f, 0.686f, 0.290f} // green
    };



};
