#pragma once

#include <vector>

#include "./FiberPoint.h"
#include "./TetMesh.h"
#include "./Arrangement.h"
#include "./ReebSpace.h"
#include "./ReebSpace2.h"

namespace fiber
{
    std::vector<FiberPoint> computeFiberFromFiberGraph(const TetMesh &tetMesh, Arrangement &arrangement, ReebSpace2 &reebSpace, const std::array<double, 2> &fiberPoint);

    std::vector<FiberPoint> computeFiber(const TetMesh &tetMesh, Arrangement &arrangement, ReebSpace &reebSpace, const std::array<double, 2> &fiberPoint, const int reebSheetIdOnly);


    std::vector<FiberPoint> computeFiberSAT(TetMesh &, Arrangement &, ReebSpace2 &, std::array<double, 2> );

    FiberGraph computeFiberGraph(TetMesh &tetMesh, Arrangement &singularArrangement, ReebSpace2 &reebSpace, std::array<double, 2> controlPoint);

    std::vector<FiberPoint> processFiberGraph(const TetMesh &tetMesh, Arrangement &arrangement, ReebSpace2 &reebSpace, const std::array<double, 2> &fiberPoint, FiberGraph &pg, const std::set<int> activeSheets);

    std::vector<FiberPoint> processFiberGraph2(const TetMesh &tetMesh, Arrangement &arrangement, ReebSpace2 &reebSpace, const std::array<double, 2> &fiberPoint, FiberGraph &pg, const std::set<int> activeSheets);

    std::vector<FiberPoint> computeFiberSurface(TetMesh &, Arrangement &, ReebSpace2 &, const std::vector<std::array<double, 2>> & );



    // Helper functions



    // Colour map for Reeb space sheets and fibers
    const std::vector<std::array<float, 3>> fiberColours = {
        {1.0f, 0.0f, 0.0f},    // Vivid Red
        {0.0f, 1.0f, 0.0f},    // Bright Green
        {0.0f, 0.0f, 1.0f},    // Pure Blue
        {1.0f, 1.0f, 0.0f},    // Bright Yellow
        {0.0f, 1.0f, 1.0f},    // Cyan
        {1.0f, 0.0f, 1.0f},    // Magenta
        {0.58f, 0.0f, 1.0f},   // Deep Purple
        {0.0f, 0.45f, 0.7f},   // Ocean Blue
        {1.0f, 0.5f, 0.0f},    // Orange
        {0.0f, 0.6f, 0.5f},    // Teal
        {1.0f, 0.84f, 0.0f},   // Gold
        {0.85f, 0.4f, 0.55f},  // Mauve
        {0.4f, 0.8f, 1.0f},    // Sky Blue
        {0.2f, 0.8f, 0.2f},    // Leaf Green
        {0.9f, 0.3f, 0.3f},    // Coral Red
        {0.6f, 0.6f, 0.0f}     // Olive
    };
};
