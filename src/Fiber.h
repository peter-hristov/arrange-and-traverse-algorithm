#pragma once

#include <vector>

#include "./TetMesh.h"
#include "./FiberPoint.h"
#include "./ReebSpace2.h"
#include "./Arrangement.h"
#include "./FiberSurface.h"





namespace fiber
{
    class FiberComponent
    {
        public:
            int sheetId;

            // All fiber segments in a fiber component
            std::vector<std::pair<int, int>> edges;

            // Gives the barycentric coordinates of the point of intersection of the fiber and the triangles
            std::unordered_map<int, std::array<float, 3>> triangleBarycentricCoordinates;

    };

    typedef std::vector<FiberComponent> Fiber;

    Fiber computeLabeledFiber(TetMesh &tetMesh, Arrangement &singularArrangement, ReebSpace2 &reebSpace, std::array<double, 2> controlPoint, const std::set<int> &selectedSheetIds);
    Fiber growSeedSet(const TetMesh &tetMesh, ReebSpace2 &reebSpace, const std::array<double, 2> &controlPoint, const std::vector<std::pair<int, int>> &fiberSeeds);

};
