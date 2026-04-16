#pragma once

#include <vector>

#include "./TetMesh.h"
#include "./ReebSpace2.h"
#include "./Arrangement.h"
#include "./FiberSurface.h"

class FiberComponent
{
    public:
        int sheetId;

        // All fiber segments in a fiber component
        std::vector<std::pair<int, int>> edges;

        // Gives the barycentric coordinates of the point of intersection of the fiber and the triangles
        std::unordered_map<int, std::array<float, 3>> triangleBarycentricCoordinates;

        // (Helper) Cheap test if a point is in a triangle, then if it is compute bary coords and them actual point coords
        static std::optional<std::array<float, 3>> tryComputePointCoordinates(const TetMesh& tetMesh, const int triangleId, const CartesianPoint& P);

        // Grow a component component from a single seed pair
        static FiberComponent growFiberComponentFromSeedPair(const TetMesh& tetMesh, const std::array<double, 2>& controlPoint, const int triangleId, const int sheetId);
};

class Fiber
{
    public:
        std::vector<FiberComponent> components;

        static Fiber computeLabeledFiber(TetMesh &tetMesh, Arrangement &singularArrangement, ReebSpace2 &reebSpace, std::array<double, 2> controlPoint, const std::set<int> &selectedSheetIds);

        // Grow all the fiber components from the Seed set
        static Fiber growFiberFromSeedSet(const TetMesh &tetMesh, ReebSpace2 &reebSpace, const std::array<double, 2> &controlPoint, const std::vector<std::pair<int, int>> &fiberSeeds);
};




