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

        // Gives the coordinates of the point of intersection of the fiber and the triangles
        std::unordered_map<int, std::array<float, 3>> trianglePointCoordinates;

        // Grow a component component from a single seed pair
        static FiberComponent growFiberComponentFromSeedPair(const TetMesh& tetMesh, const std::array<double, 2>& controlPoint, const int triangleId, const int sheetId);

        // Grow a component component from a single seed pair to determine whether it contains a tet 
        static bool doesComponentContainTet(const TetMesh &tetMesh, const Arrangement &singularArrangement, const int seedTriangleId, const std::vector<int> &tetTriangleIds, const CartesianPoint &controlPoint);
};

class Fiber
{
    public:
        std::vector<FiberComponent> components;

        static Fiber computeLabeledFiber(TetMesh &tetMesh, Arrangement &singularArrangement, ReebSpace2 &reebSpace, std::array<double, 2> controlPoint, const std::set<int> &selectedSheetIds, const bool debugPrint = false);

        // Grow all the fiber components from the Seed set
        static Fiber growFiberFromSeedSet(const TetMesh &tetMesh, ReebSpace2 &reebSpace, const std::array<double, 2> &controlPoint, const std::vector<std::pair<int, int>> &fiberSeeds);

        // Grow fiber components from a seed set to determine which one contains a tet
        static int whichComponentContainsTet(const TetMesh &tetMesh, const Arrangement &singularArrangement, const std::vector<std::pair<int, int>> &fiberSeeds, const CartesianPoint &controlPoint, const int &tetId);

};


