#pragma once

#include "./TetMesh.h"
#include "./Arrangement.h"
#include "./ReebSpace2.h"
#include "./FiberPoint.h"

// This functionality is largely depricated, related to the idea that we can form a fiber surface by stitching fiber togehter.
// This is also not complete, this only does the fiber stitching for regular and define edges, no indefinite
namespace fiber::stitching
{
    std::vector<int> extractPath(const int &start, const std::unordered_map<int, std::vector<int>> &fgAdj, std::vector<bool> &visited);
    std::vector<int> extractCycle(const int &start, const std::unordered_map<int, std::vector<int>> &fgAdj, std::vector<bool> &visited);
    std::pair<std::map<int, std::vector<int>>, std::map<int, std::vector<int>>> buildFiberGraphPathsAndCycles(const TetMesh &tetMesh, ReebSpace2 &reebSpace, FiberGraph fg);


    std::tuple<std::vector<int>, std::vector<int>, bool> getActiveTrianglesInPath(const std::vector<int> &pathA, const std::vector<int> &pathB, const std::unordered_set<int> &minusTrianglesSet, const std::unordered_set<int> &plusTrianglesSet);
    std::tuple<std::vector<int>, std::vector<int>, bool> getActiveTrianglesInCycle(const std::vector<int> &cycleA, const std::vector<int> &cycleB, const std::unordered_set<int> &minusTrianglesSet, const std::unordered_set<int> &plusTrianglesSet);

    std::vector<FiberPoint> computeTrianglesBetweenCorrespondingPaths(const std::vector<int> &pathA, const std::vector<int> &pathB, const std::unordered_set<int> &minusTrianglesSet, const std::unordered_set<int> &plusTrianglesSet, const std::vector<std::unordered_map<int, std::array<double, 3>>> &barycentricCoordinatesPerTriangle, const TetMesh &tetMesh, const int sheetId, const std::array<float, 3> &sheetColour, const int i, const std::array<float, 3> &edgePointDomain);

    std::vector<FiberPoint> computeTrianglesBetweenCorrespondingCycles(const std::vector<int> &cycleA, const std::vector<int> &cycleB, const std::unordered_set<int> &minusTrianglesSet, const std::unordered_set<int> &plusTrianglesSet, const std::vector<std::unordered_map<int, std::array<double, 3>>> &barycentricCoordinatesPerTriangle, const TetMesh &tetMesh, const int sheetId, const std::array<float, 3> &sheetColour, const int i, const std::array<float, 3> &edgePointDomain);

    std::array<double, 3> computeBarycentricCoordinates(const TetMesh &tetMesh, const int &triangleId, const std::array<double, 2> &fiberPoint);



    std::vector<FiberPoint> computeFiberSurface(TetMesh &tetMesh, Arrangement &singularArrangement, ReebSpace2 &reebSpace, const std::vector<std::array<double, 2>> &controlPoints, int _sheetId);
}
