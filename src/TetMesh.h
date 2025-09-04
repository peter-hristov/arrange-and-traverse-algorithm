#pragma once

#include <map>
#include <set>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>

#include "./DisjointSet.h"
#include "./Hashing.h"

class TetMesh
{
  public:
    TetMesh() {}

    // Domain and range coordinates
    std::vector<float> vertexCoordinatesF;
    std::vector<float> vertexCoordinatesG;
    std::vector<std::array<float, 3>> vertexDomainCoordinates;

    // Bounding box min/max for the domain and range coordinates of all vertices
    float minF, maxF, minG, maxG;
    float minX, maxX, minY, maxY, minZ, maxZ;

    // (Optional) indicative names for the two scalar fields and their units
    std::string longnameF, longnameG, units;

    // Combinatorial structure of the mesh
    std::vector<std::array<int, 4>> tetrahedra;
    // We always assume that edge vertices are in sorted order (but index)
    std::vector<std::array<int, 2>> edges;
    std::unordered_map<std::array<int, 2>, int, MyHash<std::array<int, 2>>> edgeIndices;
    std::unordered_map<std::array<int, 2>, int, MyHash<std::array<int, 2>>> edgeSingularTypes;

    std::map<std::array<int, 2>, std::set<int>> upperLink;
    std::map<std::array<int, 2>, std::set<int>> lowerLink;
    std::map<std::array<int, 2>, std::vector<int>> upperStarTriangles;
    std::map<std::array<int, 2>, std::vector<int>> lowerStarTriangles;
    std::vector<std::set<int>> triangles;
    std::unordered_map<std::set<int>, int, MyHash<std::set<int>>> triangleIndices;

    //std::vector<std::vector<int>> tetIncidentTriangles;
    std::vector<std::set<int>> tetIncidentTriangles;

    int singularEdgesNumber;
    int regularEdgesNumber;


    void computeBoundingBoxes();
    void computeDomainBoundingBox();
    void computeRangeBoundingBox();

    void sortVertices();
    void printMesh();
    void perturbRangeValues(const float &perturbationEpsilon);


    // We only add upperLink/lowerLink to data, the rest of data is unchagned
    void computeUpperLowerLinkAndStar();
    void computeSingularEdgeTypes();

    // From the tet soup, get the edges and triangles and tet-incidence of the triangles
    void computeCombinatorialStructure();

    // Give the edge (aIndex, bIndex), is the vertex vIndex from its link in the upper and lower link of the edge
    // We assume that aIndex < bIndex for consistent orientation.
    bool isUpperLinkEdgeVertex(int, int, int);
};
