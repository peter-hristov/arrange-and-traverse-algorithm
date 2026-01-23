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
    std::vector<double> vertexCoordinatesF;
    std::vector<double> vertexCoordinatesG;
    std::vector<std::array<float, 3>> vertexDomainCoordinates;

    // Bounding box min/max for the domain and range coordinates of all vertices
    double minF, maxF, minG, maxG;
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

    std::vector<std::vector<int>> upperStarTrianglesNew;
    std::vector<std::vector<int>> lowerStarTrianglesNew;

    std::vector<std::set<int>> triangles;
    std::unordered_map<std::set<int>, int, MyHash<std::set<int>>> triangleIndices;

    //std::vector<std::vector<int>> tetIncidentTriangles;
    std::vector<std::set<int>> tetIncidentTriangles;

    // given an edge Id, return it's plus/minus triangles (depending on the direction of travel)
    const std::vector<int>& getMinusTriangles(const int &edgeId, const bool &isDirectionLowerToUpper) const;
    const std::vector<int>& getPlusTriangles(const int &edgeId, const bool &isDirectionLowerToUpper) const;

    int singularEdgesNumber = 0;
    int pseudoSingularEdgesNumber = 0;
    int regularEdgesNumber = 0;

    void computeBoundingBoxes();
    void computeDomainBoundingBox();
    void computeRangeBoundingBox();

    void sortVertices();
    void printMesh();
    void perturbRangeValues(const double &perturbationEpsilon);


    // We only add upperLink/lowerLink to data, the rest of data is unchagned
    void computeUpperLowerLinkAndStar();
    void computeSingularEdgeTypes();

    // From the tet soup, get the edges and triangles and tet-incidence of the triangles
    void computeCombinatorialStructure();

    // Give the edge (aIndex, bIndex), is the vertex vIndex from its link in the upper and lower link of the edge
    // We assume that aIndex < bIndex for consistent orientation.
    bool isUpperLinkEdgeVertex(int, int, int);





    std::vector<std::set<int>> skeletonGraph;

    int getEdgeType(std::array<int, 2> edge)
    {
        if (edge[0] > edge[1])
        {
            std::swap(edge[0], edge[1]);
        }

        return this->edgeSingularTypes[edge];
    }


    const std::vector<int> findShortestPath(const std::vector<int> &source, const std::set<int> &sink) const;
    void markPseudoSingularEdges(const std::vector<int> &shortestPath);



    // Not currently used
    int computeSingularSetConnectivity();
    void singularTraversalBFS(const std::vector<int> &roots, std::vector<bool> &visited);
    int computeSingularSetConnectivity2();

};
