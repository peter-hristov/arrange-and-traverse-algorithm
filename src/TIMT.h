#pragma once

#include "./TetMesh.h"
#include "./ReebSpace2.h"
#include "./Arrangement.h"

namespace timt
{
    struct TopologyGraph
    {
        std::vector<double> vertexRangeDistances;
        std::vector<Point_2> vertexRangeCoordinates;
        std::vector<std::array<float, 3>> vertexDomainCoordinates;

        std::set<std::pair<int, int>> edges;
    };

    TopologyGraph computeExactTopologyGraph(TetMesh &tetMesh, Arrangement &singularArrangement, ReebSpace2 &reebSpace, const std::array<double, 2> &p);
    TopologyGraph computeInexactTopologyGraph(TetMesh &tetMesh, const std::array<double, 2> &p);

    void writeTopologyGraphToVTP(
            const TopologyGraph &tg,
            const std::string &filename);

};
