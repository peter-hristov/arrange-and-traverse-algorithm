#pragma once

#include "./TetMesh.h"

namespace timt
{
    void generateTopologyGraph(TetMesh &tetMesh, const std::array<double, 2> &p);

    void writeTopologyGraphToVTP(
            const std::vector<std::array<float, 3>> &domainCoordinates,
            const std::vector<Point_2> &ragneCoordinates,
            const std::vector<std::pair<int, int>> &edges,
            const std::vector<double> &vertexDistances,
            const std::string &filename);

};
