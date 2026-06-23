#pragma once

#include "./src/TetMesh.h"

namespace timt
{
    void generateTopologyGraph(TetMesh &tetMes, const std::array<double, 2> &ph);

    void writeTopologyGraphToVTP(
            const std::vector<std::array<float, 3>> &domainCoordinates,
            const std::vector<std::pair<int, int>> &edges,
            const std::vector<double> &vertexDistances,
            const std::string &filename);
}
