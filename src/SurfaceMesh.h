#pragma once

#include "./FiberPoint.h"

class SurfaceMesh
{
    public:
        std::vector<std::array<float, 3>> vertexCoordinates;
        std::vector<std::array<int, 3>> triangles;
        std::vector<double> edgeParam;


        std::vector<FiberPoint> getFiberPoints()
        {
            std::vector<FiberPoint> allFiberPoints;

            for (const auto &triangle : this->triangles)
            {
                for (const auto &vertex : triangle)
                {
                    allFiberPoints.push_back(FiberPoint(
                                vertexCoordinates[vertex],
                                {1.0, 1.0, 0.0},
                                1,
                                -1
                                ));
                }

            }

            return allFiberPoints;

        }
};
