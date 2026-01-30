#pragma once

#include <array>
#include <vector>


// The point of intersection of a fiber and a triangle in the input tet mesh
class FiberPoint{

    public:

    // (x, y, z) domain coodrinates
    std::array<float, 3> point;

    // RGB colour value
    std::array<float, 3> colour;

    int sheetId;
    int triangleId;

    // The barycentric coordinates as well as the three triangle vertices
    FiberPoint(const float alpha, const float beta, const std::array<std::array<float, 3>, 3> &triangleVertices, const std::array<float, 3> _colour, const int _sheetId = -1, const int _triangleId = -1)
    {
        this->colour = _colour;
        this->sheetId = _sheetId;
        this->triangleId = _triangleId;

        point[0] = 
            alpha * triangleVertices[0][0] +
            beta * triangleVertices[1][0] +
            (1 - alpha - beta) * triangleVertices[2][0];

        point[1] = 
            alpha * triangleVertices[0][1] +
            beta * triangleVertices[1][1] +
            (1 - alpha - beta) * triangleVertices[2][1];

        point[2] = 
            alpha * triangleVertices[0][2] +
            beta * triangleVertices[1][2] +
            (1 - alpha - beta) * triangleVertices[2][2];
    }
};
