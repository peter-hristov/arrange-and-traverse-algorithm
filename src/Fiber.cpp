#include "./CGALTypedefs.h"

#include "./Fiber.h"

#include "./io.h"
#include "./Timer.h"
#include "./SeedSet.h"
#include "./FiberSurface.h"

#include <cstdio>
#include <queue>
#include <unistd.h>
#include <unordered_map>
#include <unordered_set>



Fiber Fiber::computeLabeledFiber(TetMesh &tetMesh, Arrangement &singularArrangement, ReebSpace2 &reebSpace, std::array<double, 2> controlPoint, const std::set<int> &selectedSheetIds)
{
    const std::vector<std::pair<int, int>> fiberSeeds = SeedSet::computeFiberSeedSet(tetMesh, singularArrangement, reebSpace, controlPoint, selectedSheetIds);
    return Fiber::growFiberFromSeedSet(tetMesh, reebSpace, controlPoint, fiberSeeds);
}


std::optional<std::array<float, 3>> FiberComponent::tryComputePointCoordinates(const TetMesh& tetMesh, const int triangleId, const CartesianPoint& P)
{
    // Unpack the triangle vertices and put them it into points
    const std::set<int>& triangleVertices = tetMesh.triangles[triangleId];

    const auto itA = triangleVertices.begin();
    const CartesianPoint A(tetMesh.vertexCoordinatesF[*itA], tetMesh.vertexCoordinatesG[*itA]);

    const auto itB = std::next(itA);
    const CartesianPoint B(tetMesh.vertexCoordinatesF[*itB], tetMesh.vertexCoordinatesG[*itB]);

    const auto itC = std::next(itB);
    const CartesianPoint C(tetMesh.vertexCoordinatesF[*itC], tetMesh.vertexCoordinatesG[*itC]);

    CartesianPoint tri[3] = {A, B, C};
    if (CGAL::bounded_side_2(tri, tri + 3, P) != CGAL::ON_BOUNDED_SIDE) return std::nullopt;

    std::array<double, 3> bc;
    CGAL::Barycentric_coordinates::triangle_coordinates_2(A, B, C, P, bc.begin());

    // Fetch 3D domain coords
    const auto& A3 = tetMesh.vertexDomainCoordinates[*itA];
    const auto& B3 = tetMesh.vertexDomainCoordinates[*itB];
    const auto& C3 = tetMesh.vertexDomainCoordinates[*itC];

    // Interpolate
    std::array<float, 3> result{
        static_cast<float>(bc[0] * A3[0] + bc[1] * B3[0] + bc[2] * C3[0]),
            static_cast<float>(bc[0] * A3[1] + bc[1] * B3[1] + bc[2] * C3[1]),
            static_cast<float>(bc[0] * A3[2] + bc[1] * B3[2] + bc[2] * C3[2])
    };

return result;
}

FiberComponent FiberComponent::growFiberComponentFromSeedPair(const TetMesh& tetMesh, const std::array<double, 2>& controlPoint, const int triangleId, const int sheetId)
{
    const CartesianPoint P(controlPoint[0], controlPoint[1]);

    FiberComponent fc;
    fc.sheetId = sheetId;

    std::queue<int> bfsQueue;
    bfsQueue.push(triangleId);

    // Also our visited array
    std::unordered_map<int, int> parent;
    parent[triangleId] = triangleId;

    const auto seedTriangleBarycentricCoordinates = FiberComponent::tryComputePointCoordinates(tetMesh, triangleId, P);
    if (!seedTriangleBarycentricCoordinates) 
    {
        std::cerr << "Seed triangle does not contain control point.";
        return {};
    }

    fc.triangleBarycentricCoordinates[triangleId] = *seedTriangleBarycentricCoordinates;

    while (!bfsQueue.empty()) 
    {
        const int currentTriangleId = bfsQueue.front(); bfsQueue.pop();

        for (const int nbTriangleId : tetMesh.tetIncidentTriangles[currentTriangleId]) 
        {
            // If the neighbour has been visited (it has a parent)
            if (parent.contains(nbTriangleId)) 
            {
                // If this edge closes a loop (the neighbour is our parent), add an edge to the fiber component
                if (nbTriangleId != parent[currentTriangleId])
                {
                    fc.edges.emplace_back(currentTriangleId, nbTriangleId);
                }

                continue;
            }

            // If this triangle is active, compute it's barycentricCoordinates
            if (const auto barycentricCoordinates = FiberComponent::tryComputePointCoordinates(tetMesh, nbTriangleId, P)) 
            {
                fc.triangleBarycentricCoordinates[nbTriangleId] = *barycentricCoordinates;
                fc.edges.emplace_back(currentTriangleId, nbTriangleId);

                bfsQueue.push(nbTriangleId);
                parent[nbTriangleId] = currentTriangleId;
            }
        }
    }

    return fc;
}


Fiber Fiber::growFiberFromSeedSet(const TetMesh &tetMesh, ReebSpace2 &reebSpace, const std::array<double, 2> &controlPoint, const std::vector<std::pair<int, int>> &fiberSeeds)
{
    Fiber fiber;
    fiber.components.reserve(fiberSeeds.size());

    for (const auto &[triangleId, componentId] : fiberSeeds)
    {
        const int sheetId = reebSpace.correspondenceGraphDS.find(componentId);
        fiber.components.emplace_back(FiberComponent::growFiberComponentFromSeedPair(tetMesh, controlPoint, triangleId, sheetId));
    }

    return fiber;
}

