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



Fiber Fiber::computeLabeledFiber(TetMesh &tetMesh, Arrangement &singularArrangement, ReebSpace2 &reebSpace, std::array<double, 2> controlPoint, const std::set<int> &selectedSheetIds, const bool debugPrint)
{
    const seeds::SeedSet fiberSeeds = seeds::computeFiberSeedSet(tetMesh, singularArrangement, reebSpace, controlPoint, selectedSheetIds);

    if (debugPrint)
    {
        std::cerr << "Selected fibers :\n";
        for (const auto &[triangleId, componentId] : fiberSeeds)
        {
            const int sheetId = reebSpace.correspondenceGraphDS.find(componentId);
            std::cerr << "triangleId : " << triangleId << ", sheetId :  " << sheetId << std::endl;
        }
        std::cerr << std::endl;

        std::cerr << "Control point : (" << controlPoint[0] << ", " << controlPoint[1] << ")\n";
    }

    return Fiber::growFiberFromSeedSet(tetMesh, reebSpace, controlPoint, fiberSeeds);
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

    const auto seedTriangleBarycentricCoordinates = tetMesh.tryComputeActivePointCoordinates(triangleId, P);
    if (!seedTriangleBarycentricCoordinates) 
    {
        std::cerr << "Seed triangle does not contain control point.";
        return {};
    }

    fc.trianglePointCoordinates[triangleId] = *seedTriangleBarycentricCoordinates;

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
            if (const auto barycentricCoordinates = tetMesh.tryComputeActivePointCoordinates(nbTriangleId, P)) 
            {
                fc.trianglePointCoordinates[nbTriangleId] = *barycentricCoordinates;
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



bool FiberComponent::doesComponentContainTet(const TetMesh &tetMesh, const Arrangement &singularArrangement, const int seedTriangleId, const std::vector<int> &tetTriangleIds, const CartesianPoint &controlPoint)
{
    std::queue<int> bfsQueue;
    bfsQueue.push(seedTriangleId);

    std::vector<bool> visited(tetMesh.triangleIndices.size(), false);
    visited[seedTriangleId] = true;

    while (!bfsQueue.empty())
    {
        const int currentTriangleId = bfsQueue.front();
        bfsQueue.pop();

        // If it's the one we want, we are done
        if (std::find(tetTriangleIds.begin(), tetTriangleIds.end(), currentTriangleId) != tetTriangleIds.end())
        {
            return true;
        }

        for (const int &neighbourTriangleId : tetMesh.tetIncidentTriangles[currentTriangleId])
        {
            // Skip if it's visited
            if (visited[neighbourTriangleId]) { continue; }

            // Skip if it's not active
            if (false == tetMesh.isTriangleActive(neighbourTriangleId, controlPoint)) { continue; }

            bfsQueue.push(neighbourTriangleId);
            visited[neighbourTriangleId] = true;
        }
    }

    return false;
}

int Fiber::whichComponentContainsTet(const TetMesh &tetMesh, const Arrangement &singularArrangement, const std::vector<std::pair<int, int>> &fiberSeeds, const CartesianPoint &controlPoint, const int &tetId)
{
    // Unpack the triangles of the tet, if we intersect any one of those, we are done
    const int a = tetMesh.tetrahedra[tetId][0];
    const int b = tetMesh.tetrahedra[tetId][1];
    const int c = tetMesh.tetrahedra[tetId][2];
    const int d = tetMesh.tetrahedra[tetId][3];

    const std::vector<int> tetTriangleIds = {
        tetMesh.triangleIndices.at({a, b, c}),
        tetMesh.triangleIndices.at({a, b, d}),
        tetMesh.triangleIndices.at({a, c, d}),
        tetMesh.triangleIndices.at({b, c, d}),
    };


    for (const auto &[triangleId, componentId] : fiberSeeds)
    {
        if (FiberComponent::doesComponentContainTet(tetMesh, singularArrangement, triangleId, tetTriangleIds, controlPoint))
        {
            return componentId;
        }
    }

    return -1;
}
