#include "./Augmentation.h"
#include "./LoadingBar.hpp"
#include "src/Arrangement.h"


int containsVertexInFiber(const TetMesh& tetMesh, const Fiber& fb, const int sheetId)
{
    for (const FiberComponent& fc : fb.components)
    {
        for (const auto& [triangleId, barycentricCoordinates] : fc.trianglePointCoordinates)
        {
            for (const int triangleVertexId : tetMesh.triangles[triangleId])
            {
                if (triangleVertexId == sheetId)
                {
                    return fc.sheetId;
                }
            }
        }
    }
    return -1;
};

std::map<int, std::vector<int>> augmentation::computeRegularVerticesSheets(TetMesh &tetMesh, Arrangement &singularArrangement, ReebSpace2 &reebSpace)
{
    // For a mini perturbation for the regular vertex points
    static std::mt19937 gen(std::random_device{}());
    std::uniform_real_distribution<double> dist(-0.0000001, 0.0000001);

    std::map<int, std::vector<int>> sheetRegularVertices;

    LoadingBar bar(40, "Augmenting Reeb space...");

    //Timer::start();
    for (int vertexId = 0 ; vertexId < tetMesh.isVertexSingular.size() ; vertexId++)
    {
        if (tetMesh.isVertexSingular[vertexId]) { continue; }

        const std::array<double, 2> controlPoint = {tetMesh.vertexCoordinatesF[vertexId] + dist(gen), tetMesh.vertexCoordinatesG[vertexId] + dist(gen)};

        //std::cout << "Computing labeled fiber for vertex id " << i << std::endl;
        Fiber fb = Fiber::computeLabeledFiber(tetMesh, singularArrangement, reebSpace, controlPoint, {});

        const int sheetId = containsVertexInFiber(tetMesh, fb, vertexId);

        if (sheetId == -1)
        {
            std::cerr << "Error in computation!\n";
        }
        else
        {
            sheetRegularVertices[sheetId].push_back(vertexId);
        }

        bar.update((100 * vertexId) / tetMesh.isVertexSingular.size());
    }

    return sheetRegularVertices;
}
