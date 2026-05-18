#include "./Augmentation.h"
#include "./LoadingBar.hpp"
#include "./io.h"


int containsVertexInFiber(const TetMesh& tetMesh, const Fiber& fb, const int vertexId)
{
    for (const FiberComponent& fc : fb.components)
    {
        for (const auto& [triangleId, barycentricCoordinates] : fc.trianglePointCoordinates)
        {
            for (const int triangleVertexId : tetMesh.triangles[triangleId])
            {
                if (triangleVertexId == vertexId)
                {
                    return fc.sheetId;
                }
            }
        }
    }
    return -1;
};

std::vector<int> augmentation::computeRegularVertexSheets(TetMesh &tetMesh, Arrangement &singularArrangement, ReebSpace2 &reebSpace)
{
    // For a mini perturbation for the regular vertex points
    static std::mt19937 gen(std::random_device{}());
    std::uniform_real_distribution<double> dist(-0.00000001, 0.00000001);

    std::vector<int> regularVertexSheet(tetMesh.isVertexSingular.size(), -1);

    LoadingBar bar(40, "Augmenting Reeb space...");

    for (int vertexId = 0 ; vertexId < tetMesh.isVertexSingular.size() ; vertexId++)
    {
        if (tetMesh.isVertexSingular[vertexId]) { continue; }

        const std::array<double, 2> controlPoint = {tetMesh.vertexCoordinatesF[vertexId] + dist(gen), tetMesh.vertexCoordinatesG[vertexId] + dist(gen)};

        Fiber fb = Fiber::computeLabeledFiber(tetMesh, singularArrangement, reebSpace, controlPoint, {});

        const int sheetId = containsVertexInFiber(tetMesh, fb, vertexId);

        if (sheetId == -1)
        {
            std::cerr << "Error in computation!\n";
        }
        else
        {
            regularVertexSheet[vertexId] = sheetId;
        }

        bar.update((100 * (vertexId + 1)) / tetMesh.isVertexSingular.size());
    }

    io::saveWithSheetData(tetMesh.originalMesh, regularVertexSheet, "regularVertexSheets.vtu");

    return regularVertexSheet;
}
