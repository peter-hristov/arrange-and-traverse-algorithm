#include "./CGALTypedefs.h"

#include "./Fiber.h"

#include "./io.h"
#include "./Timer.h"

#include "./FiberPoint.h"
#include "./SurfaceMesh.h"
#include "./FiberLabeling.h"

#include <queue>
#include <unistd.h>
#include <unordered_map>
#include <unordered_set>


SurfaceMesh fiber::computeSegmentedFiberSurface(TetMesh &tetMesh, Arrangement &singularArrangement, ReebSpace2 &reebSpace, const std::vector<std::array<double, 2>> &controlPoints, const std::set<int> &selectedSheets)
{
    const Point_2 startPoint(controlPoints[0][0], controlPoints[0][1]);
    const Point_2 endPoint(controlPoints[1][0], controlPoints[1][1]);
    const Segment_2 controlSegment(startPoint, endPoint);

    //Timer::start();
    const std::vector<std::tuple<K::FT, int, int>> intersectedSegments = singularArrangement.getIntersectedSegments2(tetMesh, controlSegment, true);
    //Timer::stop("Computed Alpha intersections           :");

    //Timer::start();
    SurfaceMesh surfaceMesh = io::computeFiberSurface(tetMesh.originalMesh, controlPoints[0][0], controlPoints[0][1], controlPoints[1][0], controlPoints[1][1]);
    //Timer::stop("Computing fiber surfaces with TTK      :");

    //Timer::start();
    std::vector<double> intersectionAlpha;
    intersectionAlpha.reserve(intersectedSegments.size());

    for (const auto &[alpha, edgeId, edgeType] : intersectedSegments)
    {
        //if (edgeType == 2 || edgeType == 0)

        if (edgeType == 2)
        {
            intersectionAlpha.emplace_back(CGAL::to_double(alpha));
        }
    }

    surfaceMesh.remesh(intersectionAlpha);
    //Timer::stop("Subdivided mesh                        :");



    //Timer::start();
    surfaceMesh.labelFiberSurface(tetMesh, singularArrangement, reebSpace, intersectedSegments, controlSegment);
    //Timer::stop("Computing triangle sheets 2            :");

    //Timer::start();
    surfaceMesh.filterTriangles(selectedSheets);
    //Timer::stop("Filtering out triangles                :");

    //surfaceMesh.printSheetHistogram(reebSpace);

    return surfaceMesh;
}

std::vector<FiberPoint> fiber::computeLabeledFiber(TetMesh &tetMesh, Arrangement &singularArrangement, ReebSpace2 &reebSpace, std::array<double, 2> controlPoint, const std::set<int> &selectedSheetIds)
{
    const std::vector<std::pair<int, int>> fiberSeeds = fiber::labeling::computeSeedFibers(tetMesh, singularArrangement, reebSpace, controlPoint, selectedSheetIds);
    return growSeedSet(tetMesh, singularArrangement, reebSpace, controlPoint, fiberSeeds);
}

std::vector<FiberPoint> fiber::growSeedSet(const TetMesh &tetMesh, Arrangement &arrangement, ReebSpace2 &reebSpace, const std::array<double, 2> &fiberPoint, const std::vector<std::pair<int, int>> &fiberSeeds)
{
    Face_const_handle activeFace = arrangement.getActiveFace(fiberPoint);
    const int activeFaceId = arrangement.arrangementFacesIdices[activeFace];

    //Timer::start();

    // The sizes of these data structures are linear in the size of the fiber, not an issue
    std::queue<int> bfsQueue;
    // This also acts as the visited array
    std::unordered_map<int, int> triangleSheetId;
    // Needed to close loops closed fibers, hack because we are not visiting tets, but triangles
    std::unordered_set<std::pair<int, int>, MyHash<std::pair<int, int>>> activeAdjacentTrianglesConnected;
    // Cache barycentric coordintes, they are expensive to compute
    std::unordered_map<int, std::array<double, 3>> triangleBarycentricCoordinates;

    //vector<int> sheetIds;
    for (const auto &[triangleId, componentId] : fiberSeeds)
    {
        const int sheetId = reebSpace.correspondenceGraphDS.find(componentId);
        bfsQueue.push(triangleId);
        triangleSheetId[triangleId] = sheetId;
    }
    //std::cout << std::endl;

    CartesianPoint P(fiberPoint[0], fiberPoint[1]);
    std::vector<FiberPoint> faceFibers;

    while (false == bfsQueue.empty())
    {
        const int currentTriangleId = bfsQueue.front();
        const int currentSheeId = triangleSheetId[currentTriangleId];
        bfsQueue.pop();

        const int sheetSortId = reebSpace.sheetOrder.at(currentSheeId);
        const std::array<float, 3> sheetColour = fiber::fiberColours[sheetSortId % fiber::fiberColours.size()];

        const std::set<int> triangleUnpacked = tetMesh.triangles[currentTriangleId];
        const std::vector<int> triangleIndices = std::vector<int>(triangleUnpacked.begin(), triangleUnpacked.end());

        std::array<double, 3> barycentricCoordinatesCurrent;

        if (triangleBarycentricCoordinates.contains(currentTriangleId))
        {
            barycentricCoordinatesCurrent = triangleBarycentricCoordinates[currentTriangleId];
        }
        else
        {
            const CartesianPoint A(tetMesh.vertexCoordinatesF[triangleIndices[0]], tetMesh.vertexCoordinatesG[triangleIndices[0]]);
            const CartesianPoint B(tetMesh.vertexCoordinatesF[triangleIndices[1]], tetMesh.vertexCoordinatesG[triangleIndices[1]]);
            const CartesianPoint C(tetMesh.vertexCoordinatesF[triangleIndices[2]], tetMesh.vertexCoordinatesG[triangleIndices[2]]);
            CGAL::Barycentric_coordinates::triangle_coordinates_2(A, B, C, P, barycentricCoordinatesCurrent.begin());
            triangleBarycentricCoordinates[currentTriangleId] = barycentricCoordinatesCurrent;
        }

        // Sanity check
        assert(barycentricCoordinatesCurrent[0] > 0 && barycentricCoordinatesCurrent[1] > 0 && barycentricCoordinatesCurrent[2] > 0);

        // Look at the neighbours
        for (const int &neighbourTriagleId : tetMesh.tetIncidentTriangles[currentTriangleId])
        {
            if (neighbourTriagleId == currentTriangleId)
            {
                continue;
            }

            // We can't skip visited neighbours, because there may be a fiber between us (completing a circle)
            //if (triangleColour.contains(neighbourTriagle)) { continue; }

            const std::set<int> triangle2Unpacked = tetMesh.triangles[neighbourTriagleId];
            const std::vector<int> triangle2Indices = std::vector<int>(triangle2Unpacked.begin(), triangle2Unpacked.end());



            // The neighbour is active if we've already seen it
            bool isActive = triangleSheetId.contains(neighbourTriagleId);

            // Or if the image of the triangle contains the fiber points
            // We use a fast test to avoid having to use barycentric coordinates all the time
            if (false == isActive)
            {
                CartesianPoint A(tetMesh.vertexCoordinatesF[triangle2Indices[0]], tetMesh.vertexCoordinatesG[triangle2Indices[0]]);
                CartesianPoint B(tetMesh.vertexCoordinatesF[triangle2Indices[1]], tetMesh.vertexCoordinatesG[triangle2Indices[1]]);
                CartesianPoint C(tetMesh.vertexCoordinatesF[triangle2Indices[2]], tetMesh.vertexCoordinatesG[triangle2Indices[2]]);

                std::vector<CartesianPoint> triangle = {A, B, C};
                const auto result = CGAL::bounded_side_2(triangle.begin(), triangle.end(), P);

                isActive = (result == CGAL::ON_BOUNDED_SIDE);
            }

            // Determine if the triangle is active
            if (isActive)
            {
                // Only add the neighbour if we have not already visited it
                if (false == triangleSheetId.contains(neighbourTriagleId))
                {
                    // BFS things
                    bfsQueue.push(neighbourTriagleId);
                    triangleSheetId[neighbourTriagleId] = currentSheeId;
                }

                // Even if we have aleady added a neighbour, maybe there still isn't a fiber between us (for finishing loops)

                // At this point, we know that both us and we neighbour are active, is there already a fiber between us? Then skip
                if (activeAdjacentTrianglesConnected.contains({currentTriangleId, neighbourTriagleId}))
                {
                    continue;
                }
                else
                {
                    activeAdjacentTrianglesConnected.insert({currentTriangleId, neighbourTriagleId});
                    activeAdjacentTrianglesConnected.insert({neighbourTriagleId, currentTriangleId});
                }




                // Compute barycentric coordinates for drawing
                std::array<double, 3> barycentricCoordinatesNeighbour;
                if (triangleBarycentricCoordinates.contains(neighbourTriagleId))
                {
                    barycentricCoordinatesNeighbour = triangleBarycentricCoordinates[neighbourTriagleId];
                }
                else
                {
                    CartesianPoint A(tetMesh.vertexCoordinatesF[triangle2Indices[0]], tetMesh.vertexCoordinatesG[triangle2Indices[0]]);
                    CartesianPoint B(tetMesh.vertexCoordinatesF[triangle2Indices[1]], tetMesh.vertexCoordinatesG[triangle2Indices[1]]);
                    CartesianPoint C(tetMesh.vertexCoordinatesF[triangle2Indices[2]], tetMesh.vertexCoordinatesG[triangle2Indices[2]]);

                    CGAL::Barycentric_coordinates::triangle_coordinates_2(A, B, C, P, barycentricCoordinatesNeighbour.begin());
                    triangleBarycentricCoordinates[neighbourTriagleId] = barycentricCoordinatesNeighbour;
                }

                


                //
                // Add a fiber segment
                //
                FiberPoint fb(
                        barycentricCoordinatesCurrent[0], 
                        barycentricCoordinatesCurrent[1], 
                        {
                            tetMesh.vertexDomainCoordinates[triangleIndices[0]],
                            tetMesh.vertexDomainCoordinates[triangleIndices[1]],
                            tetMesh.vertexDomainCoordinates[triangleIndices[2]],
                        },
                        sheetColour);
                fb.sheetId = currentSheeId;
                fb.triangleId = currentTriangleId;
                faceFibers.push_back(fb);

                FiberPoint fb2(barycentricCoordinatesNeighbour[0], barycentricCoordinatesNeighbour[1], {
                        tetMesh.vertexDomainCoordinates[triangle2Indices[0]],
                        tetMesh.vertexDomainCoordinates[triangle2Indices[1]],
                        tetMesh.vertexDomainCoordinates[triangle2Indices[2]],
                        },
                        sheetColour);
                fb2.sheetId = currentSheeId;
                fb2.triangleId = neighbourTriagleId;
                faceFibers.push_back(fb2);

                //printf("Adding fiber between %d -> %d\n", currentTriangleId, neighbourTriagleId);
            }
        }
    }

    return faceFibers;

    //Timer::stop("Computed fiber in                      :");
}
