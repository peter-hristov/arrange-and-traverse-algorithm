#include "./CGALTypedefs.h"

#include "./Fiber.h"
#include "./DisjointSet.h"
#include "./FiberPoint.h"
#include "./FiberGraph.h"

#include <queue>
#include <unistd.h>



std::array<double, 3> computeBarycentricCoordinates(const TetMesh &tetMesh, const int triangleId, const std::array<double, 2> &fiberPoint)
{
    return {0.33, 0.33, 0.33};

    // Unpack the actual numerical  coordinates of the vertices of the triangle
    std::vector<CartesianPoint> triangleVertexCoordinates;

    for (const int &vertexId : tetMesh.triangles[triangleId])
    {

        triangleVertexCoordinates.push_back(
                CartesianPoint(
                    tetMesh.vertexCoordinatesF[vertexId], 
                    tetMesh.vertexCoordinatesG[vertexId]
                    ));
    }

    std::array<double, 3> barycentricCoordinates;

    CartesianPoint P(fiberPoint[0], fiberPoint[1]);

    CGAL::Barycentric_coordinates::triangle_coordinates_2(
            triangleVertexCoordinates[0], 
            triangleVertexCoordinates[0], 
            triangleVertexCoordinates[0], 
            P, barycentricCoordinates.begin());

    return barycentricCoordinates;
}


std::vector<FiberPoint> BFSFiberSearch(const int &seedTriangleId, const TetMesh &tetMesh, const FiberGraph &pg, const std::array<double, 2> &fiberPoint, const int sheetId)
{
    std::vector<FiberPoint> faceFibers;
    const std::array<float, 3> sheetColour = fiber::fiberColours[sheetId % fiber::fiberColours.size()];

    // Set up the BFs
    std::queue<int> bfsQueue;
    bfsQueue.push(seedTriangleId);

    std::vector<int> parent(tetMesh.triangles.size(), -1);
    parent[seedTriangleId] = seedTriangleId;

    // Set up the array to hold all the discovered segments
    std::set<std::pair<int, int>> fiberSegments;

    while (false == bfsQueue.empty())
    {
        const int currentTriangleId = bfsQueue.front();
        bfsQueue.pop();

        for (const int &neighbourTriangleId : tetMesh.tetIncidentTriangles[currentTriangleId])
        {
            // If the neighbour is not active, skip it
            if (false == pg.componentRoot.contains(neighbourTriangleId))
            {
                continue;
            }


            // If the neighbour has NOT been visited, visited and add fiber segment 
            if (parent[neighbourTriangleId] == -1)
            {
                parent[neighbourTriangleId] = currentTriangleId;
                bfsQueue.push(neighbourTriangleId);

                fiberSegments.insert({
                        std::min(currentTriangleId, neighbourTriangleId),
                        std::max(currentTriangleId, neighbourTriangleId)
                        });

            }

            // If the neighbour has been visited and it's 
            if (neighbourTriangleId != parent[currentTriangleId])
            {
                fiberSegments.insert({
                        std::min(currentTriangleId, neighbourTriangleId),
                        std::max(currentTriangleId, neighbourTriangleId)
                        });

            }


        }
    }


    // Set up the fiber points for each segment
    for (const auto &[triangleIdA, triangleIdB] : fiberSegments)
    {
        const std::array<double, 3> triangleABarycentricCoordinates = computeBarycentricCoordinates(tetMesh, triangleIdA, fiberPoint);
        const std::array<double, 3> triangleBBarycentricCoordinates = computeBarycentricCoordinates(tetMesh, triangleIdB, fiberPoint);

        FiberPoint fb(
                triangleABarycentricCoordinates[0], 
                triangleABarycentricCoordinates[1], 
                tetMesh.getTriangleVerticesCoordinates(triangleIdA),
                sheetColour,
                sheetId,
                triangleIdA
                );

        FiberPoint fb2(
                triangleBBarycentricCoordinates[0], 
                triangleBBarycentricCoordinates[1], 
                tetMesh.getTriangleVerticesCoordinates(triangleIdB),
                sheetColour,
                sheetId,
                triangleIdB
                );

        faceFibers.push_back(fb);
        faceFibers.push_back(fb2);

    }




    return faceFibers;
}


std::vector<FiberPoint> processFiberGraph(const TetMesh &tetMesh, Arrangement &arrangement, ReebSpace2 &reebSpace, const std::array<double, 2> &fiberPoint, FiberGraph &pg, const std::set<int> activeSheets)
{
    // We first need to sort the triangles in the fiber graph
    std::vector<FiberPoint> faceFibers;

    for (const auto &[triangleId, componentId] : pg.componentRoot)
    {
        const int sheetId = reebSpace.correspondenceGraphDS.find(componentId);

        if (false == activeSheets.contains(sheetId))
        {
            continue;
        }

        std::vector<FiberPoint> newFaceFibers = BFSFiberSearch(triangleId, tetMesh, pg, fiberPoint, sheetId);

        faceFibers.insert(
                faceFibers.end(),
                newFaceFibers.begin(),
                newFaceFibers.end()
                );

    }

    return faceFibers;

}

std::vector<FiberPoint> fiber::computeFiberFromFiberGraph(const TetMesh &tetMesh, Arrangement &arrangement, ReebSpace2 &reebSpace, const std::array<double, 2> &fiberPoint)
{
    // The output of this function
    std::vector<FiberPoint> faceFibers;

    // Determine which face has been clicked on
    Face_const_handle activeFace = arrangement.getActiveFace(fiberPoint);
    const int activeFaceId = arrangement.arrangementFacesIdices[activeFace];

    // Determine which sheets are active (contained the active face)
    std::set<int> activeSheets;
    std::cout << "The active faces are : ";
    for (const int &componentId : reebSpace.correspondenceGraph[activeFaceId])
    {
        const int sheetId = reebSpace.correspondenceGraphDS.find(componentId);
        std::cout << sheetId << " ";

        activeSheets.insert(sheetId);
    }
    std::cout << std::endl;

    // Go through all faces and pull their representative fibers for their sheets
    for (auto face = arrangement.arr.faces_begin(); face != arrangement.arr.faces_end(); ++face) 
    {
        if (face->is_unbounded() || !face->has_outer_ccb())
        {
            continue;
        }

        const int faceId = face->data();

        FiberGraph &pg = reebSpace.fiberGraphsPerFace[faceId];

        for (const int &componentId : reebSpace.correspondenceGraph[faceId])
        {
            const int sheetId = reebSpace.correspondenceGraphDS.find(componentId);

            if (activeSheets.contains(sheetId))
            {
                std::vector<FiberPoint> newFaceFibers = processFiberGraph(tetMesh, arrangement, reebSpace, fiberPoint, pg, activeSheets);

                faceFibers.insert(
                        faceFibers.end(),
                        newFaceFibers.begin(),
                        newFaceFibers.end()
                        );
            }
        }

    }


    return faceFibers;
}

std::vector<FiberPoint> fiber::computeFiber(const TetMesh &tetMesh, Arrangement &arrangement, ReebSpace &reebSpace, const std::array<double, 2> &fiberPoint, const int reebSheetIdOnly = -1)
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

    if (reebSheetIdOnly == -1)
    {
        std::cout << "There are " << reebSpace.fiberSeeds[activeFaceId].size() << " fiber components with (sheet IDs, sorted IDs): ";
    }

    //vector<int> sheetIds;
    for (const auto &[fiberComponentId, triangleId] : reebSpace.fiberSeeds[activeFaceId])
    {
        const int sheetId = reebSpace.correspondenceGraph.findElement({activeFaceId, fiberComponentId});

        if (reebSheetIdOnly == -1 || sheetId == reebSheetIdOnly)
        {
            bfsQueue.push(triangleId);
            std::cout << "Adding seed for sheet ID " << sheetId << "\n";
        }

        //const int sheetId = reebSpace.findTriangle({activeFaceId, fiberComponentId});
        //triangleColour[triangleId] = sheetConsequitiveIndices[sheetId] % fiberColours.size();
        triangleSheetId[triangleId] = sheetId;

        //sheetIds.push_back(sheetId);

        if (reebSheetIdOnly == -1)
        {
            printf("(%d, %d) ", sheetId, reebSpace.sheetConsequitiveIndices[sheetId]);
        }
    }
    //std::cout << std::endl;

    CartesianPoint P(fiberPoint[0], fiberPoint[1]);
    std::vector<FiberPoint> faceFibers;

    while (false == bfsQueue.empty())
    {
        const int currentTriangleId = bfsQueue.front();
        const int currentSheeId = triangleSheetId[currentTriangleId];
        bfsQueue.pop();

        const std::array<float, 3> sheetColour = fiber::fiberColours[reebSpace.sheetConsequitiveIndices[currentSheeId] % fiber::fiberColours.size()];

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



//
// Brute force fiber computation
//

//void Data::computeTetExitPoints(const GLfloat u, const GLfloat v, const std::vector<float> color)
//{
    //this->faceFibers.clear();
    //this->tetsWithFibers = vector<bool>(this->tetMesh.tetrahedra.size(), false);

    ////
    //// Get the ID of the face we are intersecting
    ////

    //// The query point (u, v)
    //Point_2 query_point(u, v);

    //// Locate the point in the arrangement
    //CGAL::Object result = this->arrangement.pl->locate(query_point);

    //// Try to assign to a face, edge or a vertex
    //Arrangement_2::Face_const_handle face;
    //Arrangement_2::Halfedge_const_handle edge;
    //Arrangement_2::Vertex_const_handle vertex;

    //int currentFaceID = 0;

    //if (CGAL::assign(face, result)) 
    //{
        //currentFaceID = this->arrangement.arrangementFacesIdices[face];
    //} 
    //// If we are on an edge, just grad an adjacent face
    //else if (CGAL::assign(edge, result)) 
    //{
        //face = edge->face();
        //currentFaceID = this->arrangement.arrangementFacesIdices[face];
    //} 
    //// If we are on a vertex grab an indicent edge and get its face
    //else if (CGAL::assign(vertex, result)) 
    //{
        //edge = vertex->incident_halfedges();
        //face = edge->face();
        //currentFaceID = this->arrangement.arrangementFacesIdices[face];
    //} else 
    //{
        //assert(false);
    //}

    //// For every tet, compute the two exit points
    //for(size_t tetId = 0 ; tetId < this->tetMesh.tetrahedra.size(); tetId++)
    //{
        //const auto tet = this->tetMesh.tetrahedra[tetId];

        //// For every triangle in every tet, get the fiber point in it
        //for(int i = 0 ; i < 4 ; i++)
        //{
            //for(int j = i + 1 ; j < 4 ; j++)
            //{
                //for(int k = j + 1 ; k < 4 ; k++)
                //{
                    //float x1 = this->tetMesh.vertexCoordinatesF[tet[i]];
                    //float y1 = this->tetMesh.vertexCoordinatesG[tet[i]];

                    //float x2 = this->tetMesh.vertexCoordinatesF[tet[j]];
                    //float y2 = this->tetMesh.vertexCoordinatesG[tet[j]];

                    //float x3 = this->tetMesh.vertexCoordinatesF[tet[k]];
                    //float y3 = this->tetMesh.vertexCoordinatesG[tet[k]];


                    //const float xmin = std::min({x1, x2, x3});
                    //const float xmax = std::max({x1, x2, x3});
                    //const float ymin = std::min({y1, y2, y3});
                    //const float ymax = std::max({y1, y2, y3});

                    //// This triangle is not relevant, point is outside the bounding box
                    //if (u < xmin || u > xmax || v < ymin || v > ymax) 
                    //{
                        //continue;
                    //}


                    //float det = (x1 - x3) * (y2 - y3) - (x2 - x3) * (y1 - y3);

                    //float alpha = ((y2 - y3) * (u - x3) + (x3 - x2) * (v - y3)) / det;
                    //float beta = ((y3 - y1) * (u - x3) + (x1 - x3) * (v - y3)) / det;
                    //float gamma = 1 - alpha - beta;

                    //// Are we inside the triangle. We exclude the 0 and 1 because weird things happen there
                    //if (alpha > 0 && beta > 0 && gamma > 0 && alpha < 1 && beta < 1 && gamma < 1)
                    //{
                        ////printf("In triangle %ld, %ld, %ld in tet %ld.\n", tet[i], tet[j], tet[k], tetId);
                        ////printf("In triangle (%f, %f) | (%f, %f) | (%f, %f) comparing with point (%f, %f) and alpha = %f, betta = %f, gamma = %f.\n", x1, y1, x2, y2, x3, y3, x, y, alpha, betta, gamma);

                        ////const std::set<int> triangle({tet[i], tet[j], tet[k]});

                        //const int triangleVertexA = tet[i];
                        //const int triangleVertexB = tet[j];
                        //const int triangleVertexC = tet[k];

                        ////printf("Current triangle (%d, %d, %d)\n", triangleVertexA, triangleVertexB, triangleVertexC);


                        ////for (const auto [key, value] : this->faceDisjointSets[currentFaceID].data)
                        ////{
                        ////cout << "Triangle ";
                        ////for (const auto v : key)
                        ////{
                        ////cout << v << " ";
                        ////}

                        ////cout << " with root " << this->faceDisjointSets[currentFaceID].find(value) << " ( " << this->faceDisjointSets[currentFaceID].findTriangle(key) << ") " << endl;

                        ////}


                        ////for (const auto &[key, value] : this->reebSpace.data)
                        ////{
                        ////printf("Face ID = %d, fiber component root = %d, SheetID = %d\n", key.first, key.second, this->reebSpace.findTriangle(key));

                        ////}


                        //// Which sheets does this fiber belong to?
                        //// 1. Triangle -> Face ComponentID
                        //const int triangleID = this->tetMesh.triangleIndices[std::set<int>({triangleVertexA, triangleVertexB, triangleVertexC})];
                        //const int componentID = this->reebSpace.preimageGraphs[currentFaceID].findElement(triangleID);
                        ////printf("The face ID is %d and the component ID is = %d\n", currentFaceID, componentID);

                        //// 2. Fac ComponentID -> Reeb Space Sheet
                        ////const int pairToHIndex = this->vertexHtoIndex[{currentFaceID, componentID}];
                        ////const int sheetID = this->reebSpace.findTriangle(pairToHIndex);
                        //const int sheetID = this->reebSpace.correspondenceGraph.findElement({currentFaceID, componentID});

                        ////printf("The Sheet ID is = %d\n", sheetID);

                        //const int sheetColourID = this->reebSpace.sheetConsequitiveIndices[sheetID];

                        //// 3. Get the colou of the sheet
                        //const array<float, 3> sheetColour = fiber::fiberColours[sheetColourID];


                        //FiberPoint fb(alpha, beta, {
                                //this->tetMesh.vertexDomainCoordinates[tet[i]],
                                //this->tetMesh.vertexDomainCoordinates[tet[j]],
                                //this->tetMesh.vertexDomainCoordinates[tet[k]],
                                //},
                                //sheetColour);

                        //this->faceFibers.push_back(fb);
                        //this->tetsWithFibers[tetId] = true;
                    //}
                //}
            //}
        //}
    //}
//}
