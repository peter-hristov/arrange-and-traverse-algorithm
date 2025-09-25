#include <CGAL/enum.h>
#include <cassert>
#include <utility>
#include <omp.h>

#include "./io.h"
#include "./Timer.h"
#include "./ReebSpace2.h"
#include "./DisjointSet.h"
#include "./PreimageGraph.h"
//#include "src/CGALTypedefs.h"

PolygonE_2 face_to_polygon(Arrangement_2::Face_const_handle f)
{
    PolygonE_2 poly;

    if (!f->is_unbounded())
    {
        // A face can have multiple outer CCBs (connected components of the boundary)
        auto ccb = f->outer_ccbs_begin();
        Arrangement_2::Ccb_halfedge_const_circulator curr = *ccb;

        do {
            poly.push_back(curr->source()->point());
            ++curr;
        } while (curr != *ccb);
    }

    return poly;
}


void ReebSpace2::bfs(TetMesh &tetMesh, Arrangement &singularArrangement, Face_const_handle startingFace, int sheetId, std::function<void(Face_const_handle)> callback)
{
    std::queue<Face_const_handle> traversalQueue;
    std::vector<bool> visited(singularArrangement.arr.number_of_faces(), false);

    traversalQueue.push(startingFace);
    visited[startingFace->data()] = true;

    while (false == traversalQueue.empty())
    {
        Face_const_handle currentFace = traversalQueue.front();
        traversalQueue.pop();

        auto circ = currentFace->outer_ccb();
        const auto start = circ;
        do
        {
            Face_const_handle twinFace = circ->twin()->face();
            const int &twinFaceId = twinFace->data();

            const bool doesTwinFaceHaveSheet = std::ranges::find(this->correspondenceGraph[twinFaceId], sheetId) != this->correspondenceGraph[twinFaceId].end();

            if (false == visited[twinFaceId] && doesTwinFaceHaveSheet)
            {
                callback(twinFace);
                traversalQueue.push(twinFace);
                visited[twinFaceId] = true;
            }

            ++circ;
        } while (circ != start);
    }
}

void ReebSpace2::postprocessSheets2(TetMesh &tetMesh, Arrangement &singularArrangement)
{

    std::vector<PolygonE_2> facePolygons(singularArrangement.arr.number_of_faces());

    for (auto face = singularArrangement.arr.faces_begin(); face != singularArrangement.arr.faces_end(); ++face)
    {
        facePolygons[face->data()] = face_to_polygon(face);
    }

    //std::vector<PolygonE_with_holes_2> polygonsPerSheet(PreimageGraph::componentCount);

    std::vector<PolygonE_with_holes_2> sheetPolygon(PreimageGraph::componentCount);
    std::vector<bool> processedSheets(PreimageGraph::componentCount, false);


    for (auto face = singularArrangement.arr.faces_begin(); face != singularArrangement.arr.faces_end(); ++face)
    {
        const int &faceId = face->data();
        const PolygonE_2 &facePolygon = facePolygons[faceId];

        const std::vector<int> &faceSheetIds = this->correspondenceGraph[faceId];

        for (const int &sheetId : faceSheetIds)
        {
            //std::cout << "The sheet Id is " << sheetId << " out of " << PreimageGraph::componentCount << std::endl;
            if (processedSheets[sheetId] == false)
            {
                bfs(tetMesh, singularArrangement, face, sheetId,
                        [&sheetPolygon, &facePolygon, &sheetId](Face_const_handle f) {
                        //CGAL::join(facePolygon, facePolygon, sheetPolygon[sheetId]);
                        });
            }
            else
            {
                processedSheets[sheetId] = true;
            }
        }
    }

    return;

    for (int i = 0 ; i < sheetPolygon.size() ; i++)
    {
        const PolygonE_with_holes_2 &pwh = sheetPolygon[i];

        // Outer boundary
        const PolygonE_2& outer = pwh.outer_boundary();

        std::cout << "Printing polygon " << i << std::endl;

        // Iterate over outer vertices
        for (auto vit = outer.vertices_begin(); vit != outer.vertices_end(); ++vit) 
        {
            std::cout << vit->x() << " " << vit->y() << "\n";
        }

        // Holes
        for (auto hit = pwh.holes_begin(); hit != pwh.holes_end(); ++hit) 
        {
            const PolygonE_2& hole = *hit;
            std::cout << "Hole:\n";
            for (auto vit = hole.vertices_begin(); vit != hole.vertices_end(); ++vit) 
            {
                std::cout << vit->x() << " " << vit->y() << "\n";
            }
        }

    }




}

void ReebSpace2::postprocessSheets(TetMesh &tetMesh, Arrangement &singularArrangement)
{
    // 1. Find the boundary faces of each sheet
    std::vector<std::set<int>> sheetBoundaryFaces(PreimageGraph::componentCount);

    std::vector<std::set<int>> correspondenceGraphSet(singularArrangement.arr.number_of_faces());
    for (auto face = singularArrangement.arr.faces_begin(); face != singularArrangement.arr.faces_end(); ++face)
    {
        const int &faceId = face->data();
        correspondenceGraphSet[face->data()] = std::set(this->correspondenceGraph[faceId].begin(), this->correspondenceGraph[faceId].end());
    }

    for (auto face = singularArrangement.arr.faces_begin(); face != singularArrangement.arr.faces_end(); ++face)
    {
        if (face->is_unbounded())
        {
            continue;
        }

        const int &faceId = face->data();
        const std::set<int> &faceSheets = correspondenceGraphSet[faceId];

        std::set<int> boundarySheets;

        auto circ = face->outer_ccb();
        const auto start = circ;
        do
        {
            const int &twinFaceId = circ->twin()->face()->data();
            const std::set<int> &twinFaceSheets = correspondenceGraphSet[twinFaceId];

            // If a sheet is in the face, but not in the twin, then it must be on the boundary
            // Otherwise all neighbouring faces will also have that sheet
            std::set_difference(
                    faceSheets.begin(), faceSheets.end(),
                    twinFaceSheets.begin(), twinFaceSheets.end(),
                    std::inserter(boundarySheets, boundarySheets.begin()));

            ++circ;
        } while (circ != start);

        for (const int &sheetId : boundarySheets)
        {
            sheetBoundaryFaces[sheetId].insert(faceId);
        }
    }


    // Count all faces per sheet
    std::vector<std::set<int>> sheetFaces(PreimageGraph::componentCount);
    for (auto face = singularArrangement.arr.faces_begin(); face != singularArrangement.arr.faces_end(); ++face)
    {
        if (face->is_unbounded())
        {
            continue;
        }

        const int &faceId = face->data();
        const std::set<int> &faceSheets = correspondenceGraphSet[faceId];

        for (const int &sheetId : faceSheets)
        {
            sheetFaces[sheetId].insert(faceId);
        }
    }

    //for (int i = 0 ; i < PreimageGraph::componentCount ; i++)
    //{
        //printf("Sheet with ID %d has %ld faces and %ld boundary faces.\n", i, sheetFaces[i].size(), sheetBoundaryFaces[i].size());
    //}
}
 
void ReebSpace2::loopFace(TetMesh &tetMesh, const Halfedge_const_handle &initialHalfEdge, std::queue<Halfedge_const_handle> &traversalQueue, std::vector<bool> &visited, const bool cachePreimageGraphs)
{
    PreimageGraph &pg = this->preimageGraphs[initialHalfEdge->face()->data()];
    if (cachePreimageGraphs) { preimageGraphsCached[initialHalfEdge].first = pg; }

    pg.updateComponentsRegular(tetMesh, this->edgeRegionSegments[initialHalfEdge->data()]);
    if (cachePreimageGraphs) { preimageGraphsCached[initialHalfEdge].second = pg; }

    Halfedge_const_handle currentHalfEdge = initialHalfEdge->next();
    do
    {
        pg.updateComponentsRegular(tetMesh, this->vertexRegionSegments[currentHalfEdge->prev()->data()]);
        if (cachePreimageGraphs) { preimageGraphsCached[currentHalfEdge].first = pg; }

        pg.updateComponentsRegular(tetMesh, this->edgeRegionSegments[currentHalfEdge->data()]);
        if (cachePreimageGraphs) { preimageGraphsCached[currentHalfEdge].second = pg; }

        if (false == visited[currentHalfEdge->twin()->face()->data()])
        {
            PreimageGraph &pg2 = this->preimageGraphs[currentHalfEdge->twin()->face()->data()];
            pg2 = pg;

            if (this->isHalfEdgePseudoSingular[currentHalfEdge->data()])
            {
                pg2.updateComponentsRegular(tetMesh, {this->edgeCrossingSegments[currentHalfEdge->data()]});
            }
            else
            {
                pg2.updateComponentsSingular(tetMesh, this->edgeCrossingSegments[currentHalfEdge->data()]);
            }

            traversalQueue.push(currentHalfEdge->twin());
            visited[currentHalfEdge->twin()->face()->data()] = true;
        }

        currentHalfEdge = currentHalfEdge->next();

    } while (currentHalfEdge != initialHalfEdge);

    pg.clear();
}


void ReebSpace2::traverse(TetMesh &tetMesh, Arrangement &singularArrangement, const bool cachePreimageGraphs)
{
    this->correspondenceGraph.resize(singularArrangement.arr.number_of_faces());
    this->preimageGraphs.resize(singularArrangement.arr.number_of_faces());

    // Find the outside face
    Halfedge_const_handle startingHalfedge;
    for (Face_const_iterator fit = singularArrangement.arr.faces_begin(); fit != singularArrangement.arr.faces_end(); ++fit) 
    {
        if (fit->is_unbounded()) 
        {
            startingHalfedge = *fit->holes_begin();
            break;
        }
    }

    std::queue<Halfedge_const_handle> traversalQueue;
    std::vector<bool> visited(singularArrangement.arr.number_of_faces(), false);

    traversalQueue.push(startingHalfedge->twin());
    visited[startingHalfedge->twin()->face()->data()] = true;

    // Seed the first face
    this->preimageGraphs[startingHalfedge->twin()->face()->data()].updateComponentsSingular(tetMesh, this->edgeCrossingSegments[startingHalfedge->data()]);

    // Make sure the outer face is visited as well, no need to go back
    visited[startingHalfedge->face()->data()] = true;

    while (false == traversalQueue.empty())
    {
        Halfedge_const_handle currentHalfEdge = traversalQueue.front();
        traversalQueue.pop();

        correspondenceGraph[currentHalfEdge->face()->data()] = this->preimageGraphs[currentHalfEdge->face()->data()].getUniqueComponents();
        loopFace(tetMesh, currentHalfEdge, traversalQueue, visited,  cachePreimageGraphs);
    }
}

bool ReebSpace2::doSegmentEndpointsOverlap(const Segment_2 &s1, const Segment_2 &s2)
{
    return 
        s1.source() == s2.source() || s1.source() == s2.target() ||
        s1.target() == s2.source() || s1.target() == s2.target();
}

void ReebSpace2::computeEdgeRegionSegments(const TetMesh &tetMesh, Arrangement &singularArrangement)
{
    // Set up singular half-edge segments
    //
    std::vector<Segment_2> singularSegments;
    singularSegments.reserve(singularArrangement.arr.number_of_edges());

    std::vector<int> singularSegmentsHalfEdgeIds;
    singularSegmentsHalfEdgeIds.reserve(singularArrangement.arr.number_of_edges());

    std::vector<Bbox> singularBoxes;
    singularBoxes.reserve(singularArrangement.arr.number_of_edges());

    for (auto he = singularArrangement.arr.halfedges_begin(); he != singularArrangement.arr.halfedges_end(); ++he)
    {
        // @TODO Not entirely sure this is super stable
        if (he < he->twin())
        {
            singularSegments.emplace_back(he->source()->point(), he->target()->point());
            singularBoxes.emplace_back(singularSegments.back().bbox());
            singularSegmentsHalfEdgeIds.emplace_back(he->data());
        }
    }

    // Set up regular segments
    //
    std::vector<Segment_2> regularSegments;
    regularSegments.reserve(tetMesh.regularEdgesNumber);

    std::vector<int> regularSegmentsIds;
    regularSegmentsIds.reserve(tetMesh.regularEdgesNumber);

    std::vector<Bbox> regularBoxes; 
    regularBoxes.reserve(tetMesh.regularEdgesNumber);

    for (const auto &[edge, type] : tetMesh.edgeSingularTypes) 
    {
        if (type == 1)
        {
            const std::array<int, 2> &edgeConst = edge;

            regularSegments.emplace_back(singularArrangement.arrangementPoints[edge[0]], singularArrangement.arrangementPoints[edge[1]]);
            regularBoxes.emplace_back(regularSegments.back().bbox());
            regularSegmentsIds.emplace_back(tetMesh.edgeIndices.at(edgeConst));
        }
    }

    this->edgeRegionSegments.resize(singularArrangement.arr.number_of_halfedges());

    #pragma omp parallel for schedule(dynamic)
    for (int i = 0 ; i < singularBoxes.size() ; i++)
    {
        std::map<K::FT, int> segmentRegionsOrdered;

        const Segment_2 &singularSegment = singularSegments[i];
        const K::Vector_2 v = singularSegment.target() - singularSegment.source();

        for (int j = 0 ; j < regularBoxes.size() ; j++)
        {
            // Quick reject using bounding boxes
            if (!CGAL::do_overlap(singularBoxes[i], regularBoxes[j]))
            {
                continue;
            }

            // We are only interested in interior intersections
            if (doSegmentEndpointsOverlap(singularSegments[i], regularSegments[j])) 
            { 
                continue; 
            }

            const Segment_2 &regularSegment = regularSegments[j];

            std::optional<std::variant<Point_2, Segment_2>> result = CGAL::intersection(singularSegment, regularSegment);

            if (!result)
            {
                continue;
            }

            if (std::holds_alternative<Point_2>(*result))
            {
                const Point_2 &p = std::get<Point_2>(*result);
                const K::Vector_2 wA = p - singularSegment.source();
                const K::FT tA = (wA * v) / v.squared_length();

                // Make sure the regular intersections are well ordered
                assert(false == segmentRegionsOrdered.contains(tA));

                segmentRegionsOrdered[tA] = regularSegmentsIds[j];
            }
        }

        this->edgeRegionSegments[singularSegmentsHalfEdgeIds[i]].reserve(segmentRegionsOrdered.size());

        for (const auto &[lambda, regularSegmentId] : segmentRegionsOrdered)
        {
            this->edgeRegionSegments[singularSegmentsHalfEdgeIds[i]].emplace_back(regularSegmentId, true);
        }
    }
}

void ReebSpace2::computeEdgeRegionSegments2(const TetMesh &tetMesh, Arrangement &singularArrangement)
{
    //
    // Compute the intersection point using the singular segments (not half-edges)
    // Faster, but needs postprocessing.
    //
    //std::vector<Segment_2> regularSegments;
    //regularSegments.reserve(tetMesh.regularEdgesNumber);

    //std::vector<Segment_2> singularSegments;
    //singularSegments.reserve(tetMesh.singularEdgesNumber);

    //std::vector<Box> regularBoxes; 
    //regularBoxes.reserve(tetMesh.regularEdgesNumber);

    //std::vector<Box> singularBoxes;
    //singularBoxes.reserve(tetMesh.singularEdgesNumber);

    //for (const auto &[edge, type] : tetMesh.edgeSingularTypes) 
    //{
        //if (type == 1)
        //{
            //regularSegments.emplace_back(arrangement.arrangementPoints[edge[0]], arrangement.arrangementPoints[edge[1]]);
            //regularBoxes.emplace_back(regularSegments.back().bbox(), &regularSegments.back());
        //}
        //else
        //{
            //singularSegments.emplace_back(arrangement.arrangementPoints[edge[0]], arrangement.arrangementPoints[edge[1]]);
            //singularBoxes.emplace_back(singularSegments.back().bbox(), &singularSegments.back());
        //}
    //}





    // Set up singular half-edge segments
    //
    std::vector<MySegment_2> singularSegments;
    singularSegments.reserve(singularArrangement.arr.number_of_edges());

    std::vector<Box> singularBoxes;
    singularBoxes.reserve(singularArrangement.arr.number_of_edges());

    for (auto he = singularArrangement.arr.halfedges_begin(); he != singularArrangement.arr.halfedges_end(); ++he)
    {
        // @TODO Not entirely sure this is super stable
        if (he < he->twin())
        {
            singularSegments.emplace_back(Segment_2(he->source()->point(), he->target()->point()), he);
            singularBoxes.emplace_back(singularSegments.back().seg.bbox(), &singularSegments.back());
        }
    }

    // Set up regular segments
    //
    std::vector<MySegment_2> regularSegments;
    regularSegments.reserve(tetMesh.regularEdgesNumber);

    std::vector<Box> regularBoxes; 
    regularBoxes.reserve(tetMesh.regularEdgesNumber);

    for (const auto &[edge, type] : tetMesh.edgeSingularTypes) 
    {
        if (type == 1)
        {
            const std::array<int, 2> &edgeConst = edge;

            regularSegments.emplace_back(
                    Segment_2(singularArrangement.arrangementPoints[edge[0]], singularArrangement.arrangementPoints[edge[1]]), 
                    std::nullopt, 
                    tetMesh.edgeIndices.at(edgeConst)
                    );
            regularBoxes.emplace_back(regularSegments.back().seg.bbox(), &regularSegments.back());
        }
    }

    // Callback that handles bounding box intersections
    //
    auto cb = [&](const Box& regularSegmentBox, const Box& singularSegmentBox) {
        
        const auto &regularSegmentHandle = *regularSegmentBox.handle();
        const auto &singularSegmentHandle  = *singularSegmentBox.handle();

        if (doSegmentEndpointsOverlap(regularSegmentHandle.seg, singularSegmentHandle.seg)) { return; }


        std::optional<std::variant<Point_2, Segment_2>> result = CGAL::intersection(regularSegmentHandle.seg, singularSegmentHandle.seg);

        if (result && std::holds_alternative<Point_2>(*result))
        {
            Point_2 ip = std::get<Point_2>(*result);

            Halfedge_const_handle he = *(singularSegmentHandle.originatingHalfEdge); 
            //singularArrangement.halfEdgePoints[he].push_back(ip);

            this->edgeRegionSegments[he->data()].push_back({regularSegmentHandle.originatingRegularEdgeId, true});
        }
    };

    this->edgeRegionSegments.resize(singularArrangement.arr.number_of_halfedges());

    CGAL::box_intersection_d(regularBoxes.begin(), regularBoxes.end(),
                              singularBoxes.begin(), singularBoxes.end(),
                              cb);

    // Sort all points along edges
    //
    //for (auto &[halfEdge, intersectingSegments] : edgeRegionSegments)
    for (auto halfEdge = singularArrangement.arr.halfedges_begin(); halfEdge != singularArrangement.arr.halfedges_end(); ++halfEdge) 
    {
        auto cmp = [this, &halfEdge, &tetMesh, &singularArrangement] (std::pair<int, int> intersectingSegmentA, std::pair<int, int> intersectingSegmentB)
            {
                const int segmentIdA = intersectingSegmentA.first;
                const int segmentIdB = intersectingSegmentB.first;

                if (segmentIdA == segmentIdB)
                {
                    return false;
                }

                const Segment_2 segment(halfEdge->source()->point(), halfEdge->target()->point());

                const std::array<int, 2> edgeA = tetMesh.edges[segmentIdA];
                const Segment_2 segmentA(singularArrangement.arrangementPoints[edgeA[0]], singularArrangement.arrangementPoints[edgeA[1]]);

                const std::array<int, 2> edgeB = tetMesh.edges[segmentIdB];
                const Segment_2 segmentB(singularArrangement.arrangementPoints[edgeB[0]], singularArrangement.arrangementPoints[edgeB[1]]);


                std::optional<std::variant<Point_2, Segment_2>> resultA = CGAL::intersection(segment, segmentA);
                std::optional<std::variant<Point_2, Segment_2>> resultB = CGAL::intersection(segment, segmentB);

                if (false == (resultA && std::holds_alternative<Point_2>(*resultA) && resultB && std::holds_alternative<Point_2>(*resultB)))
                {
                    throw std::runtime_error("OVERLAPPING OR NON-CROSSING SEGMENTS.");
                }

                Point_2 pA = std::get<Point_2>(*resultA);
                Point_2 pB = std::get<Point_2>(*resultB);

                K::Vector_2 v = segment.target() - segment.source();   // vector along segment

                K::Vector_2 wA = pA - segment.source();              // vector from source to point
                K::FT tA = (wA * v) / v.squared_length();       // exact dot / squared length

                K::Vector_2 wB = pB - segment.source();              // vector from source to point
                K::FT tB = (wB * v) / v.squared_length();       // exact dot / squared length

                return tA < tB;
            };

        auto &intersectingSegments = this->edgeRegionSegments[halfEdge->data()];

        std::sort(intersectingSegments.begin(), intersectingSegments.end(), cmp);
    }
}


void ReebSpace2::computeEdgeRegionMinusPlusTriangles(const TetMesh &tetMesh, Arrangement &singularArrangement)
{
    //for (auto &[halfEdge, intersectingSegments] : edgeRegionSegments)
    for (auto halfEdge = singularArrangement.arr.halfedges_begin(); halfEdge != singularArrangement.arr.halfedges_end(); ++halfEdge) 
    //for (const auto &[halfEdge, segmentIdsMap] : edgeRegionSegmentsMap)
    {
        auto &intersectingSegments = this->edgeRegionSegments[halfEdge->data()];

        // Each vertex is guaratneed to be in the singular arrangement
        const Point_2 &a = halfEdge->source()->point();
        const Point_2 &b = halfEdge->target()->point();

        //std::vector<std::vector<int>> minusTrianglesAll;
        //std::vector<std::vector<int>> plusTrianglesAll;

        // Set up handles so don't have to pay for access
        //std::vector<std::vector<int>> &edgeRegionMinusTrianglesHandle = edgeRegionMinusTriangles[halfEdge];
        //std::vector<std::vector<int>> &edgeRegionPlusTrianglesHandle = edgeRegionPlusTriangles[halfEdge];

        for (std::pair<int, bool> &intersectingSegment : intersectingSegments)
        {
            const int segmentId = intersectingSegment.first;
            const std::array<int, 2> edge = tetMesh.edges.at(segmentId);
            const int segmentSourceId = edge[0];
            const Point_2 &c = singularArrangement.arrangementPoints[segmentSourceId];

            assert(CGAL::orientation(a, b, c) != CGAL::COLLINEAR);

            //std::vector<int> minusTriangles;
            //std::vector<int> plusTriangles;

            //
            // If the orientation is the other way around note it
            //   (regular segment)
            //          d 
            //           \
            //            \
            // a ----------\--------- b (singular segment)
            //              \
            //               \
            //                c
            //
            if (CGAL::orientation(a, b, c) == CGAL::RIGHT_TURN)
            {
                intersectingSegment.second = false;
                //minusTriangles = tetMesh.upperStarTriangles.at(edge);
                //plusTriangles = tetMesh.lowerStarTriangles.at(edge);
            }

            //
            //   (regular segment)
            //          c 
            //           \
            //            \
            // a ----------\--------- b (singular segment)
            //              \
            //               \
            //                d
            //
            //else if (CGAL::orientation(a, b, c) == CGAL::LEFT_TURN)
            //{
                //minusTriangles = tetMesh.lowerStarTriangles.at(edge);
                //plusTriangles = tetMesh.upperStarTriangles.at(edge);
            //}

            //// Non-robust predicate issue
            //else
            //{
                //throw std::runtime_error("Input is degenerate, three points are collienar.");
            //}

            //printf("Minus Triangles:\n");
            //for (auto triangleId : minusTriangles)
            //{
                //printf("%d\n", triangleId);
            //}

            //printf("Plus Triangles:\n");
            //for (auto triangleId : plusTriangles)
            //{
                //printf("%d\n", triangleId);
            //}

            
            //edgeRegionMinusTrianglesHandle.push_back(minusTriangles);
            //edgeRegionPlusTrianglesHandle.push_back(plusTriangles);
        }

        auto &intersectingSegmentsTwin = this->edgeRegionSegments[halfEdge->twin()->data()];

        if (intersectingSegments.size() > 0 && intersectingSegmentsTwin.size() == 0)
        {
            // In the twin edge reverse the order
            intersectingSegmentsTwin.assign(intersectingSegments.rbegin(), intersectingSegments.rend());

            // Reverse the direction of travel as well
            for (std::pair<int, bool> &intersectingSegment : intersectingSegmentsTwin)
            {
                intersectingSegment.second = !intersectingSegment.second;
            }

        }





        //printf("\n\nThe intersecting segments are :\n");
        //for (const auto &[edgeId, isDirectionLowerToUpper] : edgeRegionSegments[halfEdge])
        //{
            //printf("%d %d\n", edgeId, isDirectionLowerToUpper);


            //const std::vector<int> &minusTriangles = tetMesh.getMinusTriangles(edgeId, isDirectionLowerToUpper);
            //const std::vector<int> &plusTriangles = tetMesh.getPlusTriangles(edgeId, isDirectionLowerToUpper);

            //printf("Minus Triangles:\n");
            //for (auto triangleId : minusTriangles)
            //{
                //printf("%d\n", triangleId);
            //}

            //printf("Plus Triangles:\n");
            //for (auto triangleId : plusTriangles)
            //{
                //printf("%d\n", triangleId);
            //}
        //}

        //printf("The intersecting segments of the twin :\n");
        //for (std::pair<int, bool> &intersectingSegment : edgeRegionSegments[halfEdge->twin()])
        //{
            //printf("%d %d\n", intersectingSegment.first, intersectingSegment.second);
        //}
        //printf("\n\n");
        //printf("------------------------------------------------------------------------------------------------------------");
        //printf("\n\n");

    }
}





























bool ReebSpace2::ifSegmentInHalfEdgeRegion(Arrangement_2::Halfedge_around_vertex_const_circulator &halfEdgeCirculator, const Segment_2 &segment)
{
    assert(halfEdgeCirculator != nullptr); 

    const auto halfEdge = halfEdgeCirculator;

    auto halfEdgeCirculatorPrevious = halfEdgeCirculator;
    const auto halfEdgeNext = ++halfEdgeCirculatorPrevious;

    // The image of the vertex is o, the endpoints of the half-edges are a and b and the other endpoints of the segment is b
    // See bellow for pictures
    const Point_2 &o = halfEdge->target()->point();
    const Point_2 &a = halfEdge->source()->point();
    const Point_2 &b = segment.source() == o ? segment.target() : segment.source();
    const Point_2 &c = halfEdgeNext->source()->point();

    //std::cout << "o = " << o << "\na = " << a << "\nb = " << b << "\nc = " << c << std::endl;

    assert(CGAL::orientation(o, a, c) != CGAL::COLLINEAR);
    assert(CGAL::orientation(o, a, b) != CGAL::COLLINEAR);

    // Case 1. 
    //
    //      a
    //      |   b
    //      |  /
    //      | /
    //      |/
    //      o-----c
    if (CGAL::orientation(o, a, c) == CGAL::RIGHT_TURN)
    {
        //printf("Returning left, left");
        //std::cout << "Return left, left " << CGAL::orientation(o, a, b) << ", " << CGAL::orientation(o, b, c) << std::endl;

        return CGAL::orientation(o, a, b) == CGAL::RIGHT_TURN && CGAL::orientation(o, b, c) == CGAL::RIGHT_TURN;
    }

    // Case 2. NOT the following
    //
    //       a
    //   b   |
    //    \  |
    //     \ |
    //      \|
    // c-----o
    else if (CGAL::orientation(o, a, c) == CGAL::LEFT_TURN)
    {
        //std::cout << "Return not right, right " << CGAL::orientation(o, a, b) << ", " << CGAL::orientation(o, b, c) << std::endl;
        return ! (CGAL::orientation(o, a, b) == CGAL::LEFT_TURN && CGAL::orientation(o, b, c) == CGAL::LEFT_TURN);
    }
    else
    {
        assert(false);
    }
}

Halfedge_const_handle ReebSpace2::getSegmentRegion(Vertex_const_handle &vertexHandle, const Segment_2 &segment)
{
    assert(segment.source() == vertexHandle->point() || segment.target() == vertexHandle->point());
    assert(false == vertexHandle->is_isolated());

    const auto first = vertexHandle->incident_halfedges();
    auto current = first;

    // Iterate in a CCW fashion to keep things consistent (faces are iterated in a CCW fashion too)
    do 
    {
        //std::cout << "  Halfedge from " << current->source()->point() << " to " << current->target()->point() << "\n";

        if (ifSegmentInHalfEdgeRegion(current, segment))
        {
            //std::cout << "Found its rightful place.\n";
            return current;
        }

    } while (++current != first);

    // Something has gone terribly wrong if we are here.
    assert(false);
}


bool ReebSpace2::compareRegularSegments(const Halfedge_const_handle &halfEdge, const Point_2& b, const Point_2& c)
{
    const Point_2 &o = halfEdge->target()->point();
    const Point_2 &a = halfEdge->source()->point();

    assert(CGAL::orientation(o, a, b) != CGAL::COLLINEAR);
    assert(CGAL::orientation(o, a, c) != CGAL::COLLINEAR);
    assert(CGAL::orientation(o, b, c) != CGAL::COLLINEAR);

    // Both are in the right half-plane
    if (CGAL::orientation(o, a, b) == CGAL::RIGHT_TURN)
    {
        if (CGAL::orientation(o, a, c) == CGAL::LEFT_TURN)
        {
            return true;
        }
        // Both are in the right half-plane
        else
        {
            return CGAL::orientation(o, b, c) == CGAL::RIGHT_TURN;
        }
    }
    // (o, a, b) Left Turn
    else
    {
        if (CGAL::orientation(o, a, c) == CGAL::RIGHT_TURN)
        {
            return false;
        }
        else
        {
            return CGAL::orientation(o, b, c) == CGAL::RIGHT_TURN;
        }
    }
}


void ReebSpace2::computeVertexRegionSegments(const TetMesh &tetMesh, Arrangement &singularArrangement)
{
    // @TODO
    // Remaking these is not efficient
    std::vector<MySegment_2> regularSegments;
    regularSegments.reserve(tetMesh.regularEdgesNumber);

    for (const auto &[edge, type] : tetMesh.edgeSingularTypes) 
    {
        if (type == 1)
        {
            const std::array<int, 2> &edgeConst = edge;

            regularSegments.emplace_back(
                    Segment_2(singularArrangement.arrangementPoints[edge[0]], singularArrangement.arrangementPoints[edge[1]]), 
                    std::nullopt, 
                    tetMesh.edgeIndices.at(edgeConst)
                    );
        }
    }

    this->vertexRegionSegments.resize(singularArrangement.arr.number_of_halfedges());

    for (const MySegment_2 &mySegment : regularSegments)
    {
        const Segment_2 &segment = mySegment.seg;

        //const int aIndex = singularArrangement.arrangementPointIndices.at(segment.source());
        //const int bIndex = singularArrangement.arrangementPointIndices.at(segment.target());
        //printf("-----------------------------------------------------  At segment [%d, %d].\n", aIndex, bIndex);
        //std::cout << "The current segment is from " << segment.source() << " to " << segment.target() << std::endl;

        const auto itSource = singularArrangement.arrangementPointHandles.find(segment.source());
        if (itSource != singularArrangement.arrangementPointHandles.end())
        {
            //std::cout << "The source point is : " << itSource->second->point() << " | " << singularArrangement.arrangementPointIndices.at(itSource->second->point()) << std::endl;
            Halfedge_const_handle regionHalfEdgeHandle = getSegmentRegion(itSource->second, segment);
            vertexRegionSegments[regionHalfEdgeHandle->data()].push_back({mySegment.originatingRegularEdgeId, true});

            //std::cout << "Added segment [" << tetMesh.edges[mySegment.originatingEdge][0] << ", " << tetMesh.edges[mySegment.originatingEdge][1] << "]\n";
            //std::cout << "Between " << singularArrangement.arrangementPoints[tetMesh.edges[mySegment.originatingEdge][0]] << " and " << singularArrangement.arrangementPoints[tetMesh.edges[mySegment.originatingEdge][0]] << std::endl;
        }


        const auto itTarget = singularArrangement.arrangementPointHandles.find(segment.target());
        if (itTarget != singularArrangement.arrangementPointHandles.end())
        {
            //std::cout << "The target point is : " << itTarget->second->point() << " | " <<  singularArrangement.arrangementPointIndices.at(itTarget->second->point()) << std::endl;
            Halfedge_const_handle regionHalfEdgeHandle = getSegmentRegion(itTarget->second, segment);
            vertexRegionSegments[regionHalfEdgeHandle->data()].push_back({mySegment.originatingRegularEdgeId, true});

            //std::cout << "Added segment [" << tetMesh.edges[mySegment.originatingEdge][0] << ", " << tetMesh.edges[mySegment.originatingEdge][1] << "]\n";
            //std::cout << "Between " << singularArrangement.arrangementPoints[tetMesh.edges[mySegment.originatingEdge][0]] << " and " << singularArrangement.arrangementPoints[tetMesh.edges[mySegment.originatingEdge][0]] << std::endl;
        }

        //printf("\n\n");
    }


    // Sort the regions
    //for (auto &[halfEdge, segmentIds] : vertexRegionSegments)
    for (auto halfEdge = singularArrangement.arr.halfedges_begin(); halfEdge != singularArrangement.arr.halfedges_end(); ++halfEdge) 
    {
        const Point_2 &vertex = halfEdge->target()->point();
        const int vertexMeshId = singularArrangement.arrangementPointIndices[vertex];

        auto cmp = [this, &halfEdge, &tetMesh, vertexMeshId, &singularArrangement] (std::pair<int, int> intersectingSegmentA, std::pair<int, int> intersectingSegmentB)
            {
                const int segmentIdA = intersectingSegmentA.first;
                const int segmentIdB = intersectingSegmentB.first;

                const std::array<int, 2> edgeA = tetMesh.edges[segmentIdA];
                const std::array<int, 2> edgeB = tetMesh.edges[segmentIdB];

                Point_2 b, c;
                if (edgeA[0] == vertexMeshId)
                    b = singularArrangement.arrangementPoints[edgeA[1]];
                else
                    b = singularArrangement.arrangementPoints[edgeA[0]];

                if (edgeB[0] == vertexMeshId)
                    c = singularArrangement.arrangementPoints[edgeB[1]];
                else
                    c = singularArrangement.arrangementPoints[edgeB[0]];

                return compareRegularSegments(halfEdge, b, c);
            };

        auto &segmentIds = this->vertexRegionSegments[halfEdge->data()];
        std::sort(segmentIds.begin(), segmentIds.end(), cmp);
    }
}


void ReebSpace2::computeVertexRegionMinusPlusTriangles(const TetMesh &tetMesh, Arrangement &singularArrangement)
{
    //for (auto &[halfEdge, intersectingSegments] : vertexRegionSegments)
    //for (auto &[halfEdge, intersectingSegments] : vertexRegionSegments)
    for (auto halfEdge = singularArrangement.arr.halfedges_begin(); halfEdge != singularArrangement.arr.halfedges_end(); ++halfEdge) 
    {
        auto &intersectingSegments = this->vertexRegionSegments[halfEdge->data()];
        // Each vertex is guaratneed to be in the singular arrangement
        const Point_2 &vertex = halfEdge->target()->point();
        const int vertexMeshId = singularArrangement.arrangementPointIndices[vertex];

        //auto &vertexRegionMinusTrianglesHandle = vertexRegionMinusTriangles[halfEdge];
        //auto &vertexRegionPlusTrianglesHandle = vertexRegionPlusTriangles[halfEdge];

        //std::set<int> minusTrianglesSet;
        //std::set<int> plusTrianglesSet;

        //std::vector<std::vector<int>> plusTrianglesAll;
        //std::vector<std::vector<int>> minusTrianglesAll;

        //for (const auto &[segmentId,  : segmentIds)
        for (std::pair<int, bool> &intersectingSegment : intersectingSegments)
        {
            //const std::array<int, 2> edge = tetMesh.edges.at(segmentId);
            //const int segmentSourceId = edge[0];
            //const Point_2 &c = singularArrangement.arrangementPoints[segmentSourceId];

            const int segmentId = intersectingSegment.first;
            const std::array<int, 2> edge = tetMesh.edges[segmentId];

            //std::vector<int> minusTriangles;
            //std::vector<int> plusTriangles;


            // Crossing ab in a CW direction goes from the lower to the upper star
            // Remember that we assume that a < b, so if a = v, then b is in the set BiggerThan(a).
            // The set BiggerThan(a) = {x \ in R^2 : x > a} has a closed boundary (|) above a.y and open boundary bellow (.) a.y.
            //
            //      |  b
            //      | /
            //      |/
            //      a
            //      .
            //      . 
            //      .
            //      
            if (edge[0] == vertexMeshId)
            {
                intersectingSegment.second = false;
                //minusTriangles = tetMesh.upperStarTriangles.at(edge);
                //plusTriangles = tetMesh.lowerStarTriangles.at(edge);
            }
            // Crossing ab in a CW direction goes from the upper to the lower star
            // Remember that we assume that a < b, so if b = v, then a is in the set SmallerThan(b).
            // The set SmallerThan(b) = {x \ in R^2 : x > a} has an open boundary (.) above b.y and closed boundary bellow (|) b.y
            //      .
            //      .
            //      .
            //      b
            //     /|
            //    / | 
            //   a  |
            //      
            //else if (edge[1] == vertexMeshId)
            //{
                //intersectingSegment.second = true;
                ////minusTriangles = tetMesh.lowerStarTriangles.at(edge);
                ////plusTriangles = tetMesh.upperStarTriangles.at(edge);
            //}
            //else
            //{
                //assert(false);
            //}

            //printf("Minus Triangles:\n");
            //for (auto triangleId : minusTriangles)
            //{
                //printf("%d\n", triangleId);
            //}

            //printf("Plus Triangles:\n");
            //for (auto triangleId : plusTriangles)
            //{
                //printf("%d\n", triangleId);
            //}

            //plusTrianglesAll.push_back(plusTriangles);
            //minusTrianglesAll.push_back(minusTriangles);

            //vertexRegionMinusTrianglesHandle.push_back(minusTriangles);
            //vertexRegionPlusTrianglesHandle.push_back(plusTriangles);

            //std::map<Halfedge_const_handle, std::vector<int>> vertexRegionMinusTriangles;
            //std::map<Halfedge_const_handle, std::vector<int>> vertexRegionPlusTriangles;

            //minusTrianglesSet.insert(minusTriangles.begin(), minusTriangles.end());
            //plusTrianglesSet.insert(plusTriangles.begin(), plusTriangles.end());
        }

        // Cancel out the plus/minus triangles and write to a vector
        //std::set_difference(
                //minusTrianglesSet.begin(), minusTrianglesSet.end(),
                //plusTrianglesSet.begin(), plusTrianglesSet.end(),
                //std::back_inserter(vertexRegionMinusTriangles[halfEdge]));

        //std::set_difference(
                //plusTrianglesSet.begin(), plusTrianglesSet.end(),
                //minusTrianglesSet.begin(), minusTrianglesSet.end(),
                //std::back_inserter(vertexRegionPlusTriangles[halfEdge]));


            


        //printf("\n\nThe intersecting segments are :\n");
        ////for (const auto &[edgeId, isDirectionLowerToUpper] : vertexRegionSegments[halfEdge])
        //for (int i = 0 ; i < vertexRegionSegments[halfEdge].size() ; i++)
        //{
            //const auto [edgeId, isDirectionLowerToUpper] = vertexRegionSegments[halfEdge][i];

            //printf("------------------------------- %d %d\n", edgeId, isDirectionLowerToUpper);

            //const std::vector<int> &minusTriangles = tetMesh.getMinusTriangles(edgeId, isDirectionLowerToUpper);
            //const std::vector<int> &plusTriangles = tetMesh.getPlusTriangles(edgeId, isDirectionLowerToUpper);


            //printf("Minus Triangles:\n");
            //for (auto triangleId : minusTriangles)
            //{
                //printf("%d\n", triangleId);
            //}

            //printf("Plus Triangles:\n");
            //for (auto triangleId : plusTriangles)
            //{
                //printf("%d\n", triangleId);
            //}

            //if (minusTriangles != minusTrianglesAll[i])
            //{
                //throw std::runtime_error("Minus Triangles differ!");
            //}

            //if (plusTriangles != plusTrianglesAll[i])
            //{
                //throw std::runtime_error("Minus Triangles differ!");
            //}
        //}
        //printf("\n\n");
        //printf("------------------------------------------------------------------------------------------------------------");
        //printf("\n\n");
    }
}








void ReebSpace2::computeEdgeCrossingMinusPlusTriangles(const TetMesh &tetMesh, Arrangement &singularArrangement)
{
    this->edgeCrossingSegments.resize(singularArrangement.arr.number_of_halfedges());

    for (auto he = singularArrangement.arr.halfedges_begin(); he != singularArrangement.arr.halfedges_end(); ++he)
    {
        //std::cout << "Edge id is " << he->data() << std::endl;
        // If we have computed this for the twin, just swap them around
        //if (edgeCrossingSegments.contains(he->twin()))
        //{
        //edgeCrossingSegments[he].first = edgeCrossingSegments[he->twin()].first;
        //edgeCrossingSegments[he].second = !edgeCrossingSegments[he->twin()].second;
        //}
        //else
        //{
        const Segment_2 &segment = *singularArrangement.arr.originating_curves_begin(he);
        //std::cout << "Half-edge   from: " << he->source()->point() << " to " << he->target()->point() << std::endl;
        //std::cout << "Source-edge from: " << segment.source() << " to " << segment.target() << std::endl;

        const int aIndex = singularArrangement.arrangementPointIndices.at(segment.source());
        const int bIndex = singularArrangement.arrangementPointIndices.at(segment.target());

        // Sanity check
        assert(aIndex < bIndex);

        const std::array<int, 2> edge = {aIndex, bIndex};

        // Check to see if the segment and half edge have the same orientation
        const bool isSegmentLeftToRight = segment.source() < segment.target(); 
        const bool isCurrentHalfEdgeLeftToRight = (he->direction() == CGAL::ARR_LEFT_TO_RIGHT);

        const int edgeId = tetMesh.edgeIndices.at(edge);

        edgeCrossingSegments[he->data()].first = edgeId;

        // The half edge has the same direction as the original edge
        if (isSegmentLeftToRight == isCurrentHalfEdgeLeftToRight)
        {
            edgeCrossingSegments[he->data()].second = false;
            //edgeCrossingMinusTriangles[he] = tetMesh.upperStarTriangles.at(edge);
            //edgeCrossingPlusTriangles[he] = tetMesh.lowerStarTriangles.at(edge);
        }
        else
        {
            edgeCrossingSegments[he->data()].second = true;
            //edgeCrossingMinusTriangles[he] = tetMesh.lowerStarTriangles.at(edge);
            //edgeCrossingPlusTriangles[he] = tetMesh.upperStarTriangles.at(edge);
        }
        //}
    }
}











bool ReebSpace2::areHalfEdgeRegionMapsEqual(const std::map<Halfedge_const_handle, std::set<int>>& a, const std::map<Halfedge_const_handle, std::set<int>>& b)
{
    if (a.size() != b.size())
    {
        std::cerr << "Size is not equal.\n";
        return false;
    }

    for (const auto& [key, setA] : a)
    {
        auto itB = b.find(key);
        if (itB == b.end())
        {
            std::cerr << "Second map does not have a key.\n";
            return false;
        }

        const std::set<int>& setB = itB->second;
        if (setA.size() != setB.size())
        {
            std::cerr << "Sets for the same key do not match in size.\n";
            return false;
        }

        if (setA != setB)
        {
            std::cerr << "Set elements for the same key do not match.\n";
            return false;
        }
    }

    return true;
}



void ReebSpace2::unitTest(const TetMesh &tetMesh, Arrangement &singularArrangement, Arrangement &regularArrangement)
{
    std::vector<MySegment_2> singularSegments;
    singularSegments.reserve(singularArrangement.arr.number_of_edges());

    std::vector<Box> singularBoxes;
    singularBoxes.reserve(singularArrangement.arr.number_of_edges());

    for (auto he = singularArrangement.arr.halfedges_begin(); he != singularArrangement.arr.halfedges_end(); ++he)
    {
        if (he < he->twin())
        {
            singularSegments.emplace_back(Segment_2(he->source()->point(), he->target()->point()), he);
            singularBoxes.emplace_back(singularSegments.back().seg.bbox(), &singularSegments.back());
        }
    }

    std::vector<MySegment_2> regularSegments;
    regularSegments.reserve(tetMesh.regularEdgesNumber);

    std::vector<Box> regularBoxes; 
    regularBoxes.reserve(tetMesh.regularEdgesNumber);

    for (const auto &[edge, type] : tetMesh.edgeSingularTypes) 
    {
        if (type == 1)
        {
            const std::array<int, 2> &edgeConst = edge;

            regularSegments.emplace_back(
                    Segment_2(singularArrangement.arrangementPoints[edge[0]], singularArrangement.arrangementPoints[edge[1]]), 
                    std::nullopt, 
                    tetMesh.edgeIndices.at(edgeConst)
                    );
            regularBoxes.emplace_back(regularSegments.back().seg.bbox(), &regularSegments.back());
        }
    }









    // Make sure we have compute all intersections correctly
    std::map<Halfedge_const_handle, std::set<int>> halfEdgeEdgeRegionSegmentsUnitTest;

    for (const auto &myRegularSegment : regularSegments)
    {
        for (const auto &mySingularSegment : singularSegments)
        {
            if (doSegmentEndpointsOverlap(myRegularSegment.seg, mySingularSegment.seg)) { continue; }

            std::optional<std::variant<Point_2, Segment_2>> result = CGAL::intersection(myRegularSegment.seg, mySingularSegment.seg);

            if (result && std::holds_alternative<Point_2>(*result))
            {
                Point_2 ip = std::get<Point_2>(*result);

                Halfedge_const_handle he = *(mySingularSegment.originatingHalfEdge); 

                halfEdgeEdgeRegionSegmentsUnitTest[he].insert(myRegularSegment.originatingRegularEdgeId);
            }
        }
    }

    //assert(areHalfEdgeRegionMapsEqual(edgeRegionSegments, halfEdgeEdgeRegionSegmentsUnitTest));



    // halfEdgeVertexRegionSegments;


    std::map<Halfedge_const_handle, std::set<int>> halfEdgeVertexRegionSegmentsUnitTest;

    // Iterate over the halfedges of every vertex
    for (auto v = regularArrangement.arr.vertices_begin(); v != regularArrangement.arr.vertices_end(); ++v)
    {
        if (v->is_isolated() || false == regularArrangement.arrangementPointIndices.contains(v->point())) 
        {
            continue;
        }

        //std::cout << "At vertex " << v->point() << std::endl;

        Arrangement_2::Halfedge_around_vertex_const_circulator first, curr;
        first = curr = v->incident_halfedges();

        bool singularEdgeFound = false;

        // Go to the first singular edge you find.
        do 
        {
            // If the edge is singular
            const Segment_2 &segment = *regularArrangement.arr.originating_curves_begin(curr);

            const int aIndex = regularArrangement.arrangementPointIndices.at(segment.source());
            const int bIndex = regularArrangement.arrangementPointIndices.at(segment.target());
            const int edgeType = tetMesh.edgeSingularTypes.at({aIndex, bIndex});
            if  (edgeType != 1)
            {
                first  = curr;
                singularEdgeFound = true;
                break;
            }

        } while (++curr != first);

        if (false == singularEdgeFound) { continue; }


        //std::cout << "Starting half-edge at " << curr->source()->point() << " to " << curr->target()->point() << "\n";

        Arrangement_2::Halfedge_const_handle currentHalfEdge;
        Arrangement_2::Halfedge_const_handle currentSingularHalfEdge;
        do 
        {
            const Segment_2 &segment = *regularArrangement.arr.originating_curves_begin(curr);
            
            const int aIndex = regularArrangement.arrangementPointIndices.at(segment.source());
            const int bIndex = regularArrangement.arrangementPointIndices.at(segment.target());
            const int edgeType = tetMesh.edgeSingularTypes.at({aIndex, bIndex});

            // If the edge is singular
            if  (edgeType != 1)
            {
                currentHalfEdge = curr;

                bool singularHalfEdgeFound = false;

                for (auto he = singularArrangement.arr.halfedges_begin(); he != singularArrangement.arr.halfedges_end(); ++he)
                {
                    const Segment_2 &singularSegment = *singularArrangement.arr.originating_curves_begin(he);

                    if (he->target()->point() == v->point() && segment == singularSegment)
                    {

                        currentSingularHalfEdge = he;
                        singularHalfEdgeFound = true;
                    }
                }

                assert(singularHalfEdgeFound);





                //std::cout << "---- Singular region start at " << curr->source()->point() << " to " << curr->target()->point() << "\n";
                //std::cout << "---- Singular region start at singular " << currentSingularHalfEdge->source()->point() << " to " << currentSingularHalfEdge->target()->point() << "\n";
                //std::cout << std::endl;
            }
            else
            {

                //std::cout << "---------- Adding regular edge " << curr->source()->point() << " to " << curr->target()->point() << "\n";

                //int nonTempSegmentId;

                //bool nonTempSegmentFound = false;
                //for (const MySegment_2 &mySegment : regularSegments)
                //{
                    //if (mySegment.seg.source() == segment.source() && mySegment.seg.target() == segment.target())
                    //{
                        //nonTempSegment = &mySegment.seg;
                        //nonTempSegmentFound = true;
                    //}
                //}

                //assert (nonTempSegmentFound);





                const int edgeId = tetMesh.edgeIndices.at({aIndex, bIndex});


                halfEdgeVertexRegionSegmentsUnitTest[currentSingularHalfEdge].insert(edgeId);
                // The current half edge is equal to which half-edge in the singular arrangement?


            }


        } while (++curr != first);
    }


    //assert(areHalfEdgeRegionMapsEqual(vertexRegionSegments, halfEdgeVertexRegionSegmentsUnitTest));
}








bool segmentsOverlap(const Segment_2& s1, const Segment_2& s2)
{
    auto result = CGAL::intersection(s1, s2);

    if (!result)
        return false;

    // Check what type of intersection we got
    if (const Segment_2* overlap = std::get_if<Segment_2>(&*result))
    {
        // Segments overlap in a positive-length segment
        return true;
    }

    // Point intersection is not considered overlap
    return false;
}

void ReebSpace2::assignHalfEdgePseudoSingular(const TetMesh &tetMesh, Arrangement &singularArrangement)
{
    this->isHalfEdgePseudoSingular.resize(singularArrangement.arr.number_of_halfedges());

    // @TODO Find a way to parallize this
    //#pragma omp parallel for schedule(dynamic)
    for (auto halfEdge = singularArrangement.arr.halfedges_begin(); halfEdge != singularArrangement.arr.halfedges_end(); ++halfEdge) 
    {
        const Segment_2 &segment = *singularArrangement.arr.originating_curves_begin(halfEdge);
        const std::array<int, 2> edge = {
            singularArrangement.arrangementPointIndices.at(segment.source()), 
            singularArrangement.arrangementPointIndices.at(segment.target())
        };

        if (tetMesh.edgeSingularTypes.at(edge) == -1)
        {
            this->isHalfEdgePseudoSingular[halfEdge->data()] = true;
        }
        else
        {
            this->isHalfEdgePseudoSingular[halfEdge->data()] = false;
        }
    }

}

std::pair<Halfedge_const_handle, Halfedge_const_handle>  findHalfEdges(const Halfedge_const_handle &singularHalfEdge,  Arrangement &singularArrangement, Arrangement &regularArrangement)
{
    std::pair<Halfedge_const_handle, Halfedge_const_handle> foundHalfEdges;

    const Segment_2 singularHalfEdgeSegment = singularHalfEdge->curve();

    Vertex_const_handle sourceVertex = regularArrangement.arrangementPointHandles.at(singularHalfEdge->source()->point());
    Arrangement_2::Halfedge_around_vertex_const_circulator first, curr;
    first = curr = sourceVertex->incident_halfedges();
    do 
    {
        const Segment_2 regularHalfEdgeSegment = curr->curve();

        if (segmentsOverlap(singularHalfEdgeSegment, regularHalfEdgeSegment))
        {
            foundHalfEdges.first = curr->twin();
        }
    } while (++curr != first);


    Vertex_const_handle targetVertex = regularArrangement.arrangementPointHandles.at(singularHalfEdge->target()->point());
    first = curr = targetVertex->incident_halfedges();
    do 
    {
        const Segment_2 regularHalfEdgeSegment = curr->curve();

        if (segmentsOverlap(singularHalfEdgeSegment, regularHalfEdgeSegment))
        {
            foundHalfEdges.second = curr;
        }
    } while (++curr != first);


    //assert (foundFaces.first == Face_const_handle() || foundFaces.second != Face_const_handle());
    assert (foundHalfEdges.first == Halfedge_const_handle() || foundHalfEdges.second != Halfedge_const_handle());

    return foundHalfEdges;
}


bool ReebSpace2::unitTestComparePreimageGraphs(const TetMesh &tetMesh, Arrangement &singularArrangement, Arrangement &regularArrangement, ReebSpace &rs)
{
    for (const auto& [halfEdge, preimageGraphs] : this->preimageGraphsCached)
    {
        //std::cout << "Checking graph equality..." << std::endl;
        std::pair<Halfedge_const_handle, Halfedge_const_handle> foundHalfEdges =  findHalfEdges(halfEdge, singularArrangement, regularArrangement);
        std::pair<Face_const_handle, Face_const_handle> foundFaces =  {foundHalfEdges.first->face(), foundHalfEdges.second->face()};
        std::pair<int, int> foundFaceIds =  {regularArrangement.arrangementFacesIdices.at(foundFaces.first), regularArrangement.arrangementFacesIdices.at(foundFaces.second)};

        std::pair<PreimageGraph, PreimageGraph> &singularPreimageGraphs = this->preimageGraphsCached[halfEdge];
        std::pair<DisjointSet<int>, DisjointSet<int>> regularPreimageGraphs = {rs.preimageGraphs[foundFaceIds.first], rs.preimageGraphs[foundFaceIds.second]};

        //assert(singularPreimageGraphs.first == regularPreimageGraphs.first);
        //assert(singularPreimageGraphs.second == regularPreimageGraphs.second);


        if (false == singularPreimageGraphs.first.areEqual(regularPreimageGraphs.first) || false == singularPreimageGraphs.second.areEqual(regularPreimageGraphs.second))
        {
            printf("\n------------------------------------------------------------------------------------\n");
            std::cout << "SingularHalf-edge is [" << halfEdge->source()->point() << "] -> [" << halfEdge->target()->point() << "]";
            printf("The IDs are %d -> %d\n", singularArrangement.arrangementPointIndices[halfEdge->source()->point()], singularArrangement.arrangementPointIndices[halfEdge->target()->point()]);
            printf("\n------------------------------------------------------------------------------------\n");


            std::cout << "First singular preimage graph is: \n";
            singularPreimageGraphs.first.printByRoot();

            std::cout << "Second singular preimage graph is: \n";
            singularPreimageGraphs.second.printByRoot();


            printf("\n------------------------------------------------------------------------------------\n");
            std::cout << "Regular-edge 1 is [" << foundHalfEdges.first->source()->point() << "] -> [" << foundHalfEdges.first->target()->point() << "]\n";
            std::cout << "Regular-edge 2 is [" << foundHalfEdges.second->source()->point() << "] -> [" << foundHalfEdges.second->target()->point() << "]";
            printf("\n------------------------------------------------------------------------------------\n");

            std::cout << "First regular preimage graph is: \n";
            regularPreimageGraphs.first.printByRoot();
            //regularPreimageGraphs.first.print([&](const int &triangleId) {
                    //io::printTriangle(tetMesh, triangleId);
                    //});

            std::cout << "Second regular preimage graph is: \n";
            regularPreimageGraphs.second.printByRoot();
            //regularPreimageGraphs.second.print([&](const int &triangleId) {
                    //io::printTriangle(tetMesh, triangleId);
                    //});

            printf("\n\n\n");

            return false;
            //throw std::runtime_error("Preimage graphs do not match!");
        }
    }
    return true;
}
