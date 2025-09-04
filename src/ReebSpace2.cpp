#include "./CGALTypedefs.h"

#include <cassert>
#include <utility>

#include "./io.h"
#include "./Timer.h"
#include "./ReebSpace2.h"
#include "./DisjointSet.h"
#include "./PreimageGraph.h"

PreimageGraph computePreimageGraph(const TetMesh &tetMesh, const std::vector<int> &minusTriangles, const std::vector<int> &plusTriangles, const PreimageGraph &preimageGraphPrevious)
{
    std::set<int> preimageGraphTriangles;
    for (const auto &triangleId : preimageGraphPrevious.triangleIds)
    {
        preimageGraphTriangles.insert(triangleId);
    }

    for (const auto &triangle: minusTriangles)
    {
        preimageGraphTriangles.erase(triangle);
    }

    for (const auto &triangle: plusTriangles)
    {
        preimageGraphTriangles.insert(triangle);
    }

    PreimageGraph pg(preimageGraphTriangles);
    pg.computeConnectedComponents(tetMesh);

    return pg;
}

// The halfEdge is in the twin face, it's second is our initial preimage graph
void ReebSpace2::loopFace(const TetMesh &tetMesh, const Halfedge_const_handle &initialHalfEdge, std::queue<Halfedge_const_handle> &traversalQueue, std::set<Face_const_handle> &visited)
{
    // We assume the the first one has already been seeded
    //
    //this->preimageGraphs[initialHalfEdge].first = computePreimageGraph(
            //tetMesh, 
            //this->edgeCrossingMinusTriangles[initialHalfEdge->twin()], 
            //this->edgeCrossingPlusTriangles[initialHalfEdge->twin()], 
            //this->preimageGraphs[initialHalfEdge->twin()].second
            //);

    this->preimageGraphs[initialHalfEdge].second = computePreimageGraph(
            tetMesh, 
            this->edgeRegionMinusTriangles[initialHalfEdge], 
            this->edgeRegionPlusTriangles[initialHalfEdge], 
            this->preimageGraphs[initialHalfEdge].first
            );

    //printf("\n------------------------------------------------------------------------------------\n");
    //std::cout << "Half-edge is [" << initialHalfEdge->source()->point() << "] -> [" << initialHalfEdge->target()->point() << "]";
    //printf("\n------------------------------------------------------------------------------------\n");
    //std::cout << "First preimage graph is: \n";
    //this->preimageGraphs[initialHalfEdge].first.print([&](const int &triangleId) {
            //io::printTriangle(tetMesh, triangleId);
            //});

    //std::cout << "Second preimage graph is: \n";
    //this->preimageGraphs[initialHalfEdge].second.print([&](const int &triangleId) {
            //io::printTriangle(tetMesh, triangleId);
            //});

    // Go the the next one
    Halfedge_const_handle currentHalfEdge = initialHalfEdge->next();

    do
    {
        Halfedge_const_handle previousHalfEdge = currentHalfEdge->prev();

        auto &currentPreimageGraphPair = this->preimageGraphs[currentHalfEdge];

        currentPreimageGraphPair.first = computePreimageGraph(
                tetMesh, 
                this->vertexRegionMinusTriangles[previousHalfEdge], 
                this->vertexRegionPlusTriangles[previousHalfEdge], 
                this->preimageGraphs[previousHalfEdge].second
                );

        currentPreimageGraphPair.second = computePreimageGraph(
                tetMesh, 
                this->edgeRegionMinusTriangles[currentHalfEdge], 
                this->edgeRegionPlusTriangles[currentHalfEdge], 
                currentPreimageGraphPair.first
                );


        if (false == visited.contains(currentHalfEdge->twin()->face()))
        {
            //std::cout<< "----------------------------------------------------------- ADDIG" << std::endl;
            //printf("\n------------------------------------------------------------------------------------\n");
            //std::cout << "Queue Half-edge is [" << currentHalfEdge->twin()->source()->point() << "] -> [" << currentHalfEdge->twin()->target()->point() << "]";
            //printf("\n------------------------------------------------------------------------------------\n");


            // Seed the first half-edge in the twin
            this->preimageGraphs[currentHalfEdge->twin()].first = computePreimageGraph(
                    tetMesh, 
                    this->edgeCrossingMinusTriangles[currentHalfEdge], 
                    this->edgeCrossingPlusTriangles[currentHalfEdge], 
                    currentPreimageGraphPair.second
                    );

            traversalQueue.push(currentHalfEdge->twin());
            visited.insert(currentHalfEdge->twin()->face());

        }




        //printf("\n------------------------------------------------------------------------------------\n");
        //std::cout << "Half-edge is [" << currentHalfEdge->source()->point() << "] -> [" << currentHalfEdge->target()->point() << "]";
        //printf("\n------------------------------------------------------------------------------------\n");

        //std::cout << "First preimage graph is: \n";
        //this->preimageGraphs[currentHalfEdge].first.print([&](const int &triangleId) {
                //io::printTriangle(tetMesh, triangleId);
                //});

        //std::cout << "Second preimage graph is: \n";
        //this->preimageGraphs[currentHalfEdge].second.print([&](const int &triangleId) {
                //io::printTriangle(tetMesh, triangleId);
                //});



        currentHalfEdge = currentHalfEdge->next();



    } while (currentHalfEdge != initialHalfEdge);

     
}


void ReebSpace2::traverse(const TetMesh &tetMesh, Arrangement &singularArrangement)
{
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
    std::set<Face_const_handle> visited;

    // Seed the first face
    this->preimageGraphs[startingHalfedge->twin()].first = computePreimageGraph(
            tetMesh, 
            this->edgeCrossingMinusTriangles[startingHalfedge], 
            this->edgeCrossingPlusTriangles[startingHalfedge], 
            PreimageGraph()
            );


    traversalQueue.push(startingHalfedge->twin());
    visited.insert(startingHalfedge->face());

    // Make sure the outer face is visited as well, no need to go back
    visited.insert(startingHalfedge->twin()->face());


    while (false == traversalQueue.empty())
    {
        Halfedge_const_handle currentHalfEdge = traversalQueue.front();
        traversalQueue.pop();

        // Process neighbours and queue them
        loopFace(tetMesh, currentHalfEdge, traversalQueue, visited);

        // Add to the correspondence graph
        //const std::vector<std::pair<int, int>> componentRepresentatives = this->preimageGraphs[currentHalfEdge].first.getComponentRepresentatives();

        //for (const auto &[componentId, triangleId] : componentRepresentatives)
        //{
            //this->correspondenceGraph.addElement({singularArrangement.arrangementFacesIdices[currentHalfEdge->face()], componentId});
        //}

        //// Clear the preimage graphs in the currect face.
        Halfedge_const_handle iterator = currentHalfEdge;
        do
        {
            this->preimageGraphsCached[iterator].first = this->preimageGraphs[iterator].first;
            this->preimageGraphsCached[iterator].second = this->preimageGraphs[iterator].second;

            this->preimageGraphs[iterator].first.clear();
            this->preimageGraphs[iterator].second.clear();
            iterator = iterator->next();

        } while (iterator != currentHalfEdge);
    }

    //printf("The correspondence graphs has %ld vertices.\n", this->correspondenceGraph..size());

    //preimageGraphFirst.print([&](const int &triangleId) {
            //io::printTriangle(tetMesh, triangleId);
            //});

    //preimageGraphSecond.print([&](const int &triangleId) {
            //io::printTriangle(tetMesh, triangleId);
            //});

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

bool ReebSpace2::doSegmentEndpointsOverlap(const Segment_2 &s1, const Segment_2 &s2)
{
    return 
        s1.source() == s2.source() || s1.source() == s2.target() ||
        s1.target() == s2.source() || s1.target() == s2.target();
}



void ReebSpace2::computeEdgeRegionSegments(const TetMesh &tetMesh, Arrangement &singularArrangement)
{
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

    int intersections = 0;
    int falsePositives = 0;

    auto cb = [&](const Box& regularSegmentBox, const Box& singularSegmentBox) {
        
        const auto &regularSegmentHandle = *regularSegmentBox.handle();
        const auto &singularSegmentHandle  = *singularSegmentBox.handle();

        if (doSegmentEndpointsOverlap(regularSegmentHandle.seg, singularSegmentHandle.seg)) { return; }

        std::optional<std::variant<Point_2, Segment_2>> result = CGAL::intersection(regularSegmentHandle.seg, singularSegmentHandle.seg);


        if (result && std::holds_alternative<Point_2>(*result))
        {
            Point_2 ip = std::get<Point_2>(*result);
            intersections++;

            Halfedge_const_handle he = *(singularSegmentHandle.originatingHalfEdge); 
            singularArrangement.halfEdgePoints[he].push_back(ip);

            edgeRegionSegments[he].insert(regularSegmentHandle.originatingEdge);

            //singularArrangement.arr.split_edge(*(singularSegmentHandle.originatingHalfEdge), ip);

            //const int aIndex = singularArrangement.arrangementPointIndices.at(singularSegmentHandle.seg.source());
            //const int bIndex = singularArrangement.arrangementPointIndices.at(singularSegmentHandle.seg.target());


            //intersectionPoints[{aIndex, bIndex}].push_back(ip);

            //std::cout << "Intersection point: " << ip << "\n";
        }
        else
        {
            //std::cout << "No intersection.\n";
            falsePositives++;
        }


        //if (CGAL::do_intersect(regularSegmentHandle, singularSegmentHandle)) 
        //{
            //intersections++;
        //}
        //else
        //{
            //falsePositives++;
        //}
    };

    CGAL::box_intersection_d(regularBoxes.begin(), regularBoxes.end(),
                              singularBoxes.begin(), singularBoxes.end(),
                              cb);




    printf("There were %d intersections with %d false positives.\n", intersections, falsePositives);



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

    // Case 2. 
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



void ReebSpace2::computeVertexRegionSegments(const TetMesh &tetMesh, Arrangement &singularArrangement)
{
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
            vertexRegionSegments[regionHalfEdgeHandle].insert(mySegment.originatingEdge);

            //std::cout << "Added segment [" << tetMesh.edges[mySegment.originatingEdge][0] << ", " << tetMesh.edges[mySegment.originatingEdge][1] << "]\n";
            //std::cout << "Between " << singularArrangement.arrangementPoints[tetMesh.edges[mySegment.originatingEdge][0]] << " and " << singularArrangement.arrangementPoints[tetMesh.edges[mySegment.originatingEdge][0]] << std::endl;
        }


        const auto itTarget = singularArrangement.arrangementPointHandles.find(segment.target());
        if (itTarget != singularArrangement.arrangementPointHandles.end())
        {
            //std::cout << "The target point is : " << itTarget->second->point() << " | " <<  singularArrangement.arrangementPointIndices.at(itTarget->second->point()) << std::endl;
            Halfedge_const_handle regionHalfEdgeHandle = getSegmentRegion(itTarget->second, segment);
            vertexRegionSegments[regionHalfEdgeHandle].insert(mySegment.originatingEdge);

            //std::cout << "Added segment [" << tetMesh.edges[mySegment.originatingEdge][0] << ", " << tetMesh.edges[mySegment.originatingEdge][1] << "]\n";
            //std::cout << "Between " << singularArrangement.arrangementPoints[tetMesh.edges[mySegment.originatingEdge][0]] << " and " << singularArrangement.arrangementPoints[tetMesh.edges[mySegment.originatingEdge][0]] << std::endl;
        }

        //printf("\n\n");
    }


}






void ReebSpace2::computeEdgeRegionMinusPlusTriangles(const TetMesh &tetMesh, Arrangement &singularArrangement)
{
    for (const auto [halfEdge, segmentIds] : edgeRegionSegments)
    {
        // Each vertex is guaratneed to be in the singular arrangement
        const Point_2 &a = halfEdge->source()->point();
        const Point_2 &b = halfEdge->target()->point();


        std::set<int> minusTrianglesSet;
        std::set<int> plusTrianglesSet;

        for (const auto &segmentId : segmentIds)
        {
            const std::array<int, 2> edge = tetMesh.edges[segmentId];
            const int segmentSourceId = edge[0];
            const Point_2 &c = singularArrangement.arrangementPoints[segmentSourceId];

            std::vector<int> minusTriangles;
            std::vector<int> plusTriangles;

            //
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
                minusTriangles = tetMesh.upperStarTriangles.at(edge);
                plusTriangles = tetMesh.lowerStarTriangles.at(edge);
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
            else if (CGAL::orientation(a, b, c) == CGAL::LEFT_TURN)
            {
                minusTriangles = tetMesh.lowerStarTriangles.at(edge);
                plusTriangles = tetMesh.upperStarTriangles.at(edge);
            }

            // Non-robust predicate issue
            else
            {
                assert(false);
            }


            minusTrianglesSet.insert(minusTriangles.begin(), minusTriangles.end());
            plusTrianglesSet.insert(plusTriangles.begin(), plusTriangles.end());
        }

        // Cancel out the plus/minus triangles and write to a vector
        std::set_difference(
                minusTrianglesSet.begin(), minusTrianglesSet.end(),
                plusTrianglesSet.begin(), plusTrianglesSet.end(),
                std::back_inserter(edgeRegionMinusTriangles[halfEdge]));

        std::set_difference(
                plusTrianglesSet.begin(), plusTrianglesSet.end(),
                minusTrianglesSet.begin(), minusTrianglesSet.end(),
                std::back_inserter(edgeRegionPlusTriangles[halfEdge]));

        edgeRegionMinusTriangles[halfEdge->twin()] = edgeRegionPlusTriangles[halfEdge];
        edgeRegionPlusTriangles[halfEdge->twin()] = edgeRegionMinusTriangles[halfEdge];



        //auto &minusVec = edgeRegionMinusTriangles[halfEdge];
        //minusVec.insert(minusVec.end(), minusTrianglesSet.begin(), minusTrianglesSet.end());

        //auto &plusVec = edgeRegionPlusTriangles[halfEdge];
        //plusVec.insert(plusVec.end(), plusTrianglesSet.begin(), plusTrianglesSet.end());

        //// These are reversed for the twin edge
        //edgeRegionMinusTriangles[halfEdge->twin()] = plusVec;
        //edgeRegionPlusTriangles[halfEdge->twin()] = minusVec;
    }
}


void ReebSpace2::computeEdgeCrossingMinusPlusTriangles(const TetMesh &tetMesh, Arrangement &singularArrangement)
{
    for (auto he = singularArrangement.arr.halfedges_begin(); he != singularArrangement.arr.halfedges_end(); ++he)
    {
        // If we have computed this for the twin, just swap them around
        if (edgeCrossingMinusTriangles.contains(he->twin()))
        {
            edgeCrossingMinusTriangles[he] = edgeCrossingPlusTriangles[he->twin()];
            edgeCrossingPlusTriangles[he] = edgeCrossingMinusTriangles[he->twin()];
        }
        else
        {
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

            // The half edge has the same direction as the original edge
            if (isSegmentLeftToRight == isCurrentHalfEdgeLeftToRight)
            {
                edgeCrossingMinusTriangles[he] = tetMesh.upperStarTriangles.at(edge);
                edgeCrossingPlusTriangles[he] = tetMesh.lowerStarTriangles.at(edge);
            }
            else
            {
                edgeCrossingMinusTriangles[he] = tetMesh.lowerStarTriangles.at(edge);
                edgeCrossingPlusTriangles[he] = tetMesh.upperStarTriangles.at(edge);
            }
        }
    }
}

void ReebSpace2::computeVertexRegionMinusPlusTriangles(const TetMesh &tetMesh, Arrangement &singularArrangement)
{
    for (const auto [halfEdge, segmentIds] : vertexRegionSegments)
    {
        // Each vertex is guaratneed to be in the singular arrangement
        const Point_2 &vertex = halfEdge->target()->point();
        const int vertexMeshId = singularArrangement.arrangementPointIndices[vertex];

        std::set<int> minusTrianglesSet;
        std::set<int> plusTrianglesSet;

        for (const auto &segmentId : segmentIds)
        {
            const std::array<int, 2> edge = tetMesh.edges[segmentId];

            std::vector<int> minusTriangles;
            std::vector<int> plusTriangles;

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
                minusTriangles = tetMesh.upperStarTriangles.at(edge);
                plusTriangles = tetMesh.lowerStarTriangles.at(edge);
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
            else if (edge[1] == vertexMeshId)
            {
                minusTriangles = tetMesh.lowerStarTriangles.at(edge);
                plusTriangles = tetMesh.upperStarTriangles.at(edge);
            }
            else
            {
                assert(false);
            }

            minusTrianglesSet.insert(minusTriangles.begin(), minusTriangles.end());
            plusTrianglesSet.insert(plusTriangles.begin(), plusTriangles.end());
        }

        // Cancel out the plus/minus triangles and write to a vector
        std::set_difference(
                minusTrianglesSet.begin(), minusTrianglesSet.end(),
                plusTrianglesSet.begin(), plusTrianglesSet.end(),
                std::back_inserter(vertexRegionMinusTriangles[halfEdge]));

        std::set_difference(
                plusTrianglesSet.begin(), plusTrianglesSet.end(),
                minusTrianglesSet.begin(), minusTrianglesSet.end(),
                std::back_inserter(vertexRegionPlusTriangles[halfEdge]));
    }
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

                halfEdgeEdgeRegionSegmentsUnitTest[he].insert(myRegularSegment.originatingEdge);
            }
        }
    }

    assert(areHalfEdgeRegionMapsEqual(edgeRegionSegments, halfEdgeEdgeRegionSegmentsUnitTest));



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


    assert(areHalfEdgeRegionMapsEqual(vertexRegionSegments, halfEdgeVertexRegionSegmentsUnitTest));
}

void ReebSpace2::checkInitialAssumptions(const TetMesh &tetMesh, Arrangement &singularArrangement)
{
    for (Face_const_iterator fit = singularArrangement.arr.faces_begin(); fit != singularArrangement.arr.faces_end(); ++fit) 
    {
        if (fit->is_unbounded()) 
        {
            if (
                    fit->number_of_inner_ccbs() != 1 ||
                    fit->number_of_outer_ccbs() != 0 ||
                    fit->number_of_holes() != 1 ||
                    fit->number_of_isolated_vertices() != 0
                    )
            {
                std::cerr << "number_of_inner_ccbs: " << fit->number_of_inner_ccbs() << std::endl;
                std::cerr << "number_of_outer_ccbs: " << fit->number_of_outer_ccbs() << std::endl;
                std::cerr << "number_of_holes: " << fit->number_of_holes() << std::endl;
                std::cerr << "number_of_isolated_vertices: " << fit->number_of_isolated_vertices() << std::endl;
                throw std::runtime_error("Outer face is degenerate!");
            }
        }
        else
        {
            if (
                    fit->number_of_inner_ccbs() != 0 ||
                    fit->number_of_outer_ccbs() != 1 ||
                    fit->number_of_holes() != 0 ||
                    fit->number_of_isolated_vertices() != 0
                    )
            {
                std::cerr << "number_of_inner_ccbs: " << fit->number_of_inner_ccbs() << std::endl;
                std::cerr << "number_of_outer_ccbs: " << fit->number_of_outer_ccbs() << std::endl;
                std::cerr << "number_of_holes: " << fit->number_of_holes() << std::endl;
                std::cerr << "number_of_isolated_vertices: " << fit->number_of_isolated_vertices() << std::endl;
                throw std::runtime_error("An inner face is degenerate!");
            }
        }
    }
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



void ReebSpace2::unitTestComparePreimageGraphs(const TetMesh &tetMesh, Arrangement &singularArrangement, Arrangement &regularArrangement, ReebSpace &rs)
{
    for (const auto& [halfEdge, preimageGraphs] : this->preimageGraphsCached)
    {
        std::pair<Halfedge_const_handle, Halfedge_const_handle> foundHalfEdges =  findHalfEdges(halfEdge, singularArrangement, regularArrangement);
        std::pair<Face_const_handle, Face_const_handle> foundFaces =  {foundHalfEdges.first->face(), foundHalfEdges.second->face()};
        std::pair<int, int> foundFaceIds =  {regularArrangement.arrangementFacesIdices.at(foundFaces.first), regularArrangement.arrangementFacesIdices.at(foundFaces.second)};

        std::pair<PreimageGraph, PreimageGraph> &singularPreimageGraphs = this->preimageGraphsCached[halfEdge];
        std::pair<DisjointSet<int>, DisjointSet<int>> regularPreimageGraphs = {rs.preimageGraphs[foundFaceIds.first], rs.preimageGraphs[foundFaceIds.second]};

        //assert(singularPreimageGraphs.first == regularPreimageGraphs.first);
        //assert(singularPreimageGraphs.second == regularPreimageGraphs.second);


        //if (false == singularPreimageGraphs.first.areEqual(regularPreimageGraphs.first) || false == singularPreimageGraphs.second.areEqual(regularPreimageGraphs.second))
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

            throw std::runtime_error("Preimage graphs do not match!");
        }
    }
}
