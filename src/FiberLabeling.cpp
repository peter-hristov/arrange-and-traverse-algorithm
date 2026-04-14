#include "./FiberLabeling.h"
#include "src/ReebSpace2.h"


std::vector<std::pair<int, int>> fiber::labeling::computeFiberSeeds(TetMesh &tetMesh, Arrangement &singularArrangement, ReebSpace2 &reebSpace, const std::array<double, 2> &controlPoint, const std::set<int> &selectedSheetIds)
{
    const Point_2 controlPointEPEC(controlPoint[0], controlPoint[1]);

    // 1. Compute the active face
    Face_const_handle activeFace = singularArrangement.getActiveFace(controlPoint);
    if (activeFace->is_unbounded()) { return {}; }

    const int activeFaceId = singularArrangement.arrangementFacesIdices[activeFace];

    // Pick the midpoint of the first half-edge, otherwise we will search through too many segments in the AABB tree
    const Point_2 endPoint = CGAL::midpoint(activeFace->outer_ccb()->source()->point(), activeFace->outer_ccb()->target()->point());

    //const Point_2 endPoint(controlPoint[0], controlPoint[1] + tetMesh.maxG + 10.0);

    const Segment_2 controlSegment(controlPointEPEC, endPoint);
    if (controlSegment.squared_length() == 0.0) { std::cerr << "The segment has zero lenght!" << std::endl; }

    int graphUpdates = 0;

    //for (const auto &p : singularArrangement.arrangementPoints)
    //{
        //if (endPoint == p)
        //{
            //std::cerr << "THE ENDPONT IS A VERTEX!";
        //}
    //}


    //Timer::start();
    std::vector<std::tuple<K::FT, int, int>> intersectedSegments = singularArrangement.getIntersectedSegments2(tetMesh, controlSegment, true, true);
    //Timer::stop("AABB search 2                          :");



    //for (const auto& [alpha, edgeId, edgeType] : intersectedSegments)
    //{
        //Point_2 p = controlSegment.source() + alpha * (controlSegment.target() - controlSegment.source());
        //std::cout << "alpha: " << CGAL::to_double(alpha) << " edge id: " << edgeId << " edge type : " << edgeType  << " and point " << p << std::endl;
    //}

    const int destinationSegmentId = std::get<1>(intersectedSegments.back());

    //std::cout << std::endl << "Destination segment is " << singularArrangement.arrangementPoints.at(tetMesh.edges[destinationSegmentId][0]) << " -> " << singularArrangement.arrangementPoints.at(tetMesh.edges[destinationSegmentId][1]) << " and ID = " << destinationSegmentId << std::endl << std::endl;


    //FiberGraph pg = reebSpace.representativeFiberGraphs[activeFace->data()];

    std::vector<std::pair<int, int>> fiberSeeds;


    // If there are selected sheets, only use their fiber components
    //
    for (const auto &[triangleId, componentId] : reebSpace.representativeFiberGraphSeeds[activeFace->data()])
    {
        if (!selectedSheetIds.empty())
        {
            const int sheetId = reebSpace.correspondenceGraphDS.find(componentId);
            if (!selectedSheetIds.contains(sheetId))
                continue;
        }

        fiberSeeds.push_back({triangleId, componentId});
    }



    // Construct a new fiber graph for that point

    //printf("This is the initial fiber graph:\n");
    //pg.printByRoot();

    // Go up to closestHalfEdgeVertexPoint
    //
    auto currentHalfEdge = activeFace->outer_ccb();
    do
    {
        // Get the segment ID of the half-edge
        //
        const Segment_2 &segment = *singularArrangement.arr.originating_curves_begin(currentHalfEdge);
        const int aIndex = singularArrangement.arrangementPointIndices.at(segment.source());
        const int bIndex = singularArrangement.arrangementPointIndices.at(segment.target());
        const std::array<int, 2> edge = {aIndex, bIndex};
        const int segmentId = tetMesh.edgeIndices.at(edge);

        //std::cout << "\nSegment with ID " << segmentId << " between " << segment.source() << " and " << segment.target() << " with desired ID " << destinationSegmentId << "\n";
        //std::cout << "The half-edge is " << currentHalfEdge->source()->point() << ", " << currentHalfEdge->target()->point() << "\n";


        if (segmentId == destinationSegmentId)
        {
            Segment_2 halfEdgeSegment(currentHalfEdge->source()->point(), currentHalfEdge->target()->point());
            if (CGAL::do_intersect(halfEdgeSegment, controlSegment))
            {
                break;
            }
        }

        //pg.updateComponentsRegular(tetMesh, reebSpace.edgeRegionSegments[currentHalfEdge->data().id]);
        FiberGraph::updateSeedsRegular(tetMesh, reebSpace.edgeRegionSegments[currentHalfEdge->data().id], fiberSeeds);
        graphUpdates++;

        //pg.updateComponentsRegular(tetMesh, reebSpace.vertexRegionSegments[currentHalfEdge->data().id]);
        FiberGraph::updateSeedsRegular(tetMesh, reebSpace.vertexRegionSegments[currentHalfEdge->data().id], fiberSeeds);
        graphUpdates++;
        ++currentHalfEdge;

    } while (true);

    //printf("\n\nFinal fiber graph after finding the desired segment:\n");
    //pg.printByRoot();
    //std::cout << "\n\n";


    //
    //
    //                          s.source()
    //                              |
    //                              |
    //                              |
    // searchSegment.source() ------x---------- searchSegment.target()
    //                              |
    //                              |
    //                              |
    //                              |
    //                          s.target()
    //
    const K::FT desiredAlpha = CGAL::Intersections::internal::s2s2_alpha(
            currentHalfEdge->target()->point().x(), 
            currentHalfEdge->target()->point().y(),
            currentHalfEdge->source()->point().x(), 
            currentHalfEdge->source()->point().y(),
            controlSegment.target().x(), 
            controlSegment.target().y(),
            controlSegment.source().x(), 
            controlSegment.source().y());

    //const double desiredAlpha = CGAL::Intersections::internal::s2s2_alpha(
            //CGAL::to_double(currentHalfEdge->target()->point().x()), 
            //CGAL::to_double(currentHalfEdge->target()->point().y()),
            //CGAL::to_double(currentHalfEdge->source()->point().x()), 
            //CGAL::to_double(currentHalfEdge->source()->point().y()),
            //CGAL::to_double(controlSegment.target().x()), 
            //CGAL::to_double(controlSegment.target().y()),
            //CGAL::to_double(controlSegment.source().x()), 
            //CGAL::to_double(controlSegment.source().y()));


    for (const std::pair<int, bool> &regionSegment :  reebSpace.edgeRegionSegments[currentHalfEdge->data().id])
    {
        const int edgeId = regionSegment.first;
        const std::array<int, 2> &edge = tetMesh.edges[edgeId];
        Segment_2 regularSegment(singularArrangement.arrangementPoints[edge[0]], singularArrangement.arrangementPoints[edge[1]]);



        K::FT alpha = CGAL::Intersections::internal::s2s2_alpha(
                currentHalfEdge->target()->point().x(), currentHalfEdge->target()->point().y(),
                currentHalfEdge->source()->point().x(), currentHalfEdge->source()->point().y(),
                regularSegment.target().x(), regularSegment.target().y(),
                regularSegment.source().x(), regularSegment.source().y());

        //const double alpha = CGAL::Intersections::internal::s2s2_alpha(
                //CGAL::to_double(currentHalfEdge->target()->point().x()), 
                //CGAL::to_double(currentHalfEdge->target()->point().y()),
                //CGAL::to_double(currentHalfEdge->source()->point().x()), 
                //CGAL::to_double(currentHalfEdge->source()->point().y()),
                //CGAL::to_double(regularSegment.target().x()), 
                //CGAL::to_double(regularSegment.target().y()),
                //CGAL::to_double(regularSegment.source().x()), 
                //CGAL::to_double(regularSegment.source().y()));


        //std::cout << "\nSegment with ID " << edgeId << " between " << regularSegment.source() << " and " << regularSegment.target() << "\n";
        //std::cout << "Reg segment alpha : " << alpha << " /  and desired alpha : " << desiredAlpha << std::endl << std::endl;

        if (alpha == desiredAlpha)
        {
            std::cerr << "DEGENERATE CASE DETECTED!!!!!!!!!!!!!";
            throw std::runtime_error("DEGEN CASE DETECTED!!!!!!!!!!!!!.");
        }

        if (alpha > desiredAlpha)
        {
            break;
        }
        
        //printf("Here the current fiber graph :\n");
        //pg.printByRoot();

        //pg.updateComponentsRegular(tetMesh, {regionSegment});
        FiberGraph::updateSeedsRegular(tetMesh, {regionSegment}, fiberSeeds);

        graphUpdates++;
    }

    //printf("Here the final fiber graph after fiding the intersection point :\n");
    //pg.printByRoot();
    //printf("\n\n\n");


    // The last one is the singular one, we don't want to cross it
    //
    for (int i = intersectedSegments.size() - 2 ; i >= 0 ; i--)
    {
        const auto [alpha, segmentId, type] = intersectedSegments[i];

        if (alpha == 1.0)
        {
            continue;
        }


        bool typicalOrientation = true;

        //std::cerr << "Intersected segment with ID " << segmentId << " and type " << tetMesh.edgeSingularTypes.at(tetMesh.edges.at(segmentId)) << " and alpha " << alpha << std::endl;

        // Change orientation in case we need to
        const std::array<int, 2> edge = tetMesh.edges.at(segmentId);
        const int segmentSourceId = edge[0];
        const int segmentTargetId = edge[1];

        const Point_2 &c = singularArrangement.arrangementPoints[segmentSourceId];
        const Point_2 &d = singularArrangement.arrangementPoints[segmentTargetId];

        //std::cout << c << " - " << d << std::endl;

        const Point_2 &a = controlPointEPEC;
        const Point_2 &b = endPoint;
        //
        // If the orientation this way around, reverse it
        //   (regular segment)
        //          d 
        //           \
        //            \
        // a ----------\--------- b (singular segment)
        //              \
        //               \
        //                c
        //
        if (CGAL::orientation(a, b, c) == CGAL::LEFT_TURN)
        {
            typicalOrientation = false;
        }


        //pg.printByRoot();


        //const std::vector<int> &minusTriangles = tetMesh.getMinusTriangles(segmentId, typicalOrientation);
        //const std::vector<int> &plusTriangles = tetMesh.getPlusTriangles(segmentId, typicalOrientation);

        //std::cout << "\n\n\n\nMinus triangles: " << std::endl;

        //for (const int &triangleId : minusTriangles)
        //{
            //std::cout << triangleId << std::endl;
        //}

        //std::cout << "\nPlus triangles: " << std::endl;

        //for (const int &triangleId : plusTriangles)
        //{
            //std::cout << triangleId << std::endl;
        //}


        //pg.updateComponentsRegular(tetMesh, {{segmentId, typicalOrientation}});
        FiberGraph::updateSeedsRegular(tetMesh, {{segmentId, typicalOrientation}}, fiberSeeds);
        graphUpdates++;

    }


    //std::cout << "We have performed " << graphUpdates << " graph updates.\n";


    return fiberSeeds;
}



std::vector<std::pair<int, int>> fiber::labeling::computeFiberSeedsGivenLine(TetMesh &tetMesh, Arrangement &singularArrangement, ReebSpace2 &reebSpace, const Segment_2 &controlSegment, const K::FT &pointAlpha, const std::vector<std::tuple<K::FT, int, int>> &intersectedSegments)
{

    const Point_2 controlPoint = CGAL::barycenter(controlSegment[0], 1.0 - pointAlpha, controlSegment[1], pointAlpha);
    if (controlSegment.squared_length() == 0.0) { std::cerr << "The segment has zero lenght!" << std::endl; }

    //std::cout << std::endl;
    //std::cout << std::endl;
    //std::cout << "------------------------------------------------------------------------------------";
    //std::cout << "Computing new flexible fiber.\n";
    //std::cout << "------------------------------------------------------------------------------------";
    //std::cout << std::endl;
    //std::cout << "The control segment is : " << controlSegment[0] << " -> " << controlSegment[1] << std::endl;
    //std::cout << "The input alpha   is   : "  << pointAlpha << std::endl;
    //std::cout << "The control piont is   : "  << controlPoint << std::endl;


    if (false == CGAL::do_intersect(controlPoint, controlSegment))
    {
        throw std::runtime_error("The control point is not on the control segment.");
    }

    // 1. Compute the active face
    Face_const_handle activeFace = singularArrangement.getActiveFace(controlPoint);
    if (activeFace->is_unbounded()) { return {}; }

    const int activeFaceId = singularArrangement.arrangementFacesIdices[activeFace];

    int graphUpdates = 0;

    //for (const auto &p : singularArrangement.arrangementPoints)
    //{
        //if (endPoint == p)
        //{
            //std::cerr << "THE ENDPONT IS A VERTEX!";
        //}
    //}

    //std::cout << "\nThese are all the intersected segments:\n";
    //for (const auto& [alpha, edgeId, edgeType] : intersectedSegments)
    //{
        //Point_2 p = controlSegment.source() + alpha * (controlSegment.target() - controlSegment.source());
        //std::cout << "alpha: " << CGAL::to_double(alpha) << " edge id: " << edgeId << " edge type : " << edgeType  << " and point " << p << std::endl;
    //}


    int destinationSegmentId = -1;
    int destinationIntersectedSegmentsId = -1;

    for (int i = 0 ; i < intersectedSegments.size() ; i++)
    {
        const auto &[alpha, edgeId, edgeType] = intersectedSegments[i];

        if (edgeType == 2 || edgeType == 0)
        {
            if (pointAlpha < alpha)
            {
                destinationSegmentId = edgeId;
                destinationIntersectedSegmentsId = i;
                break;
            }
        }
    }

    // The control segment does not intersect the boundary of this face, then just compute the fiber another way
    if (destinationSegmentId == -1)
    {
        std::array<double, 2> controlPointDouble = {CGAL::to_double(controlPoint.x()), CGAL::to_double(controlPoint.y())};
        return fiber::labeling::computeFiberSeeds(tetMesh, singularArrangement, reebSpace, controlPointDouble, {});
    }

    //const int destinationSegmentId = std::get<1>(intersectedSegments.back());



    //std::cout << std::endl << "Destination segment is " << singularArrangement.arrangementPoints.at(tetMesh.edges[destinationSegmentId][0]) << " -> " << singularArrangement.arrangementPoints.at(tetMesh.edges[destinationSegmentId][1]) << " and ID = " << destinationSegmentId << std::endl << std::endl;


    //FiberGraph pg = reebSpace.representativeFiberGraphs[activeFace->data()];
    std::vector<std::pair<int, int>> fiberSeeds = reebSpace.representativeFiberGraphSeeds[activeFace->data()];

    //printf("This is the initial fiber graph:\n");
    //pg.printByRoot();

    // Go up to closestHalfEdgeVertexPoint
    //
    auto currentHalfEdge = activeFace->outer_ccb();
    do
    {
        // Get the segment ID of the half-edge
        //
        const Segment_2 &segment = *singularArrangement.arr.originating_curves_begin(currentHalfEdge);
        const int aIndex = singularArrangement.arrangementPointIndices.at(segment.source());
        const int bIndex = singularArrangement.arrangementPointIndices.at(segment.target());
        const std::array<int, 2> edge = {aIndex, bIndex};
        const int segmentId = tetMesh.edgeIndices.at(edge);

        //std::cout << "\nSegment with ID " << segmentId << " between " << segment.source() << " and " << segment.target() << " with desired ID " << destinationSegmentId << "\n";
        //std::cout << "The half-edge is " << currentHalfEdge->source()->point() << ", " << currentHalfEdge->target()->point() << "\n";


        if (segmentId == destinationSegmentId)
        {
            Segment_2 halfEdgeSegment(currentHalfEdge->source()->point(), currentHalfEdge->target()->point());
            if (CGAL::do_intersect(halfEdgeSegment, controlSegment))
            {
                break;
            }
        }

        //pg.updateComponentsRegular(tetMesh, reebSpace.edgeRegionSegments[currentHalfEdge->data().id]);
        FiberGraph::updateSeedsRegular(tetMesh, reebSpace.edgeRegionSegments[currentHalfEdge->data().id], fiberSeeds);
        graphUpdates++;

        //pg.updateComponentsRegular(tetMesh, reebSpace.vertexRegionSegments[currentHalfEdge->data().id]);
        FiberGraph::updateSeedsRegular(tetMesh, reebSpace.vertexRegionSegments[currentHalfEdge->data().id], fiberSeeds);
        graphUpdates++;

        ++currentHalfEdge;

        // If we are back at the start but we have not foudn the destination segment
        if (currentHalfEdge == activeFace->outer_ccb())
        {
            throw std::runtime_error("Desired half-edge not found!");
        }

    } while (true);

    //printf("\n\nFinal fiber graph after finding the desired segment:\n");
    //pg.printByRoot();
    //std::cout << "\n\n";


    //
    //
    //                          s.source()
    //                              |
    //                              |
    //                              |
    // searchSegment.source() ------x---------- searchSegment.target()
    //                              |
    //                              |
    //                              |
    //                              |
    //                          s.target()
    //
    const K::FT desiredAlpha = CGAL::Intersections::internal::s2s2_alpha(
            currentHalfEdge->target()->point().x(), 
            currentHalfEdge->target()->point().y(),
            currentHalfEdge->source()->point().x(), 
            currentHalfEdge->source()->point().y(),
            controlSegment.target().x(), 
            controlSegment.target().y(),
            controlSegment.source().x(), 
            controlSegment.source().y());

    //const double desiredAlpha = CGAL::Intersections::internal::s2s2_alpha(
            //CGAL::to_double(currentHalfEdge->target()->point().x()), 
            //CGAL::to_double(currentHalfEdge->target()->point().y()),
            //CGAL::to_double(currentHalfEdge->source()->point().x()), 
            //CGAL::to_double(currentHalfEdge->source()->point().y()),
            //CGAL::to_double(controlSegment.target().x()), 
            //CGAL::to_double(controlSegment.target().y()),
            //CGAL::to_double(controlSegment.source().x()), 
            //CGAL::to_double(controlSegment.source().y()));


    for (const std::pair<int, bool> &regionSegment :  reebSpace.edgeRegionSegments[currentHalfEdge->data().id])
    {
        const int edgeId = regionSegment.first;
        const std::array<int, 2> &edge = tetMesh.edges[edgeId];
        Segment_2 regularSegment(singularArrangement.arrangementPoints[edge[0]], singularArrangement.arrangementPoints[edge[1]]);



        K::FT alpha = CGAL::Intersections::internal::s2s2_alpha(
                currentHalfEdge->target()->point().x(), currentHalfEdge->target()->point().y(),
                currentHalfEdge->source()->point().x(), currentHalfEdge->source()->point().y(),
                regularSegment.target().x(), regularSegment.target().y(),
                regularSegment.source().x(), regularSegment.source().y());

        //const double alpha = CGAL::Intersections::internal::s2s2_alpha(
                //CGAL::to_double(currentHalfEdge->target()->point().x()), 
                //CGAL::to_double(currentHalfEdge->target()->point().y()),
                //CGAL::to_double(currentHalfEdge->source()->point().x()), 
                //CGAL::to_double(currentHalfEdge->source()->point().y()),
                //CGAL::to_double(regularSegment.target().x()), 
                //CGAL::to_double(regularSegment.target().y()),
                //CGAL::to_double(regularSegment.source().x()), 
                //CGAL::to_double(regularSegment.source().y()));


        //std::cout << "\nSegment with ID " << edgeId << " between " << regularSegment.source() << " and " << regularSegment.target() << "\n";
        //std::cout << "Reg segment alpha : " << alpha << " /  and desired alpha : " << desiredAlpha << std::endl << std::endl;

        if (alpha == desiredAlpha)
        {
            std::cerr << "DEGENERATE CASE DETECTED!!!!!!!!!!!!!";
            throw std::runtime_error("DEGEN CASE DETECTED!!!!!!!!!!!!!.");
        }

        if (alpha > desiredAlpha)
        {
            break;
        }
        
        //printf("Here the current fiber graph :\n");
        //pg.printByRoot();

        //pg.updateComponentsRegular(tetMesh, {regionSegment});
        FiberGraph::updateSeedsRegular(tetMesh, {regionSegment}, fiberSeeds);
        graphUpdates++;
    }

    //printf("Here the final fiber graph after fiding the intersection point :\n");
    //pg.printByRoot();
    //printf("\n\n\n\n\n\n\n");


    // The last one is the singular one, we don't want to cross it
    //
    for (int i = destinationIntersectedSegmentsId - 1 ; i >= 0 ; i--)
    {
        const auto [alpha, segmentId, type] = intersectedSegments[i];

        if (alpha == 1.0)
        {
            continue;
        }

        // No need to go futher
        if (alpha < pointAlpha)
        {
            break;
        }

        bool typicalOrientation = true;

        //std::cerr << "Intersected segment with ID " << segmentId << " and type " << tetMesh.edgeSingularTypes.at(tetMesh.edges.at(segmentId)) << " and alpha " << alpha << " and point alpha " << pointAlpha << std::endl;

        // Change orientation in case we need to
        const std::array<int, 2> edge = tetMesh.edges.at(segmentId);
        const int segmentSourceId = edge[0];
        const int segmentTargetId = edge[1];

        const Point_2 &c = singularArrangement.arrangementPoints[segmentSourceId];
        const Point_2 &d = singularArrangement.arrangementPoints[segmentTargetId];

        //std::cout << c << " - " << d << std::endl;

        const Point_2 &a = controlPoint;
        const Point_2 &b = controlSegment[1];
        //
        // If the orientation this way around, reverse it
        //   (regular segment)
        //          d 
        //           \
        //            \
        // a ----------\--------- b (singular segment)
        //              \
        //               \
        //                c
        //
        if (CGAL::orientation(a, b, c) == CGAL::LEFT_TURN)
        {
            typicalOrientation = false;
        }




        //const std::vector<int> &minusTriangles = tetMesh.getMinusTriangles(segmentId, typicalOrientation);
        //const std::vector<int> &plusTriangles = tetMesh.getPlusTriangles(segmentId, typicalOrientation);

        //pg.printByRoot();

        //std::cout << "\n\nMinus triangles: " << std::endl;

        //for (const int &triangleId : minusTriangles)
        //{
            //std::cout << triangleId << std::endl;
        //}

        //std::cout << "\nPlus triangles: " << std::endl;

        //for (const int &triangleId : plusTriangles)
        //{
            //std::cout << triangleId << std::endl;
        //}


        //pg.updateComponentsRegular(tetMesh, {{segmentId, typicalOrientation}});
        FiberGraph::updateSeedsRegular(tetMesh, {{segmentId, typicalOrientation}}, fiberSeeds);
        graphUpdates++;

    }



    return fiberSeeds;
}

FiberGraph fiber::labeling::computeFiberGraph(TetMesh &tetMesh, Arrangement &singularArrangement, ReebSpace2 &reebSpace, std::array<double, 2> controlPoint)
{
    const Point_2 controlPointEPEC(controlPoint[0], controlPoint[1]);

    // 1. Compute the active face
    Face_const_handle activeFace = singularArrangement.getActiveFace(controlPoint);
    if (activeFace->is_unbounded()) { return {}; }

    const int activeFaceId = singularArrangement.arrangementFacesIdices[activeFace];

    // Pick the midpoint of the first half-edge, otherwise we will search through too many segments in the AABB tree
    const Point_2 endPoint = CGAL::midpoint(activeFace->outer_ccb()->source()->point(), activeFace->outer_ccb()->target()->point());

    //const Point_2 endPoint(controlPoint[0], controlPoint[1] + tetMesh.maxG + 10.0);

    const Segment_2 controlSegment(controlPointEPEC, endPoint);
    if (controlSegment.squared_length() == 0.0) { std::cerr << "The segment has zero lenght!" << std::endl; }

    int graphUpdates = 0;

    //for (const auto &p : singularArrangement.arrangementPoints)
    //{
        //if (endPoint == p)
        //{
            //std::cerr << "THE ENDPONT IS A VERTEX!";
        //}
    //}


    //Timer::start();
    std::vector<std::tuple<K::FT, int, int>> intersectedSegments = singularArrangement.getIntersectedSegments2(tetMesh, controlSegment, true, true);
    //Timer::stop("AABB search 2                          :");



    //for (const auto& [alpha, edgeId, edgeType] : intersectedSegments)
    //{
        //Point_2 p = controlSegment.source() + alpha * (controlSegment.target() - controlSegment.source());
        //std::cout << "alpha: " << CGAL::to_double(alpha) << " edge id: " << edgeId << " edge type : " << edgeType  << " and point " << p << std::endl;
    //}

    const int destinationSegmentId = std::get<1>(intersectedSegments.back());

    //std::cout << std::endl << "Destination segment is " << singularArrangement.arrangementPoints.at(tetMesh.edges[destinationSegmentId][0]) << " -> " << singularArrangement.arrangementPoints.at(tetMesh.edges[destinationSegmentId][1]) << " and ID = " << destinationSegmentId << std::endl << std::endl;



    FiberGraph pg = reebSpace.representativeFiberGraphs[activeFace->data()];

    //printf("This is the initial fiber graph:\n");
    //pg.printByRoot();

    // Go up to closestHalfEdgeVertexPoint
    //
    auto currentHalfEdge = activeFace->outer_ccb();
    do
    {
        // Get the segment ID of the half-edge
        //
        const Segment_2 &segment = *singularArrangement.arr.originating_curves_begin(currentHalfEdge);
        const int aIndex = singularArrangement.arrangementPointIndices.at(segment.source());
        const int bIndex = singularArrangement.arrangementPointIndices.at(segment.target());
        const std::array<int, 2> edge = {aIndex, bIndex};
        const int segmentId = tetMesh.edgeIndices.at(edge);

        //std::cout << "\nSegment with ID " << segmentId << " between " << segment.source() << " and " << segment.target() << " with desired ID " << destinationSegmentId << "\n";
        //std::cout << "The half-edge is " << currentHalfEdge->source()->point() << ", " << currentHalfEdge->target()->point() << "\n";


        if (segmentId == destinationSegmentId)
        {
            Segment_2 halfEdgeSegment(currentHalfEdge->source()->point(), currentHalfEdge->target()->point());
            if (CGAL::do_intersect(halfEdgeSegment, controlSegment))
            {
                break;
            }
        }

        pg.updateComponentsRegular(tetMesh, reebSpace.edgeRegionSegments[currentHalfEdge->data().id]);
        graphUpdates++;
        pg.updateComponentsRegular(tetMesh, reebSpace.vertexRegionSegments[currentHalfEdge->data().id]);
        graphUpdates++;
        ++currentHalfEdge;

    } while (true);

    //printf("\n\nFinal fiber graph after finding the desired segment:\n");
    //pg.printByRoot();
    //std::cout << "\n\n";


    //
    //
    //                          s.source()
    //                              |
    //                              |
    //                              |
    // searchSegment.source() ------x---------- searchSegment.target()
    //                              |
    //                              |
    //                              |
    //                              |
    //                          s.target()
    //
    const K::FT desiredAlpha = CGAL::Intersections::internal::s2s2_alpha(
            currentHalfEdge->target()->point().x(), 
            currentHalfEdge->target()->point().y(),
            currentHalfEdge->source()->point().x(), 
            currentHalfEdge->source()->point().y(),
            controlSegment.target().x(), 
            controlSegment.target().y(),
            controlSegment.source().x(), 
            controlSegment.source().y());

    //const double desiredAlpha = CGAL::Intersections::internal::s2s2_alpha(
            //CGAL::to_double(currentHalfEdge->target()->point().x()), 
            //CGAL::to_double(currentHalfEdge->target()->point().y()),
            //CGAL::to_double(currentHalfEdge->source()->point().x()), 
            //CGAL::to_double(currentHalfEdge->source()->point().y()),
            //CGAL::to_double(controlSegment.target().x()), 
            //CGAL::to_double(controlSegment.target().y()),
            //CGAL::to_double(controlSegment.source().x()), 
            //CGAL::to_double(controlSegment.source().y()));


    for (const std::pair<int, bool> &regionSegment :  reebSpace.edgeRegionSegments[currentHalfEdge->data().id])
    {
        const int edgeId = regionSegment.first;
        const std::array<int, 2> &edge = tetMesh.edges[edgeId];
        Segment_2 regularSegment(singularArrangement.arrangementPoints[edge[0]], singularArrangement.arrangementPoints[edge[1]]);



        K::FT alpha = CGAL::Intersections::internal::s2s2_alpha(
                currentHalfEdge->target()->point().x(), currentHalfEdge->target()->point().y(),
                currentHalfEdge->source()->point().x(), currentHalfEdge->source()->point().y(),
                regularSegment.target().x(), regularSegment.target().y(),
                regularSegment.source().x(), regularSegment.source().y());

        //const double alpha = CGAL::Intersections::internal::s2s2_alpha(
                //CGAL::to_double(currentHalfEdge->target()->point().x()), 
                //CGAL::to_double(currentHalfEdge->target()->point().y()),
                //CGAL::to_double(currentHalfEdge->source()->point().x()), 
                //CGAL::to_double(currentHalfEdge->source()->point().y()),
                //CGAL::to_double(regularSegment.target().x()), 
                //CGAL::to_double(regularSegment.target().y()),
                //CGAL::to_double(regularSegment.source().x()), 
                //CGAL::to_double(regularSegment.source().y()));


        //std::cout << "\nSegment with ID " << edgeId << " between " << regularSegment.source() << " and " << regularSegment.target() << "\n";
        //std::cout << "Reg segment alpha : " << alpha << " /  and desired alpha : " << desiredAlpha << std::endl << std::endl;

        if (alpha == desiredAlpha)
        {
            std::cerr << "DEGENERATE CASE DETECTED!!!!!!!!!!!!!";
            throw std::runtime_error("DEGEN CASE DETECTED!!!!!!!!!!!!!.");
        }

        if (alpha > desiredAlpha)
        {
            break;
        }
        
        //printf("Here the current fiber graph :\n");
        //pg.printByRoot();

        pg.updateComponentsRegular(tetMesh, {regionSegment});
        graphUpdates++;
    }

    //printf("Here the final fiber graph after fiding the intersection point :\n");
    //pg.printByRoot();
    //printf("\n\n\n");


    // The last one is the singular one, we don't want to cross it
    //
    for (int i = intersectedSegments.size() - 2 ; i >= 0 ; i--)
    {
        const auto [alpha, segmentId, type] = intersectedSegments[i];

        if (alpha == 1.0)
        {
            continue;
        }


        bool typicalOrientation = true;

        //std::cerr << "Intersected segment with ID " << segmentId << " and type " << tetMesh.edgeSingularTypes.at(tetMesh.edges.at(segmentId)) << " and alpha " << alpha << std::endl;

        // Change orientation in case we need to
        const std::array<int, 2> edge = tetMesh.edges.at(segmentId);
        const int segmentSourceId = edge[0];
        const int segmentTargetId = edge[1];

        const Point_2 &c = singularArrangement.arrangementPoints[segmentSourceId];
        const Point_2 &d = singularArrangement.arrangementPoints[segmentTargetId];

        //std::cout << c << " - " << d << std::endl;

        const Point_2 &a = controlPointEPEC;
        const Point_2 &b = endPoint;
        //
        // If the orientation this way around, reverse it
        //   (regular segment)
        //          d 
        //           \
        //            \
        // a ----------\--------- b (singular segment)
        //              \
        //               \
        //                c
        //
        if (CGAL::orientation(a, b, c) == CGAL::LEFT_TURN)
        {
            typicalOrientation = false;
        }


        //pg.printByRoot();


        const std::vector<int> &minusTriangles = tetMesh.getMinusTriangles(segmentId, typicalOrientation);
        const std::vector<int> &plusTriangles = tetMesh.getPlusTriangles(segmentId, typicalOrientation);

        //std::cout << "\n\n\n\nMinus triangles: " << std::endl;

        //for (const int &triangleId : minusTriangles)
        //{
            //std::cout << triangleId << std::endl;
        //}

        //std::cout << "\nPlus triangles: " << std::endl;

        //for (const int &triangleId : plusTriangles)
        //{
            //std::cout << triangleId << std::endl;
        //}


        pg.updateComponentsRegular(tetMesh, {{segmentId, typicalOrientation}});
        graphUpdates++;

    }


    //std::cout << "We have performed " << graphUpdates << " graph updates.\n";



    return pg;
}
