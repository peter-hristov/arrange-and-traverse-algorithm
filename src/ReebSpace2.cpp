#include <cassert>
#include <memory>
#include <unordered_set>
#include <utility>
#include <omp.h>

#include "./io.h"
#include "./Timer.h"
#include "./ReebSpace2.h"
#include "./DisjointSet.h"
#include "./LoadingBar.hpp"
#include "src/FiberGraph.h"


#include <CGAL/create_straight_skeleton_2.h>
#include <CGAL/create_offset_polygons_2.h>
#include <CGAL/Polygon_with_holes_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K_EPIC;

typedef CGAL::Polygon_2<K_EPIC>           Polygon_2_EPIC;
typedef CGAL::Straight_skeleton_2<K_EPIC> Ss;
typedef std::shared_ptr<Ss> SsPtr;



std::vector<std::vector<std::array<double, 2>>> ReebSpace2::computeSheetControlPolygons(const int sheetId)
{
    const std::vector<std::vector<Halfedge_const_handle>> &currentSheetBoundaries = this->sheetBoundaries.at(sheetId);

    std::vector<std::vector<std::array<double, 2>>> sheetControlPolygonsPoints;


    std::vector<bool> isBoundaryHole(currentSheetBoundaries.size(), false);
    std::vector<Polygon_2_EPIC> sheetBoundaryPolygons;

    // The sheet boundary is CCW, the holes are CW
    for (int i = 0 ; i <  currentSheetBoundaries.size() ; i++)
    {
        const auto &sheetBoundaryComponent = currentSheetBoundaries[i];

        Polygon_2_EPIC poly;

        for (int j = 0 ; j < sheetBoundaryComponent.size() ; j++)
        {
            const auto h1 = sheetBoundaryComponent[j];

            const double x = CGAL::to_double(h1->source()->point().x());
            const double y = CGAL::to_double(h1->source()->point().y());

            //poly.push_back({x, y});

            poly.push_back(K_EPIC::Point_2(x, y));
        }

        if (false == poly.is_counterclockwise_oriented())
        {
            poly.reverse_orientation();
            isBoundaryHole[i] = true;
        }

        sheetBoundaryPolygons.push_back(poly);
    }

    const double offset = 0.00001;

    for (int i = 0 ; i < sheetBoundaryPolygons.size() ; i++)
    {
        std::vector<std::shared_ptr<Polygon_2_EPIC>> offset_polygons;

        if (true == isBoundaryHole[i])
        {
            offset_polygons = CGAL::create_exterior_skeleton_and_offset_polygons_2(offset, sheetBoundaryPolygons[i]);
        } 
        else
        {
            offset_polygons = CGAL::create_interior_skeleton_and_offset_polygons_2(offset, sheetBoundaryPolygons[i]);
        }

        //if (offset_polygons.size() != 1)
        //{
            //std::cout << "New polygon has " << offset_polygons.size() << " components." << std::endl;
            //throw std::runtime_error("One sheet boundary polygon split into more polygons.");
        //}

        for (const auto &offset_polygon : offset_polygons)
        {
            std::vector<std::array<double, 2>> controlPolygon;

            for (auto v = offset_polygon->vertices_begin(); v != offset_polygon->vertices_end(); ++v)
            {
                const double x = CGAL::to_double(v->x());
                const double y = CGAL::to_double(v->y());

                controlPolygon.push_back(std::array<double, 2>{x, y});
            }

            sheetControlPolygonsPoints.push_back(controlPolygon);
        }

    }














    //for (int i = 0 ; i <  currentSheetBoundaries.size() ; i++)
    //for (int i = 0 ; i <  currentSheetBoundaries.size() ; i++)
    //int i = 0 ;
    //{
        //const auto &sheetBoundaryComponent = currentSheetBoundaries[i];
        //sheetControlPolygonsPoints[i].reserve(sheetBoundaryComponent.size());

        //Polygon_2_EPIC poly;

        //for (int j = 0 ; j < sheetBoundaryComponent.size() ; j++)
        //{
            //const auto h1 = sheetBoundaryComponent[j];

            //const double x = CGAL::to_double(h1->source()->point().x());
            //const double y = CGAL::to_double(h1->source()->point().y());

            //poly.push_back({x, y});
        //}

        //if (false == poly.is_counterclockwise_oriented())
        //{
            //throw std::runtime_error("Boundary Polygon is not clockwise oriented.");
        //}

        //SsPtr iss = CGAL::create_interior_straight_skeleton_2(poly.vertices_begin(), poly.vertices_end());

        //double epsilon = 0.01; // your offset distance
        //auto offset_polygons = CGAL::create_offset_polygons_2<Polygon_2_EPIC>(epsilon, *iss);
    //////auto offset_polygons = CGAL::create_offset_polygons_2<Polygon_with_holes_EPIC>(epsilon, *iss);
        //for (const auto& offset_poly : offset_polygons)
        //{
            //for (auto v = offset_poly->vertices_begin(); v != offset_poly->vertices_end(); ++v)
            //{
                //const double x = CGAL::to_double(v->x());
                //const double y = CGAL::to_double(v->y());

                //sheetControlPolygonsPoints[i].emplace_back(std::array<double, 2>{x, y});
            //}
        //}

    //}
        //sheetControlPolygonsPoints.push_back(std::move(pts)); // add as a separate polygon
    //}


    return sheetControlPolygonsPoints;
}



std::vector<Halfedge_const_handle> getBoundaryComponent(const std::unordered_set<Halfedge_const_handle> &sheetHalfEdges, const Halfedge_const_handle &startingHalfedge, std::unordered_set<Halfedge_const_handle> &visitedHalfEdges)
{
    std::vector<Halfedge_const_handle> sheetBoundaryComponent = {startingHalfedge};

    auto current = sheetBoundaryComponent[0];
    do
    {
        auto next = current->next();

        while (false == sheetHalfEdges.contains(next))
        {
            next = next->twin()->next();
        }



        // Detect pinch points
        //int counter = 0;
        //auto test = current;
        //do
        //{
            //if (sheetHalfEdges.contains(test))
            //{
                //counter++;
            //}

            //test = test->twin()->next();
        //} while (test != current);

        //if (counter != 2)
        //{
            ////std::cout << "------------------------------------- Pinch point detected! " << counter << "\n";
        //}




        const Point_2& p = current->target()->point();

        current = next;

        sheetBoundaryComponent.push_back(current);
        visitedHalfEdges.insert(current);

    } while (current->target() != sheetBoundaryComponent[0]->source());

    return sheetBoundaryComponent;
}



void ReebSpace2::computeSheetBoundaries(Arrangement &singularArrangement)
{

    // Collect the half-edge on the boundary of each sheet 
    //
    std::map<int, std::unordered_set<Halfedge_const_handle>> halfEdgePerSheet;

    for (auto he = singularArrangement.arr.halfedges_begin(); he != singularArrangement.arr.halfedges_end(); ++he)
    {
        if (true == he->face()->is_fictitious() || true == he->face()->is_unbounded())
        {
            continue;
        }

        // Collect face sheets
        //
        const int faceId = he->face()->data();
        std::set<int> faceSheets;
        for (const int &componentId : this->correspondenceGraph[faceId])
        {
            const int sheetId = this->correspondenceGraphDS.find(componentId);
            faceSheets.insert(sheetId);
        }

        // Collect twin face sheets
        //
        const int twinFaceId = he->twin()->face()->data();
        std::set<int> twinFaceSheets;
        for (const int &componentId : this->correspondenceGraph[twinFaceId])
        {
            const int sheetId = this->correspondenceGraphDS.find(componentId);
            twinFaceSheets.insert(sheetId);
        }


        // If this face has a sheet, that the twin does not, it's on the boundary
        //
        for (const int sheetId : faceSheets)
        {
            if (false == twinFaceSheets.contains(sheetId))
            {
                halfEdgePerSheet[sheetId].insert(he);
            }
        }
    }



    // Loop around the boundaries of each sheet to sort the half-edges
    //
    for (const auto &[sheetId, halfEdges] : halfEdgePerSheet)
    {
        const auto &sheetHalfEdges = halfEdgePerSheet.at(sheetId);

        std::unordered_set<Halfedge_const_handle> visitedHalfEdges;

        for (const auto &halfEdge : halfEdges)
        {
            if (false == visitedHalfEdges.contains(halfEdge))
            {
                this->sheetBoundaries[sheetId].push_back(getBoundaryComponent(sheetHalfEdges, halfEdge, visitedHalfEdges));
            }
        }
    }
}




void ReebSpace2::computeSheets(Arrangement &singularArrangement)
{
    double totalArea;

    for (auto face = singularArrangement.arr.faces_begin(); face != singularArrangement.arr.faces_end(); ++face) 
    {
        if (face->is_unbounded() || !face->has_outer_ccb())
        {
            continue;
        }

        const int faceId = face->data();
        double faceArea = 0.0;


        CartesianPolygon_2 poly;

        auto circ = face->outer_ccb();
        auto curr = circ;
        do {
            // Get point from CGAL (and convert to double )
            const double u = CGAL::to_double(curr->source()->point().x());
            const double v = CGAL::to_double(curr->source()->point().y());

            poly.push_back({u, v});
        } while (++curr != circ);


        faceArea += poly.area();
        totalArea += faceArea;


        // For each component in this face
        for (const int &componentId : this->correspondenceGraph[faceId])
        {
            const int sheetId = this->correspondenceGraphDS.find(componentId);
            sheetArea[sheetId] += faceArea;
        }
    }


    for (const auto &[sheetId, area] : this->sheetArea)
    {
        this->sheetAreaProportion[sheetId] = area /  totalArea;

    }

    // Print out the sheets in sorted order
    std::vector<std::pair<int, double>> sortedSheets;

    sortedSheets.reserve(sheetArea.size());
    for (const auto& [id, area] : sheetArea)
        sortedSheets.emplace_back(id, area);

    std::sort(sortedSheets.begin(), sortedSheets.end(),
            [](const auto& a, const auto& b)
            {
            return a.second > b.second; // ascending by area
            });

    for (int i = 0 ; i < std::min(sortedSheets.size(), (size_t)30) ; i++)
    {
        const auto &[sheetId, area] = sortedSheets[i];
        std::cout << i << ": sheet " << sheetId << " has area " << area << " which is a ratio of : " << 100.0 * this->sheetAreaProportion[sheetId] <<  std::endl;
    }


    for (int i = 0 ; i < sortedSheets.size() ; i++)
    {
        this->sheetOrder.push_back(sortedSheets[i].first);
    }

    // Can we get some fibers for each
    // How many non-empty fiber graphs are t here? should be the same as the number of faces


    int nonEmptyFiberGraphs = 0;
    for (const auto &[firstFG, secondFG] : this->fiberGraphs)
    {
        if (firstFG.componentRoot.size() != 0)
        {
            nonEmptyFiberGraphs++;
        }

    }

    // For every face face
    //for (int faceId = 0 ; faceId < this->fiberGraphsPerFace.size(); faceId++)
    //{

        //// For every sheet in the face
        //for (const int &componentId : this->correspondenceGraph[faceId])
        //{
            //const int sheetId = this->correspondenceGraphDS.find(componentId);


            //// For every thiangle in the preimage graph
            //for (const auto &[triangleId, triangleComponentId] : this->fiberGraphsPerFace[faceId].componentRoot)
            //{

                //// If this triangle is in the same component as the sheet
                //const int sheetIdComponent = this->correspondenceGraphDS.find(triangleComponentId);

                //if (sheetId == sheetIdComponent)
                //{
                    //this->trianglesPerSheet[sheetId].push_back(triangleId);
                //}

            //}
        //}
    //}


    //for (const auto &[sheetId, trianglesVector] : this->trianglesPerSheet)
    //{
        //printf("We are looking at sheet %d\n", sheetId);

        //for (const int triangleId : trianglesVector)
        //{
            //printf("%d, ", triangleId);
        //}
        //printf("\n");
    //}
    //printf("Number of non-empty fiber graphs %d and number of faces %d\n", nonEmptyFiberGraphs, singularArrangement.arr.number_of_faces());
}















void ReebSpace2::seedFace(TetMesh &tetMesh, const Halfedge_const_handle &currentHalfEdge)
{
    //printf("\n--------------------------------------------------------------------------------------\n");
    //std::cout << "Seeding face " << currentHalfEdge->twin()->face()->data();
    //printf("\n--------------------------------------------------------------------------------------\n");

    // We assume that the fiber graph of the face we have come from has been computed already
    FiberGraph &pg = this->fiberGraphs[currentHalfEdge->twin()->data().id].second;


    //printf("\n-----------------------------------------------------\n");
    //std::cout << "Half-edge is [" << currentHalfEdge->twin()->source()->point() << "] -> [" << currentHalfEdge->twin()->target()->point() << "]";
    //printf("\n-----------------------------------------------------\n");
    //std::cout << "Half-edge second (source) is: \n";
    //pg.printByRoot();

    // We wish to compute the first fiber graph for this face
    FiberGraph &pg2 = this->fiberGraphs[currentHalfEdge->data().id].first;
    pg2 = pg;

    if (currentHalfEdge->data().isPseudoSingular)
    {
        pg2.updateComponentsRegular(tetMesh, {this->edgeCrossingSegments[currentHalfEdge->twin()->data().id]});
    }
    else
    {

         const int componentsBefore = FiberGraph::componentCount;

        //std::cout << "Preimage graphs after crossing from " << currentHalfEdge->twin()->face()->data() << " to " << currentHalfEdge->twin()->face()->data() << std::endl;
        //std::cout << "Before...\n";
        //pg2.printByRoot();

        pg2.updateComponentsSingular(tetMesh, this->edgeCrossingSegments[currentHalfEdge->twin()->data().id]);

        //std::cout << "After...\n";
        //pg2.printByRoot();


        // Add the new components to the correspondence graph
        const int componentsAfter = FiberGraph::componentCount;
        this->correspondenceGraphDS.add(componentsAfter - componentsBefore);
    }


    //std::cout << "Half-edge first (seeded) is: \n";
    //pg2.printByRoot();
}




void ReebSpace2::loopFace(TetMesh &tetMesh, const Halfedge_const_handle &initialHalfEdge)
{
    //printf("\n--------------------------------------------------------------------------------------\n");
    //std::cout << "Looping face " << initialHalfEdge->face()->data();
    //printf("\n--------------------------------------------------------------------------------------\n");

    //printf("\n-----------------------------------------------------\n");
    //std::cout << "Half-edge is [" << initialHalfEdge->source()->point() << "] -> [" << initialHalfEdge->target()->point() << "]";
    //printf("\n-----------------------------------------------------\n");

    FiberGraph pg = this->fiberGraphs[initialHalfEdge->data().id].first;

    //std::cout << "First preimage graph is: \n";
    //pg.printByRoot();

    pg.updateComponentsRegular(tetMesh, this->edgeRegionSegments[initialHalfEdge->data().id]);
    this->fiberGraphs[initialHalfEdge->data().id].second = pg;

    //std::cout << "\nSecond preimage graph is: \n";
    //pg.printByRoot();

    Halfedge_const_handle currentHalfEdge = initialHalfEdge->next();
    do
    {
        pg.updateComponentsRegular(tetMesh, this->vertexRegionSegments[currentHalfEdge->prev()->data().id]);
        this->fiberGraphs[currentHalfEdge->data().id].first = pg;

        pg.updateComponentsRegular(tetMesh, this->edgeRegionSegments[currentHalfEdge->data().id]);
        this->fiberGraphs[currentHalfEdge->data().id].second = pg;


        //printf("\n-----------------------------------------------------\n");
        //std::cout << "Half-edge is [" << currentHalfEdge->source()->point() << "] -> [" << currentHalfEdge->target()->point() << "]";
        //printf("\n-----------------------------------------------------\n");
        //std::cout << "First preimage graph is: \n";
        //this->fiberGraphs[currentHalfEdge->data()].first.printByRoot();
        //std::cout << "\nSecond preimage graph is: \n";
        //this->fiberGraphs[currentHalfEdge->data()].second.printByRoot();


        currentHalfEdge = currentHalfEdge->next();

    } while (currentHalfEdge != initialHalfEdge);
}


void ReebSpace2::traverse(TetMesh &tetMesh, Arrangement &singularArrangement, const bool cachePreimageGraphs)
{
    // Set up arrays
    //
    this->correspondenceGraph.resize(singularArrangement.arr.number_of_faces());
    this->fiberGraphs.resize(singularArrangement.arr.number_of_halfedges());

    //this->representativeFiberGraphs.resize(singularArrangement.arr.number_of_faces());
    this->representativeFiberGraphSeeds.resize(singularArrangement.arr.number_of_faces());

    // Find the outside face
    //
    Halfedge_const_handle startingHalfedge;
    for (Face_const_iterator fit = singularArrangement.arr.faces_begin(); fit != singularArrangement.arr.faces_end(); ++fit) 
    {
        if (fit->is_unbounded()) 
        {
            startingHalfedge = *fit->holes_begin();
            break;
        }
    }

    // Queue for our BFS search over the faces of the arrangement
    //
    std::queue<Halfedge_const_handle> traversalQueue;
    traversalQueue.push(startingHalfedge);

    // Track how many faces we have processed so far
    //
    int orderIndex = 0;

    // Order in which singular facess are processed
    // This allows us to establish correspondence between faces in a structured way by tracking the front of the BFS
    // Also serves as a visited array in the BFS
    std::vector<int> order(singularArrangement.arr.number_of_faces(), -1);
    order[startingHalfedge->face()->data()] = ++(orderIndex);

    LoadingBar bar(40, "Computing Reeb space (SAT)...");
    int computedFaces = 1; // outside face has already been processed
    const int totalFaces = singularArrangement.arr.number_of_faces();

    while (false == traversalQueue.empty())
    {
        Halfedge_const_handle currentHalfEdge = traversalQueue.front();
        traversalQueue.pop();

        const int faceId = currentHalfEdge->face()->data();

        correspondenceGraph[currentHalfEdge->face()->data()] = this->fiberGraphs[currentHalfEdge->data().id].first.getUniqueComponents();

        // Iterate over neighbours of this face
        Halfedge_const_handle iteratorHalfEdge = currentHalfEdge;
        do
        {
            const int twinFaceId = iteratorHalfEdge->twin()->face()->data();

            // If this neighbour has not been visited
            if (-1 == order[twinFaceId])
            {
                seedFace(tetMesh, iteratorHalfEdge->twin());
                loopFace(tetMesh, iteratorHalfEdge->twin());

                traversalQueue.push(iteratorHalfEdge->twin());
                order[iteratorHalfEdge->twin()->face()->data()] = ++(orderIndex);
            }

            // @TODO should be else if, optimise
            // Compute correspondence with the neighbour
            if (order[faceId] < order[twinFaceId])
            {
                const FiberGraph &pgFace = this->fiberGraphs[iteratorHalfEdge->data().id].second;
                const FiberGraph &pgTwinFace = this->fiberGraphs[iteratorHalfEdge->twin()->data().id].first;

                const std::vector<std::pair<int, int>> componentPairs = pgFace.establishCorrespondence(
                        tetMesh, 
                        this->edgeCrossingSegments[iteratorHalfEdge->data().id],
                        pgTwinFace
                        );

                for (const auto &[componentIdA, componentIdB] : componentPairs)
                {
                    this->correspondenceGraphDS.unify(componentIdA, componentIdB);
                }
            }

            // Cache the first fiber graph for the first halfEdge of the current face
            if (false == currentHalfEdge->face()->is_unbounded() && iteratorHalfEdge->data().id == currentHalfEdge->face()->outer_ccb()->data().id)
            {
                //this->representativeFiberGraphs[iteratorHalfEdge->face()->data()] = this->fiberGraphs[iteratorHalfEdge->data().id].first;
                this->representativeFiberGraphSeeds[iteratorHalfEdge->face()->data()] = this->fiberGraphs[iteratorHalfEdge->data().id].first.getFiberGraphSeeds();

                //std::cout << "The fiber graph seeds are \n";

                //for (const auto [tId, cId] : this->representativeFiberGraphSeeds[iteratorHalfEdge->face()->data()])
                //{
                    //printf("Triangle is %d and component is %d.", tId, cId);
                //}

            }

            // We no longer need the preimage graphs, we can clear them
            if (false == cachePreimageGraphs)
            {
                this->fiberGraphs[iteratorHalfEdge->data().id].first = FiberGraph();
                this->fiberGraphs[iteratorHalfEdge->data().id].second = FiberGraph();
            }

            iteratorHalfEdge = iteratorHalfEdge->next();

        } while (iteratorHalfEdge != currentHalfEdge);


        computedFaces++;
        bar.update((100 * computedFaces) / totalFaces);
    }


    // Postprocessing
    this->correspondenceGraphDS.finalise();
    this->numberOfSheets = this->correspondenceGraphDS.countComponents();
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
            singularSegmentsHalfEdgeIds.emplace_back(he->data().id);
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

            this->edgeRegionSegments[he->data().id].push_back({regularSegmentHandle.originatingRegularEdgeId, true});
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

        auto &intersectingSegments = this->edgeRegionSegments[halfEdge->data().id];

        std::sort(intersectingSegments.begin(), intersectingSegments.end(), cmp);
    }
}



void ReebSpace2::computeEdgeRegionSegments3(const TetMesh &tetMesh, Arrangement &singularArrangement)
{
    std::vector<std::vector<int>> regularEdgeIntersections(tetMesh.edgeSingularTypes.size());

    // For each regular segment, compute its intersections with all singular segment
    //
    #pragma omp parallel for schedule(dynamic)
    for (const auto &edge : tetMesh.edges) 
    {
        if (tetMesh.edgeSingularTypes.at(edge) == 1)
        {
            const int edgeId = tetMesh.edgeIndices.at(edge);

            Arrangement_2::X_monotone_curve_2 segmentMonotoneCurve(
                    singularArrangement.arrangementPoints[edge[0]],
                    singularArrangement.arrangementPoints[edge[1]]);

            std::vector<Arrangement_2::Vertex_handle> vrts;
            std::vector<Arrangement_2::Halfedge_handle> hedges;

            CGAL::zone(singularArrangement.arr, segmentMonotoneCurve, CGAL::dispatch_or_drop_output<Arrangement_2::Vertex_handle, Arrangement_2::Halfedge_handle>(std::back_inserter(vrts), std::back_inserter(hedges)));


            for (auto h : hedges)
            {
                if (&(*h) > &(*h->twin()))
                {
                    h=h->twin();
                }

                regularEdgeIntersections[edgeId].push_back(h->data().id);
            }
        }
    }


    // Invert regularEdgeIntersections and save all regular segment intersections per half-edge in \bar{A}
    //
    this->edgeRegionSegments.resize(singularArrangement.arr.number_of_halfedges());

    for (const auto &[edge, type] : tetMesh.edgeSingularTypes) 
    {
        if (type == 1)
        {
            const int regularEdgeId = tetMesh.edgeIndices.at(edge);

            for (const auto &singularHalfEdgeId : regularEdgeIntersections[regularEdgeId])
            {
                this->edgeRegionSegments[singularHalfEdgeId].push_back({regularEdgeId, true});
            }
        }
    }


    // Temporary make a vector of half-edges so that we can use OpenMP next
    //
    std::vector<Halfedge_const_handle> halfEdgeVector;
    halfEdgeVector.reserve(singularArrangement.arr.number_of_halfedges());
    for (auto halfEdge = singularArrangement.arr.halfedges_begin(); halfEdge != singularArrangement.arr.halfedges_end(); ++halfEdge) 
    {
        halfEdgeVector.emplace_back(halfEdge);
    }

    // Sort the regular segment intersection points per half-edge
    //
    #pragma omp parallel for schedule(dynamic)
    for (const auto &halfEdge : halfEdgeVector) 
    {
        std::vector<std::pair<K::FT, int>> segmentRegionsOrdered;
        segmentRegionsOrdered.reserve(this->edgeRegionSegments[halfEdge->data().id].size());

        const Segment_2 singularSegment(halfEdge->source()->point(), halfEdge->target()->point());

        for (const auto &[regularEdgeId, boolDirection] : this->edgeRegionSegments[halfEdge->data().id])
        {
                const Segment_2 regularSegment(
                        singularArrangement.arrangementPoints[tetMesh.edges[regularEdgeId][0]],
                        singularArrangement.arrangementPoints[tetMesh.edges[regularEdgeId][1]]
                        );

                K::FT alpha = CGAL::Intersections::internal::s2s2_alpha(
                        halfEdge->target()->point().x(), halfEdge->target()->point().y(),
                        halfEdge->source()->point().x(), halfEdge->source()->point().y(),
                        regularSegment.source().x(), regularSegment.source().y(),
                        regularSegment.target().x(), regularSegment.target().y());

                segmentRegionsOrdered.emplace_back(alpha, regularEdgeId);
        }

        std::sort(segmentRegionsOrdered.begin(), segmentRegionsOrdered.end());

        for (int i = 0 ; i < segmentRegionsOrdered.size() ; i++)
        {
            edgeRegionSegments[halfEdge->data().id][i].first = segmentRegionsOrdered[i].second;
        }

    }


    // Timings on this : 2.5, 2.3 

    // Sort the regular segment intersection points per half-edge
    //
    //#pragma omp parallel for schedule(dynamic)
    //for (const auto &halfEdge : halfEdgeVector) 
    //{
        //std::map<K::FT, int> segmentRegionsOrdered;

        //const Segment_2 singularSegment(halfEdge->source()->point(), halfEdge->target()->point());

        //for (const auto &[regularEdgeId, boolDirection] : this->edgeRegionSegments3[halfEdge->data().id])
        //{
                //const Segment_2 regularSegment(
                        //singularArrangement.arrangementPoints[tetMesh.edges[regularEdgeId][0]],
                        //singularArrangement.arrangementPoints[tetMesh.edges[regularEdgeId][1]]
                        //);

                //K::FT alpha = CGAL::Intersections::internal::s2s2_alpha(
                        //halfEdge->target()->point().x(), halfEdge->target()->point().y(),
                        //halfEdge->source()->point().x(), halfEdge->source()->point().y(),
                        //regularSegment.source().x(), regularSegment.source().y(),
                        //regularSegment.target().x(), regularSegment.target().y());


                //segmentRegionsOrdered[alpha] = regularEdgeId;
        //}


        //edgeRegionSegments3[halfEdge->data().id].clear();

        //for (const auto& [alpha, regularEdgeId] : segmentRegionsOrdered)
        //{
            //edgeRegionSegments3[halfEdge->data().id].push_back({regularEdgeId, true});
        //}

    //}
}


void ReebSpace2::determineEdgeRegionSegmentsOrientation(const TetMesh &tetMesh, Arrangement &singularArrangement)
{

    for (auto halfEdge = singularArrangement.arr.halfedges_begin(); halfEdge != singularArrangement.arr.halfedges_end(); ++halfEdge) 
    {
        auto &intersectingSegments = this->edgeRegionSegments[halfEdge->data().id];

        const Point_2 &a = halfEdge->source()->point();
        const Point_2 &b = halfEdge->target()->point();

        for (std::pair<int, bool> &intersectingSegment : intersectingSegments)
        {
            const int segmentId = intersectingSegment.first;
            const std::array<int, 2> edge = tetMesh.edges.at(segmentId);
            const int segmentSourceId = edge[0];
            const int segmentTargetId = edge[1];
            const Point_2 &c = singularArrangement.arrangementPoints[segmentSourceId];
            const Point_2 &d = singularArrangement.arrangementPoints[segmentTargetId];

            // Sanity checks for degenerate cases.
            //
            if (
                    CGAL::orientation(a, b, c) == CGAL::COLLINEAR ||
                    CGAL::orientation(a, b, d) == CGAL::COLLINEAR
                    )
            {
                throw std::runtime_error("Input is degenerate, three points are collienar.");
            }

            if (
                    CGAL::orientation(a, b, c) == CGAL::orientation(a, b, d)
               )
            {
                //throw std::runtime_error("Input is degenerate, segments ab and cd do not intersect.");
            }

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
            if (CGAL::orientation(a, b, c) == CGAL::RIGHT_TURN)
            {
                intersectingSegment.second = false;
            }

            //
            // If the orientation this way around, okay, that's what we expect
            //   (regular segment)
            //          c 
            //           \
            //            \
            // a ----------\--------- b (singular segment)
            //              \
            //               \
            //                d
            //
        }

        auto &intersectingSegmentsTwin = this->edgeRegionSegments[halfEdge->twin()->data().id];

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
    }
}





























bool ReebSpace2::ifSegmentInHalfEdgeRegion(const Halfedge_const_handle &halfEdge, const Segment_2 &segment)
{
    // Previous version
    //
    //Arrangement_2::Halfedge_around_vertex_const_circulator &halfEdgeCirculator
    //const auto halfEdge = halfEdgeCirculator;
    //auto halfEdgeCirculatorPrevious = halfEdgeCirculator;
    //const auto halfEdgeNext = ++halfEdgeCirculatorPrevious;
    //assert(halfEdgeCirculator != nullptr); 

    const auto halfEdgeNext = halfEdge->next()->twin();


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
            vertexRegionSegments[regionHalfEdgeHandle->data().id].push_back({mySegment.originatingRegularEdgeId, true});

            //std::cout << "Added segment [" << tetMesh.edges[mySegment.originatingEdge][0] << ", " << tetMesh.edges[mySegment.originatingEdge][1] << "]\n";
            //std::cout << "Between " << singularArrangement.arrangementPoints[tetMesh.edges[mySegment.originatingEdge][0]] << " and " << singularArrangement.arrangementPoints[tetMesh.edges[mySegment.originatingEdge][0]] << std::endl;
        }


        const auto itTarget = singularArrangement.arrangementPointHandles.find(segment.target());
        if (itTarget != singularArrangement.arrangementPointHandles.end())
        {
            //std::cout << "The target point is : " << itTarget->second->point() << " | " <<  singularArrangement.arrangementPointIndices.at(itTarget->second->point()) << std::endl;
            Halfedge_const_handle regionHalfEdgeHandle = getSegmentRegion(itTarget->second, segment);
            vertexRegionSegments[regionHalfEdgeHandle->data().id].push_back({mySegment.originatingRegularEdgeId, true});

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

        auto &segmentIds = this->vertexRegionSegments[halfEdge->data().id];
        std::sort(segmentIds.begin(), segmentIds.end(), cmp);
    }
}


void ReebSpace2::determineVertexRegionSegmentsOrientation(const TetMesh &tetMesh, Arrangement &singularArrangement)
{
    for (auto halfEdge = singularArrangement.arr.halfedges_begin(); halfEdge != singularArrangement.arr.halfedges_end(); ++halfEdge) 
    {
        auto &intersectingSegments = this->vertexRegionSegments[halfEdge->data().id];

        // Each vertex is guaratneed to be in the singular arrangement
        const Point_2 &vertex = halfEdge->target()->point();
        const int vertexMeshId = singularArrangement.arrangementPointIndices[vertex];

        for (std::pair<int, bool> &intersectingSegment : intersectingSegments)
        {
            const int segmentId = intersectingSegment.first;
            const std::array<int, 2> edge = tetMesh.edges[segmentId];

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
            //   This is the default case, so no need to do anything
            //      
        }

    }
}








void ReebSpace2::determineEdgeCrossingSegmentsOriantation(const TetMesh &tetMesh, Arrangement &singularArrangement)
{
    this->edgeCrossingSegments.resize(singularArrangement.arr.number_of_halfedges());

    for (auto he = singularArrangement.arr.halfedges_begin(); he != singularArrangement.arr.halfedges_end(); ++he)
    {
        const Segment_2 &segment = *singularArrangement.arr.originating_curves_begin(he);

        const int aIndex = singularArrangement.arrangementPointIndices.at(segment.source());
        const int bIndex = singularArrangement.arrangementPointIndices.at(segment.target());

        // Sanity check
        assert(aIndex < bIndex);

        const std::array<int, 2> edge = {aIndex, bIndex};

        // Check to see if the segment and half edge have the same orientation
        const bool isSegmentLeftToRight = segment.source() < segment.target(); 
        const bool isCurrentHalfEdgeLeftToRight = (he->direction() == CGAL::ARR_LEFT_TO_RIGHT);

        const int edgeId = tetMesh.edgeIndices.at(edge);

        edgeCrossingSegments[he->data().id].first = edgeId;

        // The half edge has the same direction as the original edge
        if (isSegmentLeftToRight == isCurrentHalfEdgeLeftToRight)
        {
            edgeCrossingSegments[he->data().id].second = false;
        }
        else
        {
            edgeCrossingSegments[he->data().id].second = true;
        }
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


bool ReebSpace2::unitTestCompareFiberGraphs(const TetMesh &tetMesh, Arrangement &singularArrangement, Arrangement &regularArrangement, ReebSpace &rs)
{
    //for (const auto& [halfEdge, preimageGraphs] : this->preimageGraphsCached)
    for (auto halfEdge = singularArrangement.arr.halfedges_begin(); halfEdge != singularArrangement.arr.halfedges_end(); ++halfEdge)
    {
        //std::cout << "Checking graph equality..." << std::endl;
        std::pair<Halfedge_const_handle, Halfedge_const_handle> foundHalfEdges =  findHalfEdges(halfEdge, singularArrangement, regularArrangement);
        std::pair<Face_const_handle, Face_const_handle> foundFaces =  {foundHalfEdges.first->face(), foundHalfEdges.second->face()};
        std::pair<int, int> foundFaceIds =  {regularArrangement.arrangementFacesIdices.at(foundFaces.first), regularArrangement.arrangementFacesIdices.at(foundFaces.second)};

        std::pair<FiberGraph, FiberGraph> &singularFiberGraphs = this->fiberGraphs[halfEdge->data().id];
        std::pair<DisjointSet<int>, DisjointSet<int>> regularFiberGraphs = {rs.preimageGraphs[foundFaceIds.first], rs.preimageGraphs[foundFaceIds.second]};

        //assert(singularPreimageGraphs.first == regularPreimageGraphs.first);
        //assert(singularPreimageGraphs.second == regularPreimageGraphs.second);


        if (false == singularFiberGraphs.first.areEqual(regularFiberGraphs.first) || false == singularFiberGraphs.second.areEqual(regularFiberGraphs.second))
        {
            printf("\n------------------------------------------------------------------------------------\n");
            std::cout << "SingularHalf-edge is [" << halfEdge->source()->point() << "] -> [" << halfEdge->target()->point() << "]";
            printf("The IDs are %d -> %d\n", singularArrangement.arrangementPointIndices[halfEdge->source()->point()], singularArrangement.arrangementPointIndices[halfEdge->target()->point()]);
            printf("\n------------------------------------------------------------------------------------\n");


            std::cout << "First singular preimage graph is: \n";
            singularFiberGraphs.first.printByRoot();

            std::cout << "Second singular preimage graph is: \n";
            singularFiberGraphs.second.printByRoot();


            printf("\n------------------------------------------------------------------------------------\n");
            std::cout << "Regular-edge 1 is [" << foundHalfEdges.first->source()->point() << "] -> [" << foundHalfEdges.first->target()->point() << "]\n";
            std::cout << "Regular-edge 2 is [" << foundHalfEdges.second->source()->point() << "] -> [" << foundHalfEdges.second->target()->point() << "]";
            printf("\n------------------------------------------------------------------------------------\n");

            std::cout << "First regular preimage graph is: \n";
            regularFiberGraphs.first.printByRoot();

            std::cout << "Second regular preimage graph is: \n";
            regularFiberGraphs.second.printByRoot();

            printf("\n\n\n");

            return false;
            //throw std::runtime_error("Preimage graphs do not match!");
        }
        else
        {
            //printf("\n------------------------------------------------------------------------------------\n");
            //std::cout << "SingularHalf-edge is [" << halfEdge->source()->point() << "] -> [" << halfEdge->target()->point() << "]";
            //printf("The IDs are %d -> %d\n", singularArrangement.arrangementPointIndices[halfEdge->source()->point()], singularArrangement.arrangementPointIndices[halfEdge->target()->point()]);
            //printf("\n------------------------------------------------------------------------------------\n");


            //std::cout << "First singular preimage graph is: \n";
            //singularFiberGraphs.first.printByRoot();

            //std::cout << "Second singular preimage graph is: \n";
            //singularFiberGraphs.second.printByRoot();


            //printf("\n------------------------------------------------------------------------------------\n");
            //std::cout << "Regular-edge 1 is [" << foundHalfEdges.first->source()->point() << "] -> [" << foundHalfEdges.first->target()->point() << "]\n";
            //std::cout << "Regular-edge 2 is [" << foundHalfEdges.second->source()->point() << "] -> [" << foundHalfEdges.second->target()->point() << "]";
            //printf("\n------------------------------------------------------------------------------------\n");

            //std::cout << "First regular preimage graph is: \n";
            //regularFiberGraphs.first.printByRoot();

            //std::cout << "Second regular preimage graph is: \n";
            //regularFiberGraphs.second.printByRoot();

            //printf("\n\n\n");
            //std::cout << "Fiber graphs are equal!\n" << std::endl;
        }
    }
    return true;
}





std::vector<std::pair<K::FT, int>> getIntersectedSegments(TetMesh &tetMesh, Arrangement &singularArrangement, const Segment_2 &controlSegment, const bool shouldSort)
{
    //Timer::start();
    std::vector<TreeAABB::Primitive_id> intersectedSegmentsAABB;
    singularArrangement.tree.all_intersected_primitives(controlSegment, std::back_inserter(intersectedSegmentsAABB));
    //Timer::stop("Computing AABB                         :");

    std::vector<std::pair<K::FT, int>> intersectedSegments;
    intersectedSegments.reserve(intersectedSegmentsAABB.size());

    const double cx1 = CGAL::to_double(controlSegment.target().x());
    const double cy1 = CGAL::to_double(controlSegment.target().y());
    const double cx2 = CGAL::to_double(controlSegment.source().x());
    const double cy2 = CGAL::to_double(controlSegment.source().y());

    //Timer::start();
    for (auto id : intersectedSegmentsAABB)
    {
        const Segment_2& s = *id;   // dereference iterator to get the original segment

        // Get the ID of the original segment.
        const int segmentIndex = id - singularArrangement.allSegments.begin();

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
        K::FT alpha = CGAL::Intersections::internal::s2s2_alpha(
                controlSegment.target().x(), controlSegment.target().y(),
                controlSegment.source().x(), controlSegment.source().y(),
                s.source().x(), s.source().y(),
                s.target().x(), s.target().y());

        //double alpha = CGAL::Intersections::internal::s2s2_alpha(
                //cx1, cy1, cx2, cy2,
                //CGAL::to_double(s.source().x()), 
                //CGAL::to_double(s.source().y()),
                //CGAL::to_double(s.target().x()), 
                //CGAL::to_double(s.target().y())
                //);

        intersectedSegments.emplace_back(alpha, segmentIndex);
    }
    //Timer::stop("Intersections Alphas                   :");

    //Timer::stop("Computed Alpha intersections           :");
    //Timer::stop("Computed AABB intersections in         :");

    if (shouldSort)
    {
        //Timer::start();
        std::sort(intersectedSegments.begin(), intersectedSegments.end());
        //Timer::stop("Sorting alpha intersections            :");
    }


    //for (const auto &[alpha, edgeId] : intersectedSegments)
    //{
        //const int singularType = tetMesh.edgeSingularTypes.at(tetMesh.edges.at(edgeId));

        //if (alpha < 1 && singularType != 1)
        //{
            //throw std::runtime_error("Control point segment intersects other singular segments.");
        //}

        ////std::cout << "Intersected segment with ID " << edgeId << " and type " << tetMesh.edgeSingularTypes.at(tetMesh.edges.at(edgeId)) << " and alpha " << alpha << std::endl;
    //}


    return intersectedSegments;
}





















































int ReebSpace2::computeFiberGraphReverse(TetMesh &tetMesh, Arrangement &singularArrangement, std::array<double, 2> controlPoint, std::set<int> initialTriangles)
{
    //std::cout << "\n\n----------------------------------------------------------------------\n\n";
    const Point_2 controlPointEPEC(controlPoint[0], controlPoint[1]);

    // 1. Compute the active face
    //
    Face_const_handle activeFace = singularArrangement.getActiveFace(controlPoint);

    if (activeFace->is_unbounded())
    {
        return {};
    }

    const int activeFaceId = singularArrangement.arrangementFacesIdices[activeFace];

    // 2. Compute the visibel vertex
    //
    //Timer::start();
    Point_2 closestHalfEdgeVertexPoint = singularArrangement.findVisibleVertex(activeFace, controlPointEPEC);
    //Timer::stop("Found visible point  in                :");


    // Set up the control Segment
    //
    const Segment_2 controlSegment(controlPointEPEC, closestHalfEdgeVertexPoint);

    if (controlSegment.squared_length() == 0.0)
    {
        std::cerr << "The segment has zero lenght!" << std::endl;

    }



    std::set<int> currentTriangles = initialTriangles;


    //
    // Go through the NEW segment
    //
    std::vector<std::pair<K::FT, int>> intersectedSegments = getIntersectedSegments(tetMesh, singularArrangement, controlSegment, true);









    //Timer::start();

    for (int i = 0 ; i < intersectedSegments.size() ; i++)
    {
        if (intersectedSegments[i].first == 1.0)
        {
            continue;
        }


        const int segmentId = intersectedSegments[i].second;
        bool typicalOrientation = false;

        //std::cout << "Intersected segment with ID " << segmentId << " and type " << tetMesh.edgeSingularTypes.at(tetMesh.edges.at(segmentId)) << " and alpha " << intersectedSegments[i].first << std::endl;

        // Change orientation in case we need to
        const std::array<int, 2> edge = tetMesh.edges.at(segmentId);
        const int segmentSourceId = edge[0];
        const int segmentTargetId = edge[1];
        const Point_2 &c = singularArrangement.arrangementPoints[segmentSourceId];
        const Point_2 &d = singularArrangement.arrangementPoints[segmentTargetId];

        //std::cout << c << " - " << d << std::endl;

        const Point_2 &a = closestHalfEdgeVertexPoint;
        const Point_2 &b = controlPointEPEC;
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
        if (CGAL::orientation(a, b, c) == CGAL::RIGHT_TURN)
        {
            typicalOrientation = true;
        }

        const std::vector<int> &minusTriangles = tetMesh.getMinusTriangles(segmentId, typicalOrientation);
        const std::vector<int> &plusTriangles = tetMesh.getPlusTriangles(segmentId, typicalOrientation);

        for (const auto &t : minusTriangles)
        {
            if (currentTriangles.contains(t))
            {
                currentTriangles.erase(t);
                currentTriangles.insert(*plusTriangles.begin());

                break;
            }
        }
    }

    //Timer::stop("Walking along new segment              :");




    //Timer::start();

    // Go up to closestHalfEdgeVertexPoint
    auto currentHalfEdge = activeFace->outer_ccb();


    std::vector<std::vector<std::pair<int, bool>>> intersectingEdges;
    intersectingEdges.push_back(this->edgeRegionSegments[currentHalfEdge->data().id]);

    ++currentHalfEdge;
    do
    {
        intersectingEdges.push_back(this->vertexRegionSegments[currentHalfEdge->prev()->data().id]);
        intersectingEdges.push_back(this->edgeRegionSegments[currentHalfEdge->data().id]);

        // If we have reached the visible vertex AND the control segment is in its region (not always the case with poinch points)
        if (currentHalfEdge->target()->point() == closestHalfEdgeVertexPoint && this->ifSegmentInHalfEdgeRegion(currentHalfEdge, controlSegment))
        {
            break;
        }
        else
        {
            ++currentHalfEdge;
        }


    } while (true);



    // 
    // Go through the vertex region
    //
    const auto &vertexRegion = this->vertexRegionSegments[currentHalfEdge->data().id];

    const int vertexMeshId = singularArrangement.arrangementPointIndices[closestHalfEdgeVertexPoint];

    for (int i = 0 ; i < vertexRegion.size() ; i++)
    {
        const std::pair<int, int> &intersectingSegment = vertexRegion[i];

        const int segmentId = intersectingSegment.first;

        const std::array<int, 2> edge = tetMesh.edges[segmentId];

        Point_2 b;
        if (edge[0] == vertexMeshId)
        {
            b = singularArrangement.arrangementPoints[edge[1]];
        }
        else
        {
            b = singularArrangement.arrangementPoints[edge[0]];
        }

        if (this->compareRegularSegments(currentHalfEdge, b, controlPointEPEC))
        {
            intersectingEdges.push_back({intersectingSegment});
        }
        else
        {
            break;
        }
    }

    //Timer::stop("Collecting loop intersections          :");


    //Timer::start();
    for (int i = intersectingEdges.size() - 1 ; i >= 0 ; i--)
    {
        auto &intersectingEdgesInner = intersectingEdges[i];
        for (int j = intersectingEdgesInner.size() - 1 ; j >= 0 ; j--)
        {
            const std::vector<int> &minusTriangles = tetMesh.getMinusTriangles(intersectingEdgesInner[j].first, !intersectingEdgesInner[j].second);
            const std::vector<int> &plusTriangles = tetMesh.getPlusTriangles(intersectingEdgesInner[j].first, !intersectingEdgesInner[j].second);

            for (const auto &t : minusTriangles)
            {
                if (currentTriangles.contains(t))
                {
                    currentTriangles.erase(t);
                    currentTriangles.insert(*plusTriangles.begin());

                    break;
                }
            }
        }
    }
    //Timer::stop("Collecting loop intersections          :");



    FiberGraph pg = this->representativeFiberGraphs[activeFace->data()];

    for (const auto &t : currentTriangles)
    {
        if (pg.componentRoot.contains(t))
        {
            return this->correspondenceGraphDS.find(pg.componentRoot.at(t));
        }
    }

    return -1;
}



FiberGraph ReebSpace2::computeFiberGraph(TetMesh &tetMesh, Arrangement &singularArrangement, std::array<double, 2> controlPoint)
{
    const Point_2 controlPointEPEC(controlPoint[0], controlPoint[1]);

    // 1. Compute the active face
    //
    Face_const_handle activeFace = singularArrangement.getActiveFace(controlPoint);

    if (activeFace->is_unbounded())
    {
        return {};
    }

    const int activeFaceId = singularArrangement.arrangementFacesIdices[activeFace];

    // 2. Compute the visibel vertex
    //
    //Timer::start();
    Point_2 closestHalfEdgeVertexPoint = singularArrangement.findVisibleVertex(activeFace, controlPointEPEC);
    //Timer::stop("Found visible point  in                :");


    // Set up the control Segment
    //
    const Segment_2 controlSegment(controlPointEPEC, closestHalfEdgeVertexPoint);

    if (controlSegment.squared_length() == 0.0)
    {
        std::cerr << "The segment has zero lenght!" << std::endl;

    }

    int graphUpdates = 0;

    //Timer::start();

    // Compute the fiber graph
    //

    // Go up to closestHalfEdgeVertexPoint
    auto currentHalfEdge = activeFace->outer_ccb();

    //std::cout << "Starting half-edge is " << currentHalfEdge->source()->point() << " -> " << currentHalfEdge->target()->point() << std::endl;

    FiberGraph pg = this->representativeFiberGraphs[activeFace->data()];
    //pg.printByRoot();
    pg.updateComponentsRegular(tetMesh, this->edgeRegionSegments[currentHalfEdge->data().id]);
    //pg.printByRoot();

    ++currentHalfEdge;
    do
    {
        //std::cout << "Next half-edge is " << currentHalfEdge->source()->point() << " -> " << currentHalfEdge->target()->point() << std::endl;
        pg.updateComponentsRegular(tetMesh, this->vertexRegionSegments[currentHalfEdge->prev()->data().id]);
        //pg.printByRoot();
        pg.updateComponentsRegular(tetMesh, this->edgeRegionSegments[currentHalfEdge->data().id]);
        //pg.printByRoot();

        graphUpdates++;
        graphUpdates++;

        // If we have reached the visible vertex AND the control segment is in its region (not always the case with poinch points)
        if (currentHalfEdge->target()->point() == closestHalfEdgeVertexPoint && this->ifSegmentInHalfEdgeRegion(currentHalfEdge, controlSegment))
        {
            break;
        }
        else
        {
            ++currentHalfEdge;
        }


    } while (true);


    //std::cout << "Final half-edge is " << currentHalfEdge->source()->point() << " -> " << currentHalfEdge->target()->point() << std::endl;


    // Go to contorlSegment in the region of the vertex
    //
    const auto &vertexRegion = this->vertexRegionSegments[currentHalfEdge->data().id];

    const int vertexMeshId = singularArrangement.arrangementPointIndices[closestHalfEdgeVertexPoint];

    for (int i = 0 ; i < vertexRegion.size() ; i++)
    {
        const std::pair<int, int> &intersectingSegment = vertexRegion[i];

        const int segmentId = intersectingSegment.first;

        const std::array<int, 2> edge = tetMesh.edges[segmentId];

        Point_2 b;
        if (edge[0] == vertexMeshId)
        {
            b = singularArrangement.arrangementPoints[edge[1]];
        }
        else
        {
            b = singularArrangement.arrangementPoints[edge[0]];
        }

        if (this->compareRegularSegments(currentHalfEdge, b, controlPointEPEC))
        {
            pg.updateComponentsRegular(tetMesh, {intersectingSegment});
            graphUpdates++;
            //std::cout << "It's happening!" << std::endl;
        }
        else
        {
            break;
        }
    }


    //
    // Got along the control segment 


    //std::cout << std::endl;
    //std::cout << std::endl;



    //Timer::start();
    std::vector<std::pair<K::FT, int>> intersectedSegments = getIntersectedSegments(tetMesh, singularArrangement, controlSegment, true);
    //Timer::stop("AABB search 1                          :");

    for (int i = intersectedSegments.size() - 1 ; i >= 0 ; i--)
    {
        if (intersectedSegments[i].first == 1.0)
        {
            continue;
        }


        const int segmentId = intersectedSegments[i].second;
        bool typicalOrientation = true;

        //std::cout << "Intersected segment with ID " << segmentId << " and type " << tetMesh.edgeSingularTypes.at(tetMesh.edges.at(segmentId)) << " and alpha " << intersectedSegments[i].first << std::endl;

        // Change orientation in case we need to
        const std::array<int, 2> edge = tetMesh.edges.at(segmentId);
        const int segmentSourceId = edge[0];
        const int segmentTargetId = edge[1];
        const Point_2 &c = singularArrangement.arrangementPoints[segmentSourceId];
        const Point_2 &d = singularArrangement.arrangementPoints[segmentTargetId];

        //std::cout << c << " - " << d << std::endl;

        const Point_2 &a = closestHalfEdgeVertexPoint;
        const Point_2 &b = controlPointEPEC;
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
        if (CGAL::orientation(a, b, c) == CGAL::RIGHT_TURN)
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






FiberGraph ReebSpace2::computeFiberGraph2(TetMesh &tetMesh, Arrangement &singularArrangement, std::array<double, 2> controlPoint)
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



    FiberGraph pg = this->representativeFiberGraphs[activeFace->data()];

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

        pg.updateComponentsRegular(tetMesh, this->edgeRegionSegments[currentHalfEdge->data().id]);
        graphUpdates++;
        pg.updateComponentsRegular(tetMesh, this->vertexRegionSegments[currentHalfEdge->data().id]);
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


    for (const std::pair<int, bool> &regionSegment :  this->edgeRegionSegments[currentHalfEdge->data().id])
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




FiberGraph ReebSpace2::computeFiberGraph3(TetMesh &tetMesh, Arrangement &singularArrangement, const Segment_2 &controlSegment, const K::FT &pointAlpha, const std::vector<std::tuple<K::FT, int, int>> &intersectedSegments)
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
        return this->computeFiberGraph2(tetMesh, singularArrangement, controlPointDouble);
    }


    //const int destinationSegmentId = std::get<1>(intersectedSegments.back());



    //std::cout << std::endl << "Destination segment is " << singularArrangement.arrangementPoints.at(tetMesh.edges[destinationSegmentId][0]) << " -> " << singularArrangement.arrangementPoints.at(tetMesh.edges[destinationSegmentId][1]) << " and ID = " << destinationSegmentId << std::endl << std::endl;


    FiberGraph pg = this->representativeFiberGraphs[activeFace->data()];

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

        pg.updateComponentsRegular(tetMesh, this->edgeRegionSegments[currentHalfEdge->data().id]);
        graphUpdates++;
        pg.updateComponentsRegular(tetMesh, this->vertexRegionSegments[currentHalfEdge->data().id]);
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


    for (const std::pair<int, bool> &regionSegment :  this->edgeRegionSegments[currentHalfEdge->data().id])
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




        const std::vector<int> &minusTriangles = tetMesh.getMinusTriangles(segmentId, typicalOrientation);
        const std::vector<int> &plusTriangles = tetMesh.getPlusTriangles(segmentId, typicalOrientation);



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


        pg.updateComponentsRegular(tetMesh, {{segmentId, typicalOrientation}});
        graphUpdates++;

    }


    //std::cout << "We have performed " << graphUpdates << " graph updates.\n";



    return pg;
}




std::vector<std::pair<int, int>> ReebSpace2::computeSeedFibers(TetMesh &tetMesh, Arrangement &singularArrangement, std::array<double, 2> controlPoint, const std::set<int> &selectedSheetIds)
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


    //FiberGraph pg = this->representativeFiberGraphs[activeFace->data()];

    std::vector<std::pair<int, int>> fiberSeeds;


    // If there are selected sheets, only use their fiber components
    //
    for (const auto &[triangleId, componentId] : this->representativeFiberGraphSeeds[activeFace->data()])
    {
        if (!selectedSheetIds.empty())
        {
            const int sheetId = this->correspondenceGraphDS.find(componentId);
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

        //pg.updateComponentsRegular(tetMesh, this->edgeRegionSegments[currentHalfEdge->data().id]);
        FiberGraph::updateSeedsRegular(tetMesh, this->edgeRegionSegments[currentHalfEdge->data().id], fiberSeeds);
        graphUpdates++;

        //pg.updateComponentsRegular(tetMesh, this->vertexRegionSegments[currentHalfEdge->data().id]);
        FiberGraph::updateSeedsRegular(tetMesh, this->vertexRegionSegments[currentHalfEdge->data().id], fiberSeeds);
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


    for (const std::pair<int, bool> &regionSegment :  this->edgeRegionSegments[currentHalfEdge->data().id])
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






std::vector<std::pair<int, int>> ReebSpace2::computeSeedFibersGivenLine(TetMesh &tetMesh, Arrangement &singularArrangement, const Segment_2 &controlSegment, const K::FT &pointAlpha, const std::vector<std::tuple<K::FT, int, int>> &intersectedSegments)
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
        return this->computeSeedFibers(tetMesh, singularArrangement, controlPointDouble, {});
    }


    //const int destinationSegmentId = std::get<1>(intersectedSegments.back());



    //std::cout << std::endl << "Destination segment is " << singularArrangement.arrangementPoints.at(tetMesh.edges[destinationSegmentId][0]) << " -> " << singularArrangement.arrangementPoints.at(tetMesh.edges[destinationSegmentId][1]) << " and ID = " << destinationSegmentId << std::endl << std::endl;


    //FiberGraph pg = this->representativeFiberGraphs[activeFace->data()];
    std::vector<std::pair<int, int>> fiberSeeds = this->representativeFiberGraphSeeds[activeFace->data()];

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

        //pg.updateComponentsRegular(tetMesh, this->edgeRegionSegments[currentHalfEdge->data().id]);
        FiberGraph::updateSeedsRegular(tetMesh, this->edgeRegionSegments[currentHalfEdge->data().id], fiberSeeds);
        graphUpdates++;

        //pg.updateComponentsRegular(tetMesh, this->vertexRegionSegments[currentHalfEdge->data().id]);
        FiberGraph::updateSeedsRegular(tetMesh, this->vertexRegionSegments[currentHalfEdge->data().id], fiberSeeds);
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


    for (const std::pair<int, bool> &regionSegment :  this->edgeRegionSegments[currentHalfEdge->data().id])
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


    //std::cout << "We have performed " << graphUpdates << " graph updates.\n";



    return fiberSeeds;
}
