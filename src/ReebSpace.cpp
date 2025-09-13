#include "./CGALTypedefs.h"

#include <cassert>
#include <cstddef>
#include <map>
#include <queue>
#include <set>
#include <iterator>
#include <unordered_map>

#include "./Timer.h"
#include "./DisjointSet.h"
#include "./ReebSpace.h"

std::pair<std::vector<int>, std::vector<int>> ReebSpace::getMinusPlusTrianglesIndex(const TetMesh &tetMesh, const Arrangement &arrangement, const Arrangement_2::Halfedge_const_handle currentHalfEdge)
{
    const Segment_2 &segment = *arrangement.arr.originating_curves_begin(currentHalfEdge);
    //std::cout << "Half-edge   from: " << currentHalfEdge->source()->point() << " to " << currentHalfEdge->target()->point() << std::endl;
    //std::cout << "Source-edge from: " << segment.source() << " to " << segment.target() << std::endl;

    const int aIndex = arrangement.arrangementPointIndices.at(segment.source());
    const int bIndex = arrangement.arrangementPointIndices.at(segment.target());

    // Sanity check
    assert(aIndex < bIndex);

    const std::array<int, 2> edge = {aIndex, bIndex};

    // Check to see if the segment and half edge have the same orientation
    const bool isSegmentLeftToRight = segment.source() < segment.target(); 
    const bool isCurrentHalfEdgeLeftToRight = (currentHalfEdge->direction() == CGAL::ARR_LEFT_TO_RIGHT);

    // The half edge has the same direction as the original edge
    if (isSegmentLeftToRight == isCurrentHalfEdgeLeftToRight)
    {
        return {tetMesh.upperStarTriangles.at(edge), tetMesh.lowerStarTriangles.at(edge)};
    }
    else
    {
        return {tetMesh.lowerStarTriangles.at(edge), tetMesh.upperStarTriangles.at(edge)};
    }
}



void ReebSpace::testTraverseArrangement(const Arrangement &arrangement)
{
    // Find the unbounded face (hold the boundary of the arrangement)
    Face_const_handle outerFace;

    // Iterate over all faces and find the unbounded one
    for (Face_const_iterator fit = arrangement.arr.faces_begin(); fit != arrangement.arr.faces_end(); ++fit) 
    {
        if (fit->is_unbounded()) 
        {
            // If this has already been set, we have two outerFace, this should never happen.
            assert(outerFace == Face_const_handle());
            outerFace = fit;
            break;
        }
    }

    // Sanity check, make sure the outer face is simple
    assert(outerFace->number_of_inner_ccbs() == 1);
    assert(outerFace->number_of_outer_ccbs() == 0);
    assert(outerFace->number_of_holes() == 1);
    assert(outerFace->number_of_isolated_vertices() == 0);

    // Get the the first half edge of the outerFace (could be any edge, this is a matter of convention)
    Halfedge_const_handle outerHalfEdge = *outerFace->holes_begin();

    // Make sure there is only one originating curve, something has gone wrong otherwise (edge overlap)
    assert(std::distance(arrangement.arr.originating_curves_begin(outerHalfEdge), arrangement.arr.originating_curves_end(outerHalfEdge)) == 1);

    // Starting from an outer half edge
    std::queue<Arrangement_2::Halfedge_const_handle> traversalQueue;
    traversalQueue.push(outerHalfEdge);

    // Which faces are visited
    std::set<Arrangement_2::Face_const_handle> visited;
    visited.insert(outerHalfEdge->face());

    while (false == traversalQueue.empty())
    {
        // Pop an half edge out
        Arrangement_2::Halfedge_const_handle currentHalfEdge = traversalQueue.front();
        traversalQueue.pop();
        Arrangement_2::Face_const_handle currentFace = currentHalfEdge->face();

        // Get the twin
        Arrangement_2::Halfedge_const_handle twin = currentHalfEdge->twin();
        Arrangement_2::Face_const_handle twinFace = twin->face();

        // Get ids of the current face and the twin face
        int currentFaceID = arrangement.arrangementFacesIdices.at(currentFace);
        int twinFaceID = arrangement.arrangementFacesIdices.at(twinFace);

        // If we have never visited this face, then we have never visited any of the half edges.
        if (visited.find(twinFace) == visited.end())
        {
            //printf("NEW FACE ------------------------------------------ %d -> %d \n", data->arrangementFacesIdices[currentFace], data->arrangementFacesIdices[twinFace]);
            visited.insert(twinFace);

            // Sanity check we should only have onbounded face, the outside face.
            assert(false == twinFace->is_unbounded());


            Arrangement_2::Ccb_halfedge_const_circulator start = twinFace->outer_ccb();
            Arrangement_2::Ccb_halfedge_const_circulator curr = start;

            do {
                traversalQueue.push(curr);

                // Make sure there is only one originating curve (sanity check)
                //const Segment_2 &segment = *data->arr.originating_curves_begin(curr);
                //std::cout << "Half-edge   from: " << curr->source()->point() << " to " << curr->target()->point() << std::endl;
                //std::cout << "Source-edge from: " << segment.source() << " to " << segment.target() << std::endl;
                //printf("The original indices are %d and %d", data->arrangementPointsIdices[segment.source()], data->arrangementPointsIdices[segment.target()]);
                //printf("\n\n");

                ++curr;
            } while (curr != start);
        }
    }
}


void ReebSpace::computeTwinFacePreimageGraph(const TetMesh &tetMesh, const Arrangement &arrangement, const Arrangement_2::Halfedge_const_handle &currentHalfEdge)
{
    Face_const_handle currentFace = currentHalfEdge->face();
    const int currentFaceID = arrangement.arrangementFacesIdices.at(currentFace);

    Face_const_handle twinFace = currentHalfEdge->twin()->face();
    const int twinFaceID = arrangement.arrangementFacesIdices.at(twinFace);

    const auto [minusTriangles, plusTriangles] = ReebSpace::getMinusPlusTrianglesIndex(tetMesh, arrangement, currentHalfEdge);

    std::set<int> preimageGraph;
    for (const auto &[triangleId, internalIndex] : preimageGraphs[currentFaceID].data)
    {
        preimageGraph.insert(triangleId);
    }

    for (const auto &triangle: minusTriangles)
    {
        preimageGraph.erase(triangle);
    }

    for (const auto &triangle: plusTriangles)
    {
        preimageGraph.insert(triangle);
    }

    this->preimageGraphs[twinFaceID].initialize(preimageGraph);
    for (const auto &[t1, id1] : this->preimageGraphs[twinFaceID].data)
    {
        for (const auto &t2 : tetMesh.tetIncidentTriangles[t1])
        {
            if (this->preimageGraphs[twinFaceID].data.contains(t2))
            {
                this->preimageGraphs[twinFaceID].unionElements(t1, t2);
            }
        }
    }

    // Finaly make sure every element in the disjoint set points to their root
    this->preimageGraphs[twinFaceID].finalise();
}



void ReebSpace::determineCorrespondence(const TetMesh &tetMesh, const Arrangement &arrangement, const Arrangement_2::Halfedge_const_handle &halfEdge)
{
    Face_const_handle face = halfEdge->face();
    const int faceId = arrangement.arrangementFacesIdices.at(face);

    Face_const_handle twinFace = halfEdge->twin()->face();
    const int twinFaceId = arrangement.arrangementFacesIdices.at(twinFace);

    const auto [minusTriangles, plusTriangles] = ReebSpace::getMinusPlusTrianglesIndex(tetMesh, arrangement, halfEdge);

    std::set<int> affectedComponents;
    for (const int &triangle : minusTriangles)
    {
        affectedComponents.insert(this->preimageGraphs[faceId].findElement(triangle));
    }

    std::set<int> twinAffectedComponents;
    for (const int &triangle : plusTriangles)
    {
        twinAffectedComponents.insert(preimageGraphs[twinFaceId].findElement(triangle));

    }

    // Definite edge
    //if (activeRootsFace.size() == 0 || activeRootsTwinFace.size() == 0)
    //{

        //data->jacobiType[originatingEdge] = 0;
    //}
    //// Reeb-regular (regualr of indefinite of type 1-1)
    //else if (activeRootsFace.size() == 1 && activeRootsTwinFace.size() == 1)
    //{

        //data->jacobiType[originatingEdge] = 1;
    //}
    //// Indefinite edge
    //else
    //{
        //data->jacobiType[originatingEdge] = 2;

    //}

    // This is a regular edge then connect the two active components, for singular edge we skip this, they are considered new
    if (affectedComponents.size() == 1 && twinAffectedComponents.size() == 1)
    {
        const int affectedComponentId = *affectedComponents.begin();
        const int twinAffectedComponentId = *twinAffectedComponents.begin();
        this->correspondenceGraph.unionElements({faceId, affectedComponentId}, {twinFaceId, twinAffectedComponentId});
    }

    // Link together all other connected components, there is a one-to-one correspondence between them
    for (const auto &[triangleId, internalId] : this->preimageGraphs[faceId].data)
    {
        // The root of the triangle in the face
        const int componentId = this->preimageGraphs[faceId].findElement(triangleId);

        // We have already deal with the active fiber components
        if (affectedComponents.contains(componentId)) { continue; }

        // The root of the triangle in the twin face
        const int twinComponentId = this->preimageGraphs[twinFaceId].findElement(triangleId);

        this->correspondenceGraph.unionElements({faceId, componentId}, {twinFaceId, twinComponentId});
    }
}

void ReebSpace::computeTraversal(const TetMesh &tetMesh, const Arrangement &arrangement, const bool discardFiberSeedsSets)
{
    Face_const_handle outerFace;
    for (Face_const_iterator fit = arrangement.arr.faces_begin(); fit != arrangement.arr.faces_end(); ++fit) 
    {
        if (fit->is_unbounded()) 
        {
            assert(outerFace == Face_const_handle());
            outerFace = fit;
            break;
        }
    }

    // Sanity check, make sure the outer face is simple
    assert(outerFace->number_of_inner_ccbs() == 1);
    assert(outerFace->number_of_outer_ccbs() == 0);
    assert(outerFace->number_of_holes() == 1);
    assert(outerFace->number_of_isolated_vertices() == 0);

    std::queue<Face_const_handle> traversalQueue;
    // This is the order in which faces are added, also servers as a visited flag
    std::vector<int> order(arrangement.arrangementFacesIdices.size(), -1);


    traversalQueue.push(outerFace);
    // This is the order in which a face has been processed, note this is different than level
    // This is used as an index for faces, so that we don't do double work when checking edges for correspondence
    int orderIndex = 0;
    order[arrangement.arrangementFacesIdices.at(outerFace)] = orderIndex;

    this->preimageGraphs.resize(arrangement.arrangementFacesIdices.size());

    // If we want the fiber seeds, initialize them
    if (false == discardFiberSeedsSets)
    {
        this->fiberSeeds.resize(arrangement.arrangementFacesIdices.size());
    }

    auto bar = Timer::getLoadingBar();
    int graphsInMemory = 0;
    float averageAraphsInMemory = 0;
    int barTickThreshold = arrangement.arrangementFacesIdices.size() / 50;

    while (false == traversalQueue.empty())
    {
        Face_const_handle currentFace = traversalQueue.front();
        traversalQueue.pop();

        int currentFaceID = arrangement.arrangementFacesIdices.at(currentFace);

        // Iterate over all neighbouring faces
        Halfedge_const_handle startingHalfedge = currentFace->is_unbounded() ? *currentFace->holes_begin() : *currentFace->outer_ccbs_begin();

        Arrangement_2::Ccb_halfedge_const_circulator currentHalfedge = startingHalfedge;
        do {

            Face_const_handle twinFace = currentHalfedge->twin()->face();
            const int twinFaceID = arrangement.arrangementFacesIdices.at(twinFace);

            // If the neighbour has not been visited, we queue it and compute its preimage graph
            if (-1 == order[twinFaceID])
            //if (-1 == order[twinFaceID] && arrangement.singularFaces.contains(twinFace))
            {
                traversalQueue.push(twinFace);
                order[twinFaceID] = ++orderIndex;

                ReebSpace::computeTwinFacePreimageGraph(tetMesh, arrangement, currentHalfedge);
                graphsInMemory++;

                // Initialize the vertices of H with the connected components of the twin graph
                const std::vector<std::pair<int, int>> componentRepresentatives = this->preimageGraphs[twinFaceID].getComponentRepresentatives();

                for (const auto &[componentId, triangleId] : componentRepresentatives)
                {
                    this->correspondenceGraph.addElement({twinFaceID, componentId});
                }

                if (false == discardFiberSeedsSets)
                {
                    this->fiberSeeds[twinFaceID] = componentRepresentatives;
                }
            }

            if (order[currentFaceID] < order[twinFaceID])
            {
                ReebSpace::determineCorrespondence(tetMesh, arrangement, currentHalfedge);
            }

            ++currentHalfedge;

        } while (currentHalfedge != startingHalfedge);

        averageAraphsInMemory = averageAraphsInMemory + ((float)graphsInMemory - (float)averageAraphsInMemory) / (float)orderIndex;

        // Update the loading bar
        if (barTickThreshold > 0 && order[currentFaceID] % barTickThreshold == 0)
        {
            if (false == bar->is_completed()) {  bar->tick(); }
            if (false == bar->is_completed()) {  bar->tick(); }
        }

        // Dispose of the preimage graph we will no longer need it
        //this->preimageGraphs[currentFaceID].clear();
        //graphsInMemory--;
    }

    this->correspondenceGraph.finalise();

    bar->set_progress(100); // all done
    printf("\n\nThere is an average of %f / %ld active preimage graphs.\n", averageAraphsInMemory, this->preimageGraphs.size());
    printf("The correspondence graphs has %ld vertices and the Reeb space has %ld sheets.\n\n", this->correspondenceGraph.data.size(), this->correspondenceGraph.getComponentRepresentatives().size());
}

void ReebSpace::computeSheetGeometry(const TetMesh &tetMesh, const Arrangement &arrangement)
{
    // For faster lookups cache the sheetd IDs for each face
    std::vector<std::unordered_set<int>> faceSheets(this->fiberSeeds.size());
    for (int i = 0 ; i < this->fiberSeeds.size() ; i++)
    {
        for (const auto &[fiberComponentId, triangleId] : this->fiberSeeds[i])
        {
            const int sheetId = this->correspondenceGraph.findElement({i, fiberComponentId});
            faceSheets[i].insert(sheetId);
        }

    }

    Timer::start();
    // In order to compute the polygon of each sheet, first obtain a halfEdge of the arrangement that is on the boundary of the sheet
    std::unordered_map<int, Arrangement_2::Halfedge_const_handle> sheetSeeds;

    for (auto currentFaceIterator = arrangement.arr.faces_begin(); currentFaceIterator != arrangement.arr.faces_end(); ++currentFaceIterator) 
    {
        Arrangement_2::Face_const_handle currentFace = currentFaceIterator;
        if (currentFace->is_unbounded()) { continue; }

        const int currentFaceID = arrangement.arrangementFacesIdices.at(currentFace);

        // Which sheets contain the current face
        const std::unordered_set<int> currentFaceSheetIds = faceSheets[currentFaceID];


        Arrangement_2::Ccb_halfedge_const_circulator start = currentFace->outer_ccb();
        Arrangement_2::Ccb_halfedge_const_circulator curr = start;
        do {
            Face_const_handle twinFace = curr->twin()->face();
            if (twinFace->is_unbounded()) { curr++; continue; }

            const int twinFaceID = arrangement.arrangementFacesIdices.at(twinFace);
            const std::unordered_set<int> twinFaceSheetIds = faceSheets[twinFaceID];

            // Which sheets are in the currentFace, but NOT in the twin face
            std::vector<int> diff;
            for (const int &sheetId : currentFaceSheetIds) 
            {
                if (false == twinFaceSheetIds.contains(sheetId))
                {
                    diff.push_back(sheetId);
                }
            }

            for (const int &sheetId : diff)
            {
                //std::cout << "Half-edge   from: " << curr->source()->point() << " to " << curr->target()->point() << " is a seed for sheet " << sheetId << std::endl;

                if (false == sheetSeeds.contains(sheetId))
                {
                    sheetSeeds[sheetId] = curr;
                }
            }

            ++curr;
        } while (curr != start);
    }
    Timer::stop("Computing sheet seeds                  :");





    //
    // At this point we have a half edge seed for every sheet, we use it to traverse the boundary of the sheet to obtain its geometry
    // We rely on the half edge data structure provided by the arrangement
    // if ab is an edge between vertices a and b in A and, then we can get the next edge bc as the one of the neighours of b that has
    // the sheet in its face and it doesn't have the sheet in the twin face
    //


    //printf("Done with sheet seeds\n");


    Timer::start();
    for (const auto &[sheetId, halfEdge]: sheetSeeds)
    {
        // Sanity check, make sure the seed we have picked actually is on the boundary of the sheet
        Face_const_handle faceSeedEdge = halfEdge->face();
        const int faceSeedEdgeId = arrangement.arrangementFacesIdices.at(faceSeedEdge);
        if (false == faceSheets[faceSeedEdgeId].contains(sheetId))
        {
            throw std::runtime_error("The face of the seed edge does is not in the sheet.");
        }

        Face_const_handle faceSeedTwinEdge = halfEdge->twin()->face();
        const int faceSeedTwinEdgeId = arrangement.arrangementFacesIdices.at(faceSeedTwinEdge);
        if (true == faceSheets[faceSeedTwinEdgeId].contains(sheetId))
        {
            throw std::runtime_error("The face of the twin of the seed edge does is no in the sheet.");
        }


        //printf("---------------------------------------------------------------------------------------------       At sheet %d \n", sheetId);

        // Starting location, we don't make it visited on purpose, we want to discover it latex to finish the loop
        Vertex_const_handle startVertex = halfEdge->source();

        // Time to find the next vertex
        Vertex_const_handle currentVertex = startVertex;
        std::unordered_set<Vertex_const_handle> visited;
        do 
        {

            // Add the current vertex to the polygon
            this->sheetPolygon[sheetId].push_back({
                    CGAL::to_double(currentVertex->point().x()), 
                    CGAL::to_double(currentVertex->point().y())
                    });


            //printf("-------------------------- The current vertex (%f, %f) \n ", CGAL::to_double(currentVertex->point().x()), CGAL::to_double(currentVertex->point().y()));



            //
            // Find the next vertex on the boundary of the sheet, it will be adjacent to the current vertex
            //
            int sheetBoundaryNeighbours = 0;
            int allNeighbours = 0;
            visited.insert(currentVertex);

            // This loop should always finish
            const auto begin = currentVertex->incident_halfedges();
            auto circ = begin;
            do {

                allNeighbours++;

                // Too many neighourss
                if (allNeighbours > arrangement.arr.number_of_vertices())
                {
                    break;
                }

                const auto twinHalfEdge = circ->twin();
                const auto nextVertex = twinHalfEdge->target();
                //printf("Looking at neighbour (%f, %f) \n ", CGAL::to_double(twinHalfEdge->target()->point().x()), CGAL::to_double(twinHalfEdge->target()->point().y()));

                Face_const_handle faceA = circ->face();
                const int faceAId = arrangement.arrangementFacesIdices.at(faceA);

                Face_const_handle faceB = twinHalfEdge->face();
                const int faceBId = arrangement.arrangementFacesIdices.at(faceB);

                // If one of the face contains the sheet
                const bool faceHalfEdgeContainsSheet = faceSheets[faceAId].contains(sheetId);
                const bool faceTwinHalfEdgeContainsEdge = faceSheets[faceBId].contains(sheetId);

                if (faceHalfEdgeContainsSheet == true && faceTwinHalfEdgeContainsEdge == false)
                {
                    currentVertex = nextVertex;
                    sheetBoundaryNeighbours++;
                }

                ++circ;
            } while (circ != begin);

            // A the polygon of a sheet sheet cannot have more vertices than the number of the vertices in the arrangement, something went wrong.
            if (this->sheetPolygon[sheetId].size() > arrangement.arr.number_of_vertices())
            {
                printf("The boundary of sheet %d is degenerate, it has more vertices %ld than the arrangement.\n", sheetId, this->sheetPolygon[sheetId].size());
                this->incompleteSheets.insert(sheetId);
                break;
            }
            
            if (sheetBoundaryNeighbours !=1)
            {
                printf("The boundary of sheet %d is degenerate, a vertex has %d boundary neighours, should only be one.\n", sheetId, sheetBoundaryNeighbours);
                this->incompleteSheets.insert(sheetId);
                break;
            }

            // I'm not sure if this will ever get triggered, since if we loop back to someone else than the start 
            // that should mean 
            if (currentVertex != startVertex && visited.contains(currentVertex))
            {
                printf("The boundary of sheet %d is degenerate, a vertex loops back to an already visited one.\n", sheetId);
                this->incompleteSheets.insert(sheetId);
                break;

            }

            if (allNeighbours > arrangement.arr.number_of_vertices())
            {
                printf("The boundary of sheet %d is degenerate, the inner loop kept going for too long.\n", sheetId);
                this->incompleteSheets.insert(sheetId);
                break;
            }

        } while (currentVertex != startVertex);
    }
    Timer::stop("Computing sheet boundary polygons      :");

    if (this->incompleteSheets.size() == 0)
    {
        std::cout << "\nThe boundaries of all sheets are okay! No degeneracy.\n" << std::endl;
    }
}


void ReebSpace::computeSheetArea(const TetMesh &tetMesh, const Arrangement &arrangement)
{
    Timer::start();
    for (const auto &[sheetId, polygon] : this->sheetPolygon)
    {
        double area = 0.0;

        // If the sheet is incomplete sum up the faces that make it
        if (this->incompleteSheets.contains(sheetId))
        {
            // Loop through all faces to see which ones are in the sheet
            for (auto f = arrangement.arr.faces_begin(); f != arrangement.arr.faces_end(); ++f) 
            {
                const int currentFaceID = arrangement.arrangementFacesIdices.at(f);

                // For each fiber component in the face, see if one of those is in our sheet
                for (const auto &[fiberComponentId, triangleId] : this->fiberSeeds[currentFaceID])
                {
                    const int componentSheetId = this->correspondenceGraph.findElement({currentFaceID, fiberComponentId});

                    // Now we can add the polygon
                    if (componentSheetId == sheetId)
                    {
                        CartesianPolygon_2 poly;

                        typename Arrangement_2::Ccb_halfedge_const_circulator circ = f->outer_ccb();
                        typename Arrangement_2::Ccb_halfedge_const_circulator curr = circ;
                        do {
                            typename Arrangement_2::Halfedge_const_handle e = curr;

                            // Get point from CGAL (and convert to double )
                            const float u = CGAL::to_double(e->source()->point().x());
                            const float v = CGAL::to_double(e->source()->point().y());

                            poly.push_back({u, v});
                        } while (++curr != circ);

                        area += poly.area();
                    }
                }
            }
            printf("Incomplete sheet %d has area %f\n", sheetId, area);
        }
        else
        {
            area = abs(polygon.area());
        }

        this->sheetArea[sheetId] = area;
    }
    Timer::stop("Computing sheet areas                  :");
}


void ReebSpace::printTopSheets(const TetMesh &tetMesh, const Arrangement &arrangement, const int &numberToPrint)
{
    Timer::start();
    // Transfer the map entries into a vector of pairs so that we can sort
    std::vector<std::pair<int, double>> sheetAreaSortVector(this->sheetArea.begin(), this->sheetArea.end());

    //printf("The sort vector has size %ld\n", sheetAreaSortVector.size());

    // Sort by area
    std::sort(sheetAreaSortVector.begin(), sheetAreaSortVector.end(), [](const std::pair<int, double>& a, const std::pair<int, double>& b) {
            return a.second > b.second; 
            });

    // Assign a sorted order to each sheet
    for (int i = 0 ; i < sheetAreaSortVector.size() ; i++)
    {
        const auto &[sheetId, area] = sheetAreaSortVector[i];
        this->sheetConsequitiveIndices[sheetId] = i;
    }
    Timer::stop("Sorting sheets and labeling them       :");

    const int actualNumberToPrint = std::min((size_t)30, sheetAreaSortVector.size());

    std::cout << "\nHere are the top " << actualNumberToPrint << " sheets sorted by range area.\n";
    // Print to debug, at least the first new
    //for (int i = 0 ; i < sheetAreaSortVector.size() ; i++)
    for (int i = 0 ; i < actualNumberToPrint ; i++)
    {
        std::cout << i << " -- sheet " << sheetAreaSortVector[i].first << " has area " << sheetAreaSortVector[i].second << std::endl;
    }
}
