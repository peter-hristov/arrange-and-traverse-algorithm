#include "./Timer.h"
#include "./Arrangement.h"



void Arrangement::buildAABBtree(const TetMesh &tetMesh)
{
    for (int edgeId = 0 ; edgeId < tetMesh.edges.size() ; edgeId++)
    {
        const std::array<int, 2> &edge = tetMesh.edges[edgeId];
        const int &edgeType = tetMesh.edgeSingularTypes.at(edge);

        this->allSegments.push_back(Segment_2(this->arrangementPoints[edge[0]], this->arrangementPoints[edge[1]]));

    }

    this->tree = TreeAABB(this->allSegments.begin(), this->allSegments.end());
    this->tree.build();



    for (const auto &[edge, type] : tetMesh.edgeSingularTypes) 
    {
        if (type != 1)
        {
            this->singularSegments.push_back(Segment_2(this->arrangementPoints[edge[0]], this->arrangementPoints[edge[1]]));
        }
    }

    this->treeSingular = TreeAABB(this->singularSegments.begin(), this->singularSegments.end());
    this->treeSingular.build();

}

Face_const_handle Arrangement::getActiveFace(const Point_2 &query_point)
{
    // Locate the point in the arrangement
    CGAL::Object result = pl->locate(query_point);

    // Try to assign to a face, edge or a vertex
    Arrangement_2::Face_const_handle face;
    Arrangement_2::Halfedge_const_handle edge;
    Arrangement_2::Vertex_const_handle vertex;

    if (CGAL::assign(face, result)) 
    {
    } 
    // If we are on an edge, just grad an adjacent face
    else if (CGAL::assign(edge, result)) 
    {
        face = edge->face();
    } 
    // If we are on a vertex grab an indicent edge and get its face
    else if (CGAL::assign(vertex, result)) 
    {
        edge = vertex->incident_halfedges();
        face = edge->face();
    } else 
    {
        assert(false);
    }

    return face;
}

Face_const_handle Arrangement::getActiveFace(const std::array<double, 2> fiberPoint)
{
    // The query point (u, v)
    Point_2 query_point(fiberPoint[0], fiberPoint[1]);

    return this->getActiveFace(query_point);
}

void Arrangement::computeArrangement(const TetMesh &tetMesh, const SegmentMode &segmentMode) 
{
    // Add in the vertices of the mesh 
    this->arrangementPoints.resize(tetMesh.vertexCoordinatesF.size());
    for (int i = 0 ; i < tetMesh.vertexCoordinatesF.size() ; i++)
    {
        const double u = tetMesh.vertexCoordinatesF[i];
        const double v = tetMesh.vertexCoordinatesG[i];
        const Point_2 point(u, v);

        this->arrangementPoints[i] = point;
        this->arrangementPointIndices[point] = i;

    };

    // Add the unique edges as setments to the arrangement
    std::vector<Segment_2> segments;
    //for (const auto& edge : tetMesh.edges) 
    for (const auto &[edge, type] : tetMesh.edgeSingularTypes) 
    {
        if (segmentMode == SegmentMode::UseSingularSegments)
        {
            if (type != 1)
            {
                segments.push_back(Segment_2(this->arrangementPoints[edge[0]], this->arrangementPoints[edge[1]]));
            }

        }
        else
        {
            segments.push_back(Segment_2(this->arrangementPoints[edge[0]], this->arrangementPoints[edge[1]]));
        }
    }

    CGAL::insert(this->arr, segments.begin(), segments.end());


    // Sequential arrangement computationa
    //Timer::start();
    //Arrangement_2 arr;
    //for (const auto& segment : segments) {
        //CGAL::insert(arr, Curve_2(segment));
    //}
    //Timer::stop("Computed arrangement sequantially      :");

    std::cout << "The arrangement size:"
        << "   |V| = " << this->arr.number_of_vertices()
        << ",  |E| = " << this->arr.number_of_edges()
        << ",  |F| = " << this->arr.number_of_faces() << std::endl << std::endl;

    // Set up the indices and their reverse lookup for all faces and vertices
    int counter = 0;
    this->arrangementIndexToFace.resize(this->arr.number_of_faces());
    for (auto f = this->arr.faces_begin(); f != this->arr.faces_end(); ++f) 
    {
        Arrangement_2::Face_const_handle a = f;
        this->arrangementFacesIdices[a] = counter;
        this->arrangementIndexToFace[counter] = f;
        counter++;
    }

    for (auto vit = this->arr.vertices_begin(); vit != this->arr.vertices_end(); ++vit)
    {
        Vertex_const_handle vertexHandle = vit;
        arrangementPointHandles[vertexHandle->point()] = vertexHandle;
    }
}

void Arrangement::computePointLocationDataStructure()
{
    this->pl = std::make_unique<Point_location>(this->arr);
}

void Arrangement::connectNestedFaces(TetMesh &tetMesh)
{
    int faceBoundarySizeTotal = 0;
    int innerFaces = 0;
    for (Face_const_iterator fit = this->arr.faces_begin(); fit != this->arr.faces_end(); ++fit) 
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
                //throw std::runtime_error("Outer face is degenerate!");
            }
        }
        else
        {
            std::vector<int> ounterBoundaryVertices;

            int faceBoundarySize = 0;
            auto circ = fit->outer_ccb();
            auto start = circ;
            do
            {
                if (this->arrangementPointIndices.contains(circ->source()->point()))
                {
                    ounterBoundaryVertices.push_back(
                            this->arrangementPointIndices.at(circ->source()->point())
                            );

                }

                ++circ;
                faceBoundarySize++;
            } while (circ != start);

            innerFaces++;
            faceBoundarySizeTotal += faceBoundarySize;

            if (
                    fit->number_of_inner_ccbs() != 0 ||
                    fit->number_of_outer_ccbs() != 1 ||
                    fit->number_of_holes() != 0 ||
                    fit->number_of_isolated_vertices() != 0
               )
            {

                std::vector<std::set<int>> innerBoundaries;

                int innerFaceBoundarySize = 0;
                for (auto icit = fit->inner_ccbs_begin(); icit != fit->inner_ccbs_end(); ++icit) {
                    auto circ = *icit;   // circulator around this hole
                    auto start = circ;
                    std::set<int> innerBoundary;
                    do {
                        if (this->arrangementPointIndices.contains(circ->source()->point()))
                        {
                            innerBoundary.insert(
                                    this->arrangementPointIndices.at(circ->source()->point())
                                    );
                        }

                        ++circ;
                        innerFaceBoundarySize++;
                    } while (circ != start);

                    innerBoundaries.push_back(innerBoundary);
                }

                // count inner face boundary

                std::cerr << "\nFound a face with holes...\n";
                std::cerr << "number_of_inner_ccbs: " << fit->number_of_inner_ccbs() << std::endl;
                std::cerr << "number_of_outer_ccbs: " << fit->number_of_outer_ccbs() << std::endl;
                std::cerr << "number_of_holes: " << fit->number_of_holes() << std::endl;
                std::cerr << "number_of_isolated_vertices: " << fit->number_of_isolated_vertices() << std::endl;
                std::cerr << "Size of the outer boundary: " << faceBoundarySize << std::endl;
                std::cerr << "Size of the inner boundary: " << innerFaceBoundarySize << std::endl;
                std::cerr << "Is fictitious?: " << fit->is_fictitious() << std::endl;
                //throw std::runtime_error("An inner face is degenerate!");


                for (const auto &innerBoundary : innerBoundaries)
                {
                    const std::vector<int> shortestPath = tetMesh.findShortestPath(ounterBoundaryVertices, innerBoundary);
                    tetMesh.markPseudoSingularEdges(shortestPath);
                }
            }
        }
    }

    double averageFaceBoundarySize = (double)faceBoundarySizeTotal / (double)innerFaces;
    std::cerr << "There is this number of average faces in a boundary " << averageFaceBoundarySize << std::endl;



    for (auto he = this->arr.halfedges_begin(); he != this->arr.halfedges_end(); ++he)
    {
        // If you want to not iterate over the same edge twice
        if (he->face() == he->twin()->face())
        {
            std::cerr << "\nA half-edge has the same face as it's twin!!!\n";
        }
    }
}



void Arrangement::assignIndices()
{
    int faceId = 0;
    std::set<int> uniqueFaceIds;
    for (auto fit = arr.faces_begin(); fit != arr.faces_end(); ++fit, ++faceId) 
    {
        fit->set_data(faceId);           // store ID in the face
        uniqueFaceIds.insert(faceId);
    }

    if (uniqueFaceIds.size() != arr.number_of_faces())
    {
        throw std::runtime_error("Not all faces have unique IDs!");
    }

    int hedgeId = 0;
    std::set<int> uniqueHalfEdgeIds;
    for (auto heit = arr.halfedges_begin(); heit != arr.halfedges_end(); ++heit, ++hedgeId) 
    {
        heit->set_data({hedgeId, false});
        uniqueHalfEdgeIds.insert(hedgeId);
    }

    if (uniqueHalfEdgeIds.size() != arr.number_of_halfedges())
    {
        throw std::runtime_error("Not all faces have unique IDs!");
    }
}


void Arrangement::assignHalfEdgePseudoSingular(const TetMesh &tetMesh, Arrangement &singularArrangement)
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
            halfEdge->data().isPseudoSingular = true;
        }
        else
        {
            halfEdge->data().isPseudoSingular = false;
        }
    }

}

// Given a point p within a face F, find a vertex v in the polygon such that the segment pv is entirely within F
// This relies on that fact given a point in a simple, non-intersecting polygon there at least one visible vertex from p.
//
// In practise we use the number of intersected halfEdge, which sould be zero, and only one vertex should be intersected.
//
Point_2 Arrangement::findVisibleVertex(const Face_const_handle activeFace, const Point_2 p)
{
    auto circ = activeFace->outer_ccb();
    auto start = circ;
    do
    {
        Arrangement_2::X_monotone_curve_2 segmentMonotoneCurve(
                p,
                circ->source()->point());

            std::vector<Arrangement_2::Vertex_handle> vertices;
            std::vector<Arrangement_2::Halfedge_handle> halfEdges;

            CGAL::zone(
                    this->arr, 
                    segmentMonotoneCurve, 
                    CGAL::dispatch_or_drop_output<Arrangement_2::Vertex_handle, Arrangement_2::Halfedge_handle>(
                        std::back_inserter(vertices), 
                        std::back_inserter(halfEdges))
                    );

            if (halfEdges.size() == 0 && vertices.size() == 1)
            {
                return circ->source()->point();
            }


        ++circ;
    } while (circ != start);

    throw std::runtime_error("No arrangement vertices are visible from the control point.");
}




std::vector<std::tuple<K::FT, int, int>> Arrangement::getIntersectedSegments2(TetMesh &tetMesh, const Segment_2 &controlSegment, const bool shouldSort)
{
    //Timer::start();
    std::vector<TreeAABB::Primitive_id> intersectedSegmentsAABB;
    this->tree.all_intersected_primitives(controlSegment, std::back_inserter(intersectedSegmentsAABB));
    //Timer::stop("Computing AABB                         :");

    std::vector<std::tuple<K::FT, int, int>> intersectedSegments;
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
        const int segmentIndex = id - this->allSegments.begin();

        const int type = tetMesh.edgeSingularTypes.at(tetMesh.edges.at(segmentIndex));

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

        intersectedSegments.emplace_back(alpha, segmentIndex, type);
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


    // Only keep the part until a singular segment
    //
    //int counter = 0;

    //for (const auto &[alpha, edgeId, edgeType] : intersectedSegments)
    //{
        //if (edgeType == 2 || edgeType == 0)
        //{
            //break;
        //}
        //counter++;
    //}

    //intersectedSegments.erase(intersectedSegments.begin() + counter + 1, intersectedSegments.end());

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
