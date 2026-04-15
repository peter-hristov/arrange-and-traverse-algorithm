#include "FiberSurface.h"
#include "FiberLabeling.h"
#include "io.h"

void FiberSurface::print()
{
    const auto edgeParamMap = this->edgeParam();
    const auto tetIdMap = this->tetId();
    const auto sheetIdMap = this->sheetId();

    std::cout << "Vertices:\n";
    for (auto v : mesh.vertices())
    {
        const auto &p = mesh.point(v);
        double e = edgeParamMap[v];
        std::cout << "  Vertex " << v << ": ("
            << p[0] << ", " << p[1] << ", " << p[2]
            << "), edgeParam = " << e << "\n";
    }

    std::cout << "\nFaces:\n";
    for (auto f : mesh.faces())
    {
        std::cout << "  Face " << f << ": [";
        bool first = true;
        for (auto v : vertices_around_face(mesh.halfedge(f), mesh))
        {
            if (!first) std::cout << ", ";
            std::cout << v;
            first = false;
        }
        std::cout << "], tetId = " << tetIdMap[f]
            << ", sheetId = " << sheetIdMap[f] << "\n";
    }
}

void FiberSurface::printSheetHistogram(ReebSpace2 &reebSpace)
{
    const auto sheetIdMap = this->sheetId();

    std::set<int> sheetIds;
    for (auto f : mesh.faces())
    {
        sheetIds.insert(sheetIdMap[f]);
    }

    std::vector<std::tuple<double, double, int>> sheetsAndAreas;
    sheetsAndAreas.reserve(sheetIds.size());

    for (const int id : sheetIds)
    {
        const std::tuple<double, double, int> sheetAndArea = {reebSpace.sheetArea.at(id), reebSpace.sheetAreaProportion.at(id), id};
        sheetsAndAreas.emplace_back(sheetAndArea);
    }

    std::sort(sheetsAndAreas.begin(), sheetsAndAreas.end(), std::greater<>());

    std::cout << "\nThe fiber surface has the following histogram of intersected sheets.\n";
    printf("%-10s %-8s %-12s %-8s %-12s\n", "Sheet ID", "|", "Area", "|", "Percentage");
    printf("----------------------------------------------------\n");
    for (const auto &[area, areaProportion, id] : sheetsAndAreas)
    {
        printf("%-10d %-8s %-12.4f %-8s %.4f%%\n", id, "|", area, "|", areaProportion);
    }
}

void FiberSurface::remesh(const std::vector<double> &isovalues)
{
    for (int i = 0 ; i < isovalues.size() ; i++)
    {
        this->remeshOnce(isovalues[i]);
    }

    this->triangulate();
}

bool FiberSurface::bfsComponentFromSeed(const TetMesh &tetMesh, const Arrangement &singularArrangement, const int seedTriangleId, const std::vector<int> &tetTriangleIds, const CartesianPoint &controlPoint, std::vector<bool> &visited)
{
    std::queue<int> bfsQueue;
    bfsQueue.push(seedTriangleId);
    visited[seedTriangleId] = true;

    while (!bfsQueue.empty())
    {
        const int currentTriangleId = bfsQueue.front();
        bfsQueue.pop();

        // If it's the one we want, we are done
        if (std::find(tetTriangleIds.begin(), tetTriangleIds.end(), currentTriangleId) != tetTriangleIds.end())
        {
            return true;
        }

        for (const int &neighbourTriangleId : tetMesh.tetIncidentTriangles[currentTriangleId])
        {
            // Skip if it's visited
            if (visited[neighbourTriangleId]) { continue; }

            // Skip if it's not active
            if (false == tetMesh.isTriangleActive(neighbourTriangleId, controlPoint)) { continue; }

            bfsQueue.push(neighbourTriangleId);
            visited[neighbourTriangleId] = true;
        }
    }

    return false;
}

int FiberSurface::findFiberPointComponent(const TetMesh &tetMesh, const Arrangement &singularArrangement, const std::vector<std::pair<int, int>> &fiberSeeds, const std::vector<int> &tetTriangleIds, const Segment_2 &controlSegment, const double pointAlpha)
{
    const Point_2 controlPoint = CGAL::barycenter(controlSegment[0], 1.0 - pointAlpha, controlSegment[1], pointAlpha);
    const CartesianPoint controlPointCartesian(CGAL::to_double(controlPoint.x()), CGAL::to_double(controlPoint.y()));

    std::vector<bool> visited(tetMesh.triangleIndices.size(), false);

    for (const auto &[triangleId, componentId] : fiberSeeds)
    {
        if (bfsComponentFromSeed(tetMesh, singularArrangement, triangleId, tetTriangleIds, controlPointCartesian, visited))
        {
            return componentId;
        }
    }

    return -1;
}










int FiberSurface::labelTriangle(TetMesh &tetMesh, Arrangement &singularArrangement, ReebSpace2 &reebSpace, const CGALMesh::Face_index &triangle, const std::vector<std::tuple<K::FT, int, int>> &intersectedSegments, const Segment_2 &controlSegment)
{
    // 1. Compute the edgePara at the center of the triangle
    double midPointAlpha = 0.0;
    for (const auto v : mesh.vertices_around_face(mesh.halfedge(triangle))) 
    {
        midPointAlpha += this->edgeParam()[v];
    }
    midPointAlpha /= 3.0;

    // 2. Compute the fiber graph at the alpha in the range
    const std::vector<std::pair<int, int>> fiberSeeds = fiber::labeling::computeFiberSeedsGivenLine(tetMesh, singularArrangement, reebSpace, controlSegment, midPointAlpha, intersectedSegments);

    // 3. Determine which fiber component contains a triangle from the tet
    const int tetId = this->tetId()[triangle];

    const int a = tetMesh.tetrahedra[tetId][0];
    const int b = tetMesh.tetrahedra[tetId][1];
    const int c = tetMesh.tetrahedra[tetId][2];
    const int d = tetMesh.tetrahedra[tetId][3];

    const std::vector<int> tetTriangleIds = {
        tetMesh.triangleIndices.at({a, b, c}),
        tetMesh.triangleIndices.at({a, b, d}),
        tetMesh.triangleIndices.at({a, c, d}),
        tetMesh.triangleIndices.at({b, c, d}),
    };

    const int componentId = findFiberPointComponent(tetMesh, singularArrangement, fiberSeeds, tetTriangleIds, controlSegment, midPointAlpha);

    if (componentId != -1)
    {
        return reebSpace.correspondenceGraphDS.find(componentId);
    }

    std::cerr << "Triangle ID has not been found in the interactive fiber surface computation!\n";

    return -1;
}



void FiberSurface::labelFiberSurface(TetMesh &tetMesh, Arrangement &singularArrangement, ReebSpace2 &reebSpace, const std::vector<std::tuple<K::FT, int, int>> &intersectedSegments, const Segment_2 &controlSegment)
{
    auto sheetIdMap = this->sheetId();
    auto componentIdMap = this->componentId();
    auto isImpassableMap = this->isImpassable();

    std::size_t numComponents2 = CGAL::Polygon_mesh_processing::connected_components(
            this->mesh,
            componentIdMap,
            CGAL::parameters::edge_is_constrained_map(isImpassableMap)
            );

    // 2. Get a number of representative triangles per component (in case just one doesn't work)
    //
    std::vector<std::vector<CGALMesh::Face_index>> componentRepresentatives(numComponents2);

    for (auto face : mesh.faces())
    {
        const int componentId = componentIdMap[face];

        if (componentRepresentatives[componentId].size() < 100)
        {
            componentRepresentatives[componentId].push_back(face);
        }
    }

    // 3. Compute one flexible fiber per representative triangle
    //
    std::vector<int> componentSheets(componentRepresentatives.size(), -1);

#pragma omp parallel for schedule(dynamic)
    for (int componentId = 0 ; componentId < componentRepresentatives.size() ; componentId++)
    {
        if (componentRepresentatives[componentId].empty())
        {
            //std::cerr << "Component " << componentId << " has not been assigned a representative triangle.\n";
        }

        for (CGALMesh::Face_index face : componentRepresentatives[componentId])
        {
            try
            {
                componentSheets[componentId] = this->labelTriangle(tetMesh, singularArrangement, reebSpace, face, intersectedSegments, controlSegment);

                if (componentSheets[componentId] != -1)
                {
                    break;
                }
            }
            catch (...) { }

            //std::cerr << "Could not compute fiber fabeling for representative triangle trying another ...\n";
        }

        if (componentSheets[componentId] == -1)
        {
            //std::cerr << "Could not label component " << componentId << " ...\n";
        }
    }


    // 4. Set up the sheetIds of each triangle based on the connected component
    //
    for (auto face : mesh.faces())
    {
        const int componentId = componentIdMap[face];
        sheetIdMap[face] =  componentSheets[componentId];
    }
}



CartesianPoint_3 FiberSurface::interpolateVertex(CGALMesh::Vertex_index v0, CGALMesh::Vertex_index v1, double isovalue)
{
    const auto edgeParamMap = this->edgeParam();

    const double e0 = edgeParamMap[v0];
    const double e1 = edgeParamMap[v1];

    const double t = (isovalue - e0) / (e1 - e0);

    const CartesianPoint_3 &p0 = mesh.point(v0);
    const CartesianPoint_3 &p1 = mesh.point(v1);

    //return CartesianPoint_3(p0 + t * (p1 - p0));
    return CartesianPoint_3(
            (1.0 - t) * p0.x() + t * p1.x(),
            (1.0 - t) * p0.y() + t * p1.y(),
            (1.0 - t) * p0.z() + t * p1.z()
            );
}










std::pair<CGALMesh::Property_map<CGALMesh::Vertex_index, int>, std::vector<CGALMesh::Vertex_index>> FiberSurface::getVertexColours(CGALMesh &cgalMesh, const double isovalue)
{
    const auto edgeParamMap = this->edgeParam();

    auto [vertexColour, created] = cgalMesh.add_property_map<CGALMesh::Vertex_index, int>("v:colour", -2);
    if (!created)
        throw std::runtime_error("Could not make mesh colour array.");

    std::vector<CGALMesh::Vertex_index> grayVertices;

    for (auto v : cgalMesh.vertices())
    {
        const double e = edgeParamMap[v];
        const double diff = e - isovalue;

        if (std::abs(diff) <= FiberSurface::epsilon)
        {
            vertexColour[v] = 0;       // on the isovalue
            grayVertices.push_back(v);
        }
        else if (diff < 0.0)
            vertexColour[v] = -1;      // below
        else
            vertexColour[v] = +1;      // above
    }

    return {vertexColour, grayVertices};
}

void FiberSurface::remeshOnce(const double isovalue)
{

    auto edgeParamMap = this->edgeParam();
    auto tetIdMap = this->tetId();
    auto sheetIdMap = this->sheetId();
    auto componentIdMap = this->componentId();
    auto isImpassableMap = this->isImpassable();

    // 1. Compute the colours of all the vertices
    //
    auto [vertexColour, grayVertices] = getVertexColours(this->mesh, isovalue);


    // 2. Fill in the active edges
    //
    std::unordered_set<CGALMesh::Edge_index> activeEdges;

    // If no gray vertices were found, just brute force it, something went wrong with the numerics
    if (grayVertices.size() == 0)
    {
        for (const auto &e : this->mesh.edges())
        {
            const auto h = this->mesh.halfedge(e);
            const auto v0 = this->mesh.source(h);
            const auto v1 = this->mesh.target(h);
            if (vertexColour[v0] * vertexColour[v1] == -1)
            {
                activeEdges.insert(e);
            }
        }
    }
    // If there were gray vertices search only from them to fill in the contour
    else
    {
        activeEdges = getActiveEdges(grayVertices, vertexColour);
    }


    std::unordered_set<CGALMesh::Face_index> facesToSplit;

    // 3. Split the edges
    //
    for (const auto &e : activeEdges)
    {
        auto h = this->mesh.halfedge(e);

        const auto v0 = this->mesh.source(h);
        const auto v1 = this->mesh.target(h);
        const double val0 = edgeParamMap[v0];
        const double val1 = edgeParamMap[v1];

        const CartesianPoint_3 edgeVertex = interpolateVertex(v0, v1, isovalue);
        const auto hNew = CGAL::Euler::split_edge(h, mesh);
        const auto vNew = mesh.target(hNew);

        // Copy over the vertex values for the new points
        mesh.point(vNew) = edgeVertex;
        edgeParamMap[vNew] = isovalue;
        vertexColour[vNew] = 0;

        // Make whether this is impassable
        const auto isEdgeImpassable = isImpassableMap[e];
        isImpassableMap[mesh.edge(hNew)] = isEdgeImpassable;
        isImpassableMap[mesh.edge(mesh.next(hNew))] = isEdgeImpassable;

        // Save the adjacent faces that now need splitting
        if (mesh.face(hNew) != CGALMesh::null_face())
        {
            facesToSplit.insert(mesh.face(hNew));
        }
        if (mesh.face(mesh.opposite(hNew)) != CGALMesh::null_face())
        {
            facesToSplit.insert(mesh.face(mesh.opposite(hNew)));
        }
    }

    // 4. Split faces and save theones that need to be triangulated later
    //
    std::unordered_set<CGALMesh::Face_index> facesToTriangulate;

    for (const auto f : facesToSplit)
    {
        const int faceTetId = tetIdMap[f];

        if (mesh.degree(f) == 3)
        {
            throw std::runtime_error("Face to split has no new edge points and has degree " + std::to_string(mesh.degree(f)));
        }

        // Collect the two half-edge where we want to split
        std::vector<CGALMesh::Halfedge_index> activeHalfEdges;
        for (auto h : halfedges_around_face(mesh.halfedge(f), mesh))
        {
            if (vertexColour[mesh.target(h)] == 0)
            {
                activeHalfEdges.push_back(h);
            }
        }

        if (activeHalfEdges.size() != 2)
        {
            throw std::runtime_error("There should be exactly 2 active half-edges per active face.");
        }

        // Split
        const auto h0 = activeHalfEdges[0];
        const auto h1 = activeHalfEdges[1];
        const auto newH = CGAL::Euler::split_face(h0, h1, this->mesh);

        // Make the new edge impassable
        isImpassableMap[mesh.edge(newH)] = true;

        // Write tet array
        const auto newF0 = mesh.face(newH);
        const auto newF1 = mesh.face(mesh.opposite(newH));

        tetIdMap[newF0] = faceTetId;
        tetIdMap[newF1] = faceTetId;

        // Save non-triangle faces to later triangulation
        if (mesh.degree(newF0) > 3)
        {
            facesToTriangulate.insert(newF0);
        }
        if (mesh.degree(newF1) > 3)
        {
            facesToTriangulate.insert(newF1);
        }
    }


    // 6. Cleanup
    //
    this->mesh.remove_property_map(vertexColour);
}


void FiberSurface::triangulate()
{
    auto tetIdMap = this->tetId();

    // This visitor sets the tetId to the triangulated faces to that of their parent
    struct Visitor : CGAL::Polygon_mesh_processing::Triangulate_faces::Default_visitor<CGALMesh>
    {
        int tetId_val;
        CGALMesh::Property_map<CGALMesh::Face_index, int>& tetId;
        Visitor(int tetId, CGALMesh::Property_map<CGALMesh::Face_index, int>& t) : tetId_val(tetId), tetId(t) {}
        void after_subface_created(CGALMesh::Face_index f) { tetId[f] = tetId_val; }
    }; 

    for (const auto f : mesh.faces())
    {
        const int faceTetId = tetIdMap[f];

        if (mesh.degree(f) == 3)
        {
            continue;
        }

        Visitor visitor(faceTetId, tetIdMap);

        CGAL::Polygon_mesh_processing::triangulate_face(f, mesh, CGAL::parameters::visitor(visitor));
    }
}

void FiberSurface::repairMesh()
{
    // Finish up with some postprocessing
    mesh.collect_garbage(); // before calling connected_components
    CGAL::Polygon_mesh_processing::orient(this->mesh);

    // 1. Remove degenerate faces (zero or near-zero area)
    CGAL::Polygon_mesh_processing::remove_degenerate_faces(mesh);

    // 2. Remove degenerate edges (edges shorter than a threshold)
    CGAL::Polygon_mesh_processing::remove_degenerate_edges(mesh);
}

void FiberSurface::validateMesh()
{
    // Make sure the edge is valid
    if (false == CGAL::is_valid_polygon_mesh(this->mesh))
    {
        throw std::runtime_error("New mesh is not valid.");
    }

    for (auto f : this->mesh.faces())
    {
        if (mesh.degree(f) != 3)
        {
            std::cerr << "A face of the mesh is not a triangle, it has degree " + std::to_string(mesh.degree(f)) + ".";

            //throw std::runtime_error("A face of the mesh is not a triangle, it has degree " + std::to_string(mesh.degree(f)) + ".");
        }

        // check if face area is zero (or near zero)
        const auto area = CGAL::Polygon_mesh_processing::face_area(f, mesh);

        if (area < 1e-12)
        {
            //std::cerr << "\n\nFace " << f.idx() << " of the mesh has a tiny area:" << std::setprecision(17) << area << " \n";

            //// Print edge lengths
            //auto h = mesh.halfedge(f);
            //for (int i = 0; i < 3; ++i, h = mesh.next(h))
            //{
            //const auto p0 = mesh.point(mesh.source(h));
            //const auto p1 = mesh.point(mesh.target(h));
            //const double lenSq = CGAL::to_double(CGAL::squared_distance(p0, p1));
            //std::cerr << "Face " << f.idx() << " edge " << i << " length: " << std::setprecision(17) << lenSq << "\n";
            //}


            //throw std::runtime_error("A face of the mesh has a tiny area:" + std::to_string(area) + ".");
        }
    }
}

std::unordered_set<CGALMesh::Edge_index> FiberSurface::getActiveEdges(const std::vector<CGALMesh::Vertex_index> &grayVertices, const CGALMesh::Property_map<CGALMesh::Vertex_index, int> &vertexColour)
{
    std::unordered_set<CGALMesh::Edge_index> activeEdges;

    auto [visited, created] = this->mesh.add_property_map<CGALMesh::Face_index, bool>("v:visited", false);

    if (false == created)
    {
        throw std::runtime_error("Could not ad visited array.");
    }

    std::queue<CGALMesh::Face_index> queue;

    for (const auto &v : grayVertices)
    {
        // Double check we have nice vertices
        if (!mesh.halfedge(v).is_valid()) 
        {
            throw std::runtime_error("Initial vertex not valud.");
        }
        if (mesh.halfedge(v) == CGALMesh::null_halfedge()) 
        {
            throw std::runtime_error("Initial vertex not valud.");
        }

        for (auto f : faces_around_target(mesh.halfedge(v), mesh))
        {
            if (f == CGALMesh::null_face()) { continue; }
            visited[f] = true;
            queue.push(f);
        }
    }

    while (!queue.empty())
    {
        auto f = queue.front();
        queue.pop();

        for (auto h : halfedges_around_face(mesh.halfedge(f), mesh))
        {
            // If the edge is not intersected, skip it
            const auto v0 = this->mesh.source(h);
            const auto v1 = this->mesh.target(h);
            if (vertexColour[v0] * vertexColour[v1] != -1) continue;

            // Add the active edge regardless
            activeEdges.insert(mesh.edge(h));

            // If it's on the boundary skip it
            auto opp = mesh.opposite(h);
            if (mesh.is_border(opp)) continue;

            // Only queue the neighbor if not yet visited
            auto neighbor = mesh.face(opp);
            if (visited[neighbor]) continue;

            visited[neighbor] = true;
            queue.push(neighbor);
        }
    }

    // At the end, clean up:
    this->mesh.remove_property_map(visited);

    return activeEdges;
}


void FiberSurface::filterTriangles(const std::set<int> &selectedSheets)
{
    if (selectedSheets.empty()) { return; }

    const auto sheetIdMap = this->sheetId();

    std::vector<CGALMesh::Face_index> facesToRemove;
    for (const auto f : mesh.faces())
    {
        const int sheetId = sheetIdMap[f];

        if (false == selectedSheets.contains(sheetId))
        {
            facesToRemove.push_back(f);
        }
    }

    for (auto fd : facesToRemove)
    {
        mesh.remove_face(fd);
    }

    mesh.collect_garbage();
}


FiberSurface::FiberSurface(const std::vector<std::array<double, 3>> &vertexCoordinates, const std::vector<std::array<int, 3>> &triangles, const std::vector<double> &vertexEdgePara, const std::vector<int> &tetId)
{
    // Unpack the points
    std::vector<CartesianPoint_3> points;
    for (auto& p : vertexCoordinates)
    {
        points.push_back(CartesianPoint_3(p[0],p[1],p[2]));
    }

    auto polygons = triangles;

    // We assume that the mesh is already cleaned up, otherwuse this can be used, but its doesn't work now, fix if you need it, merging points makes issues with adding the maps to the vertices and faces
    //
    // Merge duplicate vertices
    //std::vector<std::size_t> old_to_new;
    //CGAL::Polygon_mesh_processing::merge_duplicate_points_in_polygon_soup(points, polygons,
    //CGAL::parameters::vertex_to_vertex_map(boost::make_iterator_property_map(
    //old_to_new.begin(), boost::identity_property_map(), std::size_t(0)
    //))
    //);


    //// Merge duplicate polygons
    //std::vector<std::size_t> old_to_new2;
    //CGAL::Polygon_mesh_processing::merge_duplicate_polygons_in_polygon_soup(points, polygons,
    //CGAL::parameters::vertex_to_vertex_map(boost::make_iterator_property_map(
    //old_to_new2.begin(), boost::identity_property_map(), std::size_t(0)
    //))
    //);


    // Orient triangles
    CGAL::Polygon_mesh_processing::orient_polygon_soup(points, polygons);

    // 3. Build Surface_mesh
    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, polygons, this->mesh);

    auto resultEdgeParam = this->mesh.add_property_map<CGALMesh::Vertex_index,double>(this->EDGE_PARAM_KEY, -1.0);

    if (false == resultEdgeParam.second)
    {
        throw std::runtime_error("EdgeParam property could not be added to the mesh.");
    }

    for (std::size_t i = 0; i < vertexEdgePara.size(); ++i)
    {
        resultEdgeParam.first[CGALMesh::Vertex_index(i)] = vertexEdgePara[i];
    }

    auto resultTetId = this->mesh.add_property_map<CGALMesh::Face_index, int>(this->TET_ID_KEY, -1);

    if (false == resultTetId.second)
    {
        throw std::runtime_error("TetId property could not be added to the mesh.");
    }

    for (std::size_t i = 0; i < tetId.size(); ++i)
    {
        resultTetId.first[CGALMesh::Face_index(i)] = tetId[i];
    }

    auto resultSheetId = this->mesh.add_property_map<CGALMesh::Face_index, int>(this->SHEET_ID_KEY, -1);

    if (false == resultSheetId.second)
    {
        throw std::runtime_error("SheetId property could not be added to the mesh.");
    }

    auto resultComponentId = this->mesh.add_property_map<CGALMesh::Face_index, int>(this->COMPONENT_ID_KEY, -1);

    if (false == resultComponentId.second)
    {
        throw std::runtime_error("ComponentId property could not be added to the mesh.");
    }

    auto resultIsImpassable = this->mesh.add_property_map<CGALMesh::Edge_index, bool>(IMPASSABLE_KEY, false);

    if (false == resultIsImpassable.second)
    {
        throw std::runtime_error("isImpassable property could not be added to the mesh.");
    }
}

FiberSurface::FiberSurface() 
{
    mesh.add_property_map<CGALMesh::Vertex_index, double>(EDGE_PARAM_KEY,    -1.0);
    mesh.add_property_map<CGALMesh::Face_index,   int>   (TET_ID_KEY,        -1);
    mesh.add_property_map<CGALMesh::Face_index,   int>   (SHEET_ID_KEY,      -1);
    mesh.add_property_map<CGALMesh::Face_index,   int>   (COMPONENT_ID_KEY,  -1);
    mesh.add_property_map<CGALMesh::Edge_index,   bool>  (IMPASSABLE_KEY,    false);
}



FiberSurface FiberSurface::constructSegmentedFiberSurface(TetMesh &tetMesh, Arrangement &singularArrangement, ReebSpace2 &reebSpace, const std::vector<std::array<double, 2>> &controlPoints, const std::set<int> &selectedSheets)
{
    const Point_2 startPoint(controlPoints[0][0], controlPoints[0][1]);
    const Point_2 endPoint(controlPoints[1][0], controlPoints[1][1]);
    const Segment_2 controlSegment(startPoint, endPoint);

    //Timer::start();
    const std::vector<std::tuple<K::FT, int, int>> intersectedSegments = singularArrangement.getIntersectedSegments2(tetMesh, controlSegment, true);
    //Timer::stop("Computed Alpha intersections           :");

    //Timer::start();
    FiberSurface surfaceMesh = io::computeFiberSurface(tetMesh.originalMesh, controlPoints[0][0], controlPoints[0][1], controlPoints[1][0], controlPoints[1][1]);
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


// Helpers to get the maps of the mesh
CGALMesh::Property_map<CGALMesh::Vertex_index, double> FiberSurface::edgeParam()
{
    auto r = mesh.property_map<CGALMesh::Vertex_index, double>(EDGE_PARAM_KEY);
    if (!r.has_value()) throw std::runtime_error("edgeParam not initialized");
    return r.value();
}

CGALMesh::Property_map<CGALMesh::Face_index,   int> FiberSurface::tetId()
{
    auto r = mesh.property_map<CGALMesh::Face_index, int>(TET_ID_KEY);
    if (!r.has_value()) throw std::runtime_error("tetId not initialized");
    return r.value();
}

CGALMesh::Property_map<CGALMesh::Face_index,   int> FiberSurface::sheetId()
{
    auto r = mesh.property_map<CGALMesh::Face_index, int>(SHEET_ID_KEY);
    if (!r.has_value()) throw std::runtime_error("sheetId not initialized");
    return r.value();
}

CGALMesh::Property_map<CGALMesh::Face_index,   int> FiberSurface::componentId()
{
    auto r = mesh.property_map<CGALMesh::Face_index, int>(COMPONENT_ID_KEY);
    if (!r.has_value()) throw std::runtime_error("componentId not initialized");
    return r.value();
}

CGALMesh::Property_map<CGALMesh::Edge_index,   bool> FiberSurface::isImpassable()
{
    auto r = mesh.property_map<CGALMesh::Edge_index, bool>(IMPASSABLE_KEY);
    if (!r.has_value()) throw std::runtime_error("isImpassable not initialized");
    return r.value();
}


CGALMesh::Property_map<CGALMesh::Vertex_index, double> FiberSurface::edgeParam() const
{
    auto r = mesh.property_map<CGALMesh::Vertex_index, double>(EDGE_PARAM_KEY);
    if (!r.has_value()) throw std::runtime_error("edgeParam not initialized");
    return r.value();
}

CGALMesh::Property_map<CGALMesh::Face_index,   int> FiberSurface::tetId() const
{
    auto r = mesh.property_map<CGALMesh::Face_index, int>(TET_ID_KEY);
    if (!r.has_value()) throw std::runtime_error("tetId not initialized");
    return r.value();
}

CGALMesh::Property_map<CGALMesh::Face_index,   int> FiberSurface::sheetId() const
{
    auto r = mesh.property_map<CGALMesh::Face_index, int>(SHEET_ID_KEY);
    if (!r.has_value()) throw std::runtime_error("sheetId not initialized");
    return r.value();
}

CGALMesh::Property_map<CGALMesh::Face_index,   int> FiberSurface::componentId() const
{
    auto r = mesh.property_map<CGALMesh::Face_index, int>(COMPONENT_ID_KEY);
    if (!r.has_value()) throw std::runtime_error("componentId not initialized");
    return r.value();
}

CGALMesh::Property_map<CGALMesh::Edge_index,   bool> FiberSurface::isImpassable() const
{
    auto r = mesh.property_map<CGALMesh::Edge_index, bool>(IMPASSABLE_KEY);
    if (!r.has_value()) throw std::runtime_error("isImpassable not initialized");
    return r.value();
}
