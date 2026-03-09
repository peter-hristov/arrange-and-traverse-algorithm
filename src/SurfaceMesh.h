#pragma once

#include "./CGALTypedefs.h"

#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <CGAL/enum.h>
#include <iostream>
#include <stack>
#include <string>
#include <unordered_map>
#include <unordered_set>

#include "./FiberPoint.h"
#include "./TetMesh.h"
#include "./ReebSpace2.h"
#include "./Arrangement.h"
#include "./FiberGraph.h"

#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>

#include <CGAL/Polygon_mesh_processing/repair.h>
#include <CGAL/Polygon_mesh_processing/repair_degeneracies.h>

class SurfaceMesh
{
    public:
        CGALMesh mesh;
        CGALMesh::Property_map<CGALMesh::Vertex_index, double> edgeParam;

        CGALMesh::Property_map<CGALMesh::Face_index, int> tetId;
        CGALMesh::Property_map<CGALMesh::Face_index, int> sheetId;
        CGALMesh::Property_map<CGALMesh::Face_index, int> componentId;

        CGALMesh::Property_map<CGALMesh::Edge_index, bool> isImpassable;

        const double epsilon = 1e-10;

        SurfaceMesh()
        {

        }


        void print()
        {
            std::cout << "Vertices:\n";
            for (auto v : mesh.vertices())
            {
                const auto &p = mesh.point(v);
                double e = edgeParam[v];
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
                std::cout << "], tetId = " << tetId[f]
                    << ", sheetId = " << sheetId[f] << "\n";
            }

        }

        // Create mesh from a from a triangle soup
        SurfaceMesh(const std::vector<std::array<double, 3>> &vertexCoordinates, const std::vector<std::array<int, 3>> &triangles, const std::vector<double> &vertexEdgePara, const std::vector<int> &tetId)
        {
            // Unpack the points
            std::vector<CartesianPoint_3> points;
            for (auto& p : vertexCoordinates)
            {
                points.push_back(CartesianPoint_3(p[0],p[1],p[2]));
            }

            auto polygons = triangles;

            // Assume that the mesh is already cleaned up
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

            bool created;
            std::tie(this->edgeParam, created) = this->mesh.add_property_map<CGALMesh::Vertex_index,double>("v:edgeParam", -1.0);

            if (false == created)
            {
                throw std::runtime_error("EdgeParam property could not be added to the mesh.");
            }

            for (std::size_t i = 0; i < vertexEdgePara.size(); ++i)
            {
                this->edgeParam[CGALMesh::Vertex_index(i)] = vertexEdgePara[i];
            }


            std::tie(this->tetId, created) = this->mesh.add_property_map<CGALMesh::Face_index, int>("f:tetId", -1);
            
            if (false == created)
            {
                throw std::runtime_error("TetId property could not be added to the mesh.");
            }

            for (std::size_t i = 0; i < tetId.size(); ++i)
            {
                this->tetId[CGALMesh::Face_index(i)] = tetId[i];
            }

            std::tie(this->sheetId, created) = this->mesh.add_property_map<CGALMesh::Face_index, int>("f:sheetId", -1);

            if (false == created)
            {
                throw std::runtime_error("SheetId property could not be added to the mesh.");
            }

            std::tie(this->componentId, created) = this->mesh.add_property_map<CGALMesh::Face_index, int>("f:componentId", -1);

            if (false == created)
            {
                throw std::runtime_error("ComponentId property could not be added to the mesh.");
            }

            std::tie(this->isImpassable, created) = this->mesh.add_property_map<CGALMesh::Edge_index, bool>("f:isImpassable", false);

            if (false == created)
            {
                throw std::runtime_error("isImpassable property could not be added to the mesh.");
            }


        }









        // Generated by chat gpt
        std::array<double, 4> computeBarycentricCoordinates(const TetMesh &tetMesh, const int &tetId, const std::array<double, 3> &point3)
        {
            // Unpack tetrahedron vertices
            std::array<CartesianPoint_3,4> tetVertices;
            for(int i=0;i<4;i++)
            {
                const int vid = tetMesh.tetrahedra[tetId][i];
                tetVertices[i] = CartesianPoint_3(
                        tetMesh.vertexDomainCoordinates[vid][0],
                        tetMesh.vertexDomainCoordinates[vid][1],
                        tetMesh.vertexDomainCoordinates[vid][2]
                        );
            }

            const double x0=tetVertices[0].x(), y0=tetVertices[0].y(), z0=tetVertices[0].z();
            const double x1=tetVertices[1].x(), y1=tetVertices[1].y(), z1=tetVertices[1].z();
            const double x2=tetVertices[2].x(), y2=tetVertices[2].y(), z2=tetVertices[2].z();
            const double x3=tetVertices[3].x(), y3=tetVertices[3].y(), z3=tetVertices[3].z();
            const double xp=point3[0], yp=point3[1], zp=point3[2];

            // Determinant of the tetrahedron
            const double detT = 
                (x1-x0)*((y2-y0)*(z3-z0) - (y3-y0)*(z2-z0)) -
                (x2-x0)*((y1-y0)*(z3-z0) - (y3-y0)*(z1-z0)) +
                (x3-x0)*((y1-y0)*(z2-z0) - (y2-y0)*(z1-z0));

            std::array<double,4> w;

            w[0] = ((x1-xp)*((y2-yp)*(z3-zp) - (y3-yp)*(z2-zp)) -
                    (x2-xp)*((y1-yp)*(z3-zp) - (y3-yp)*(z1-zp)) +
                    (x3-xp)*((y1-yp)*(z2-zp) - (y2-yp)*(z1-zp))) / detT;

            w[1] = ((xp-x0)*((y2-y0)*(z3-z0) - (y3-y0)*(z2-z0)) -
                    (x2-x0)*((yp-y0)*(z3-z0) - (y3-y0)*(zp-z0)) +
                    (x3-x0)*((yp-y0)*(z2-z0) - (y2-y0)*(zp-z0))) / detT;

            w[2] = ((x1-x0)*((yp-y0)*(z3-z0) - (y3-y0)*(zp-z0)) -
                    (xp-x0)*((y1-y0)*(z3-z0) - (y3-y0)*(z1-z0)) +
                    (x3-x0)*((y1-y0)*(zp-z0) - (yp-y0)*(z1-z0))) / detT;

            w[3] = 1.0 - w[0] - w[1] - w[2];

            return w;
        }


        // Function to compute midpoint of a triangular face
        std::array<double, 3> triangleMidpoint(const CGALMesh::Face_index &f)
        {
            std::array<double, 3> midpoint = {0.0, 0.0, 0.0};
            int count = 0;

            // Iterate over the three vertices of the face
            for (auto v : vertices_around_face(mesh.halfedge(f), mesh))
            {
                const auto &p = mesh.point(v);
                midpoint[0] += p[0];
                midpoint[1] += p[1];
                midpoint[2] += p[2];
                ++count;
            }

            // Average coordinates
            midpoint[0] /= count;
            midpoint[1] /= count;
            midpoint[2] /= count;

            return midpoint;
        }









        int computeTriangleSheetId(TetMesh &tetMesh, Arrangement &singularArrangement, ReebSpace2 &reebSpace, const CGALMesh::Face_index &triangle)
        {
            const int tetId = this->tetId[triangle];

            const std::array<double, 3> midpoint = this->triangleMidpoint(triangle);
            const std::array<double, 4> barycentricCoordinates = this->computeBarycentricCoordinates(tetMesh, tetId, midpoint);

            const int a = tetMesh.tetrahedra[tetId][0];
            const int b = tetMesh.tetrahedra[tetId][1];
            const int c = tetMesh.tetrahedra[tetId][2];
            const int d = tetMesh.tetrahedra[tetId][3];

            const double ua = tetMesh.vertexCoordinatesF[a];
            const double va = tetMesh.vertexCoordinatesG[a];

            const double ub = tetMesh.vertexCoordinatesF[b];
            const double vb = tetMesh.vertexCoordinatesG[b];

            const double uc = tetMesh.vertexCoordinatesF[c];
            const double vc = tetMesh.vertexCoordinatesG[c];

            const double ud = tetMesh.vertexCoordinatesF[d];
            const double vd = tetMesh.vertexCoordinatesG[d];


            const double u = barycentricCoordinates[0] * ua + barycentricCoordinates[1] * ub + barycentricCoordinates[2] * uc + barycentricCoordinates[3] * ud;
            const double v = barycentricCoordinates[0] * va + barycentricCoordinates[1] * vb + barycentricCoordinates[2] * vc + barycentricCoordinates[3] * vd;



            //const std::set<int> tetTriangleIds = {
                //tetMesh.triangleIndices.at({a, b, c}),
                //tetMesh.triangleIndices.at({a, b, d}),
                //tetMesh.triangleIndices.at({a, c, d}),
                //tetMesh.triangleIndices.at({b, c, d}),
            //};

            //return reebSpace.computeFiberGraphReverse(tetMesh, singularArrangement, {u, v}, tetTriangleIds);


            // Compute the fiber graph
            //
            const auto fg = reebSpace.computeFiberGraph(tetMesh, singularArrangement, {u, v});


            // Find which componnt we are in
            //
            const std::vector<int> tetTriangleIds = {
                tetMesh.triangleIndices.at({a, b, c}),
                tetMesh.triangleIndices.at({a, b, d}),
                tetMesh.triangleIndices.at({a, c, d}),
                tetMesh.triangleIndices.at({b, c, d}),
            };

            for (const int triangleId : tetTriangleIds)
            {
                if (fg.componentRoot.contains(triangleId))
                {
                    return reebSpace.correspondenceGraphDS.find(fg.componentRoot.at(triangleId));
                }

                //printf("The barycentric coordinate of triangle id %d with midpoint (%f, %f, %f) are (%f, %f, %f, %f).\nThe range value is (%f, %f) and the sheet is %d\n", i, midpoint[0], midpoint[1], midpoint[2], barycentricCoordinates[0], barycentricCoordinates[1], barycentricCoordinates[2], barycentricCoordinates[3], u, v, triangleSheet[i]);
            }

            std::cerr << "Sheet for fiber surface triangle could not be found!\n";
            return -1;
        }




        // Back for debuggin connected component compututaion, otherwise CGAL's CGAL::Polygon_mesh_processing::connected_components
        std::size_t computeConnectedComponentsBFS(CGALMesh::Property_map<CGALMesh::Face_index, int>& componentMap)
        {
            const CGALMesh& mesh = this->mesh;

            // Initialize all faces as unvisited (-1)
            for (auto f : mesh.faces())
                componentMap[f] = -1;

            int componentId = 0;

            for (auto startFace : mesh.faces())
            {
                if (componentMap[startFace] != -1) continue; // already visited

                printf("----------------------------------------------- At Component %d\n\n", componentId);

                // BFS from startFace
                std::queue<CGALMesh::Face_index> queue;
                queue.push(startFace);
                componentMap[startFace] = componentId;

                while (!queue.empty())
                {
                    auto f = queue.front();
                    queue.pop();

                    printf("------------------------- At Face %d\n", f.idx());

                    // Iterate over the 3 halfedges of this face
                    auto h = mesh.halfedge(f);
                    for (int i = 0; i < 3; ++i, h = mesh.next(h))
                    {
                        // Skip if this edge is impassable
                        auto e = mesh.edge(h);
                        if (this->isImpassable[e]) continue;

                        // Get the opposite face
                        auto opp = mesh.opposite(h);
                        if (mesh.is_border(opp)) continue; // no face on the other side

                        auto neighbor = mesh.face(opp);
                        if (componentMap[neighbor] != -1) continue; // already visited

                        printf("--- Adding neighbour %d\n", neighbor.idx());

                        componentMap[neighbor] = componentId;
                        queue.push(neighbor);
                    }
                }

                componentId++;
            }

            return static_cast<std::size_t>(componentId);
        }




        void computeTriangleSheets(TetMesh &tetMesh, Arrangement &singularArrangement, ReebSpace2 reebSpace)
        {

            std::size_t numComponents2 = CGAL::Polygon_mesh_processing::connected_components(
                    this->mesh,
                    this->componentId,
                    CGAL::parameters::edge_is_constrained_map(this->isImpassable)
                    );

            //std::size_t numComponents2 = computeConnectedComponentsBFS(this->componentId);


            // 2. Get a representative triangle per component
            //
            std::set<int> usedComponents;
            std::vector<std::pair<CGALMesh::Face_index, int>> componentRepresentatives;
            for (auto face : mesh.faces())
            {
                const int componentId = this->componentId[face];

                if (false == usedComponents.contains(componentId))
                {
                    usedComponents.insert(componentId);
                    componentRepresentatives.push_back({face, componentId});
                }
            }

            std::cout << "Computing " << componentRepresentatives.size() << " fiber graphs...\n";

            // 3. Compute one flexible fiber per representative triangle
            //
            std::vector<int> componentSheets(componentRepresentatives.size());
            //#pragma omp parallel for schedule(dynamic)
            for (auto &[face, componentId] : componentRepresentatives)
            {
                componentSheets[componentId] = this->computeTriangleSheetId(tetMesh, singularArrangement, reebSpace, face);
            }

            // 4. Set up the sheetIds of each triangle based on the connected component
            //
            for (auto face : mesh.faces())
            {
                const int componentId = this->componentId[face];
                this->sheetId[face] =  componentSheets[componentId];



                //this->sheetId[face] = this->computeTriangleSheetId(tetMesh, singularArrangement, reebSpace, face);
                
                //this->sheetId[face] =  componentId;
                //const int realSheetId = this->computeTriangleSheetId(tetMesh, singularArrangement, reebSpace, face);

                //if (this->sheetId[face] != realSheetId)
                //{
                    //throw std::runtime_error("Triangle sheet Id not the same as its component id.");
                //}
            }
        }



        int computeTriangleSheetId2(TetMesh &tetMesh, Arrangement &singularArrangement, ReebSpace2 &reebSpace, const CGALMesh::Face_index &triangle, const std::vector<std::tuple<K::FT, int, int>> &intersectedSegments, const Segment_2 &controlSegment)
        {
            // 1. Compute the edgePara at the center of the triangle
            double midPointAlpha = 0.0;
            for (const auto v : mesh.vertices_around_face(mesh.halfedge(triangle))) 
            {
                midPointAlpha += this->edgeParam[v];
            }
            midPointAlpha /= 3.0;

            // 2. Compute the fiber graph at the alpha in the range
            const FiberGraph fg = reebSpace.computeFiberGraph3(tetMesh, singularArrangement, controlSegment, midPointAlpha, intersectedSegments);

            // 3. Determine which fiber component contains a triangle from the tet
            const int tetId = this->tetId[triangle];

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

            for (const int triangleId : tetTriangleIds)
            {
                if (fg.componentRoot.contains(triangleId))
                {
                    return reebSpace.correspondenceGraphDS.find(fg.componentRoot.at(triangleId));
                }

                //printf("The barycentric coordinate of triangle id %d with midpoint (%f, %f, %f) are (%f, %f, %f, %f).\nThe range value is (%f, %f) and the sheet is %d\n", i, midpoint[0], midpoint[1], midpoint[2], barycentricCoordinates[0], barycentricCoordinates[1], barycentricCoordinates[2], barycentricCoordinates[3], u, v, triangleSheet[i]);
            }

            return -1;
        }




        void computeTriangleSheets2(TetMesh &tetMesh, Arrangement &singularArrangement, ReebSpace2 reebSpace, const std::vector<std::tuple<K::FT, int, int>> &intersectedSegments, const Segment_2 &controlSegment)
        {

            std::size_t numComponents2 = CGAL::Polygon_mesh_processing::connected_components(
                    this->mesh,
                    this->componentId,
                    CGAL::parameters::edge_is_constrained_map(this->isImpassable)
                    );

            //std::size_t numComponents2 = computeConnectedComponentsBFS(this->componentId);


            // 2. Get a representative triangle per component
            //
            std::set<int> usedComponents;
            std::vector<std::pair<CGALMesh::Face_index, int>> componentRepresentatives;
            for (auto face : mesh.faces())
            {
                const int componentId = this->componentId[face];

                if (false == usedComponents.contains(componentId))
                {
                    usedComponents.insert(componentId);
                    componentRepresentatives.push_back({face, componentId});
                }
            }

            // 3. Compute one flexible fiber per representative triangle
            //
            std::vector<int> componentSheets(componentRepresentatives.size());

            #pragma omp parallel for schedule(dynamic)
            for (auto &[face, componentId] : componentRepresentatives)
            {
                componentSheets[componentId] = this->computeTriangleSheetId2(tetMesh, singularArrangement, reebSpace, face, intersectedSegments, controlSegment);
            }

            // 4. Set up the sheetIds of each triangle based on the connected component
            //
            for (auto face : mesh.faces())
            {
                const int componentId = this->componentId[face];
                this->sheetId[face] =  componentSheets[componentId];


                //this->sheetId[face] = this->computeTriangleSheetId(tetMesh, singularArrangement, reebSpace, face);
                
                //this->sheetId[face] =  componentId;
                //const int realSheetId = this->computeTriangleSheetId(tetMesh, singularArrangement, reebSpace, face);

                //if (this->sheetId[face] != realSheetId)
                //{
                    //throw std::runtime_error("Triangle sheet Id not the same as its component id.");
                //}
            }
        }










        // Helper: interpolate two points and scalar values
        CartesianPoint_3 interpolate_vertex(CGALMesh::Vertex_index v0, CGALMesh::Vertex_index v1, double isovalue)
        {
            const double e0 = edgeParam[v0];
            const double e1 = edgeParam[v1];

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


        CGALMesh::Property_map<CGALMesh::Vertex_index, int> getVertexColours(CGALMesh &cgalMesh, const double isovalue)
        {
            auto [vertexColour, created] = cgalMesh.add_property_map<CGALMesh::Vertex_index, int>("v:colour", -2);
            if (!created)
                throw std::runtime_error("Could not make mesh colour array.");

            for (auto v : cgalMesh.vertices())
            {
                const double e = edgeParam[v];
                const double diff = e - isovalue;

                if (std::abs(diff) <= this->epsilon)
                    vertexColour[v] = 0;       // on the isovalue
                else if (diff < 0.0)
                    vertexColour[v] = -1;      // below
                else
                    vertexColour[v] = +1;      // above
            }

            return vertexColour;
        }


        void subdivideMesh(const std::vector<double> &isovalues)
        {
            for (int i = 0 ; i < isovalues.size() ; i++)
            {
                this->subdivideMeshOnce(isovalues[i]);
            }

            this->triangulateMesh();
        }

        void subdivideMeshOnce(const double isovalue)
        {
            // 1. Compute the colours of all the vertices
            //
            CGALMesh::Property_map<CGALMesh::Vertex_index, int> vertexColour = getVertexColours(this->mesh, isovalue);

            // 2. Collect crossing edges first (don't modify while iterating!)
            //
            std::unordered_set<CGALMesh::Edge_index> activeEdges;

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

            std::unordered_set<CGALMesh::Face_index> facesToSplit;

            // 3. Split the edges
            //
            for (const auto &e : activeEdges)
            {
                auto h = this->mesh.halfedge(e);

                const auto v0 = this->mesh.source(h);
                const auto v1 = this->mesh.target(h);
                const double val0 = edgeParam[v0];
                const double val1 = edgeParam[v1];

                const CartesianPoint_3 edgeVertex = interpolate_vertex(v0, v1, isovalue);
                const auto hNew = CGAL::Euler::split_edge(h, mesh);
                const auto vNew = mesh.target(hNew);

                // Copy over the vertex values for the new points
                mesh.point(vNew) = edgeVertex;
                edgeParam[vNew] = isovalue;
                vertexColour[vNew] = 0;

                // Make whether this is impassable
                const auto isEdgeImpassable = this->isImpassable[e];
                this->isImpassable[mesh.edge(hNew)] = isEdgeImpassable;
                this->isImpassable[mesh.edge(mesh.next(hNew))] = isEdgeImpassable;

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
                const int faceTetId = this->tetId[f];

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
                    throw std::runtime_error("There should be 2 active half-edges.");
                }

                // Split
                const auto h0 = activeHalfEdges[0];
                const auto h1 = activeHalfEdges[1];
                const auto newH = CGAL::Euler::split_face(h0, h1, this->mesh);

                // Make the new edge impassable
                this->isImpassable[mesh.edge(newH)] = true;

                // Write tet array
                const auto newF0 = mesh.face(newH);
                const auto newF1 = mesh.face(mesh.opposite(newH));

                this->tetId[newF0] = faceTetId;
                this->tetId[newF1] = faceTetId;

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

        void triangulateMesh()
        {
            struct Visitor : CGAL::Polygon_mesh_processing::Triangulate_faces::Default_visitor<CGALMesh>
            {
                int tetId_val;
                CGALMesh::Property_map<CGALMesh::Face_index, int>& tetId;
                Visitor(int tetId, CGALMesh::Property_map<CGALMesh::Face_index, int>& t) : tetId_val(tetId), tetId(t) {}
                void after_subface_created(CGALMesh::Face_index f) { tetId[f] = tetId_val; }
            }; 

            for (const auto f : mesh.faces())
            {
                const int faceTetId = this->tetId[f];

                if (mesh.degree(f) == 3)
                {
                    continue;
                }

                Visitor visitor(faceTetId, this->tetId);
                // Vistor that sets the tetId of the newly created faces

                CGAL::Polygon_mesh_processing::triangulate_face(f, mesh, CGAL::parameters::visitor(visitor));
            }
        }

        void repairMesh()
        {
            // Finish up with some postprocessing
            mesh.collect_garbage(); // before calling connected_components
            CGAL::Polygon_mesh_processing::orient(this->mesh);

            // 1. Remove degenerate faces (zero or near-zero area)
            CGAL::Polygon_mesh_processing::remove_degenerate_faces(mesh);

            // 2. Remove degenerate edges (edges shorter than a threshold)
            CGAL::Polygon_mesh_processing::remove_degenerate_edges(mesh);
        }

        void validateMesh()
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

























        // Marching triangles
        //SurfaceMesh splitSingularTriangles(const double isovalue)
        //{
            //// Compute the colours of all the vertices
            ////
            //CGALMesh::Property_map<CGALMesh::Vertex_index, int> vertexColour = getVertexColours(this->mesh, isovalue);

            //// Set up the new mesh
            ////
            //SurfaceMesh newMesh;

            //bool created;
            //std::tie(newMesh.edgeParam, created) = newMesh.mesh.add_property_map<CGALMesh::Vertex_index,double>("v:edgeParam", -1.0);
            //if (false == created)
            //{
                //throw std::runtime_error("Could not make new mesh edgeParam array.");
            //}

            //std::map<CGALMesh::Vertex_index, CGALMesh::Vertex_index> oldToNewVertexMap;

            //// Copy over the previous vertices and their properties and build oldToNewVertexMap
            ////
            //for (auto v : this->mesh.vertices())
            //{
                //const auto &point = this->mesh.point(v);
                //const double e = this->edgeParam[v];

                //CGALMesh::Vertex_index newVertexId = newMesh.mesh.add_vertex(CartesianPoint_3(point[0], point[1], point[2]));
                //newMesh.edgeParam[newVertexId] = e;

                //oldToNewVertexMap[v] = newVertexId;
            //}


            //// Old Mesh ID and NEW mesh vertex
            //std::map<CGALMesh::Edge_index, CGALMesh::Vertex_index> edgeVertexMap;

            //// Find all edges crossing the isovalue
            ////
            //for (auto e : this->mesh.edges())
            //{
                //const auto h = this->mesh.halfedge(e);

                //const auto v0 = this->mesh.source(h);
                //const auto v1 = this->mesh.target(h);

                //const double val0 = edgeParam[v0];
                //const double val1 = edgeParam[v1];

                //if (vertexColour[v0] * vertexColour[v1] == -1) // edge crosses isovalue
                //{
                    //const CartesianPoint_3 edgeVertex = interpolate_vertex(v0, v1, isovalue);

                    //CGALMesh::Vertex_index edgeVertexIndex = newMesh.mesh.add_vertex(edgeVertex);
                    //edgeVertexMap[e] = edgeVertexIndex;

                    //newMesh.edgeParam[edgeVertexIndex] = isovalue;
                //}
            //}


            //// Add new faces to the mesh
            ////
            //std::tie(newMesh.tetId, created) = newMesh.mesh.add_property_map<CGALMesh::Face_index, int>("f:tetId", -1);
            //std::tie(newMesh.sheetId, created) = newMesh.mesh.add_property_map<CGALMesh::Face_index, int>("f:sheetId", -1);
            //std::vector<std::pair<CGALMesh::Vertex_index, CGALMesh::Vertex_index>> constraintEdges = buildNewFaces(*this, newMesh, oldToNewVertexMap, edgeVertexMap, vertexColour);

            //// Mark immpassable edges
            ////
            //std::tie(newMesh.isImpassable, created) = newMesh.mesh.add_property_map<CGALMesh::Edge_index, bool>("f:isImpassable", false);
            //buildImpassableEdges( *this, newMesh, oldToNewVertexMap, edgeVertexMap, constraintEdges, isImpassable);

            //// Finish up with some postprocessing
            ////
            //mesh.collect_garbage(); // before calling connected_components
            //CGAL::Polygon_mesh_processing::orient(newMesh.mesh);

            //// Make sure the edge is valid
            ////
            //if (false == CGAL::is_valid_polygon_mesh(newMesh.mesh))
            //{
                //throw std::runtime_error("New mesh is not valid.");
            //}

            ////std::cout << "Previous number of faces : " << mesh.num_faces() << std::endl;
            ////std::cout << "New number of faces : " << newMesh.mesh.num_faces() << std::endl;

            //return newMesh;
        //}



            //std::vector<std::pair<CGALMesh::Vertex_index, CGALMesh::Vertex_index>> buildNewFaces(
                    //const SurfaceMesh &oldMesh, 
                    //SurfaceMesh &newMesh, 
                    //const std::map<CGALMesh::Vertex_index, CGALMesh::Vertex_index> &oldToNewVertexMap, 
                    //const std::map<CGALMesh::Edge_index, CGALMesh::Vertex_index> &edgeVertexMap, 
                    //const CGALMesh::Property_map<CGALMesh::Vertex_index, int> &vertexColour
                    //)
            //{

                //// Return value
                //std::vector<std::pair<CGALMesh::Vertex_index, CGALMesh::Vertex_index>> constraintEdges;

                //std::vector<int> newTetIds;
                //std::vector<std::array<CGALMesh::Vertex_index, 3>> newFaces;

                //// 3. Collect faces that are affected by splits
                ////
                //for (auto f : this->mesh.faces())
                //{
                    //const int tetId = this->tetId[f];

                    //std::vector<CGALMesh::Vertex_index> faceVertices;
                    //std::vector<CGALMesh::Halfedge_index> activeHalfEdges;

                    //for (auto h : halfedges_around_face(mesh.halfedge(f), mesh))
                    //{
                        //auto e = mesh.edge(h);

                        //faceVertices.push_back(mesh.source(h));

                        //if (edgeVertexMap.contains(e))
                        //{
                            //activeHalfEdges.push_back(h);
                        //}
                    //}

                    //if (activeHalfEdges.size() == 0)
                    //{
                        //const CGALMesh::Vertex_index a = oldToNewVertexMap.at(faceVertices[0]);
                        //const CGALMesh::Vertex_index b = oldToNewVertexMap.at(faceVertices[1]);
                        //const CGALMesh::Vertex_index c = oldToNewVertexMap.at(faceVertices[2]);

                        //newFaces.push_back({a, b, c});
                        //newTetIds.push_back(tetId);
                    //}


                    ////
                    ////           c
                    ////          /|\     /\
                    ////         / | \     \  next(h)
                    ////        /  |  \     \
                    ////       /   |   \     \
                    ////      /    |    \
                    ////   a /_____|_____\ b
                    ////           d
                    ////       
                    ////       -------->
                    ////           h
                    ////

                    //else if (activeHalfEdges.size() == 1)
                    //{
                        //const auto h = activeHalfEdges[0];

                        //const CGALMesh::Vertex_index a = oldToNewVertexMap.at(mesh.source(h));
                        //const CGALMesh::Vertex_index b = oldToNewVertexMap.at(mesh.target(h));
                        //const CGALMesh::Vertex_index c = oldToNewVertexMap.at(mesh.target(mesh.next(h)));

                        //const CGALMesh::Vertex_index d = edgeVertexMap.at(mesh.edge(h));

                        //newFaces.push_back({a, d, c});
                        //newTetIds.push_back(tetId);

                        //newFaces.push_back({b, c, d});
                        //newTetIds.push_back(tetId);

                        //constraintEdges.push_back({c, d});

                        //if (vertexColour[c] != 0)
                        //{
                            //throw std::runtime_error("Vertec c should be gray.");
                        //}

                        //if (vertexColour[a] * vertexColour[b] != -1)
                        //{
                            //throw std::runtime_error("Vertices a and b should have different colours, neither gray.");
                        //}
                    //}



                    ////
                    ////                c
                    ////               /|
                    ////              / |
                    ////             /  |      ^
                    ////      /     /   |      |
                    ////  h1 /  d1 /----| d0   | h0
                    ////    /     /\    |      | 
                    ////   \/    /  \   |      
                    ////        /    \  |
                    ////       /      \ |
                    ////      /________\|
                    ////     a           b
                    ////
                    ////

                    //else if (activeHalfEdges.size() == 2)
                    //{
                        //// Find the 3rd vertex

                        //auto h0 = activeHalfEdges[0];
                        //auto h1 = activeHalfEdges[1];

                        //// Swap to make sure that next(h1) = h2
                        //if (mesh.next(h1) == h0)
                        //{
                            //std::swap(h0, h1);
                        //}

                        //const CGALMesh::Vertex_index a = oldToNewVertexMap.at(mesh.target(h1));
                        //const CGALMesh::Vertex_index b = oldToNewVertexMap.at(mesh.source(h0));
                        //const CGALMesh::Vertex_index c = oldToNewVertexMap.at(mesh.target(h0));

                        //const CGALMesh::Vertex_index d0 = edgeVertexMap.at(mesh.edge(h0));
                        //const CGALMesh::Vertex_index d1 = edgeVertexMap.at(mesh.edge(h1));

                        //newFaces.push_back({d1, a, b});
                        //newTetIds.push_back(tetId);

                        //newFaces.push_back({b, d0, d1});
                        //newTetIds.push_back(tetId);

                        //newFaces.push_back({d0, c, d1});
                        //newTetIds.push_back(tetId);

                        //constraintEdges.push_back({d0, d1});

                        //if (vertexColour[a] * vertexColour[b] * vertexColour[c] == 0)
                        //{
                            //throw std::runtime_error("Neither of a, b, c should be gray.");
                        //}

                        //if (vertexColour[a] * vertexColour[c] != -1)
                        //{
                            //throw std::runtime_error("Vertices a and c should have different colours, neither gray.");
                        //}

                        //if (vertexColour[a] != vertexColour[b])
                        //{
                            //throw std::runtime_error("Vertices a and b should have the same colour.");
                        //}
                    //}
                    //else
                    //{
                        //throw std::runtime_error("Degenerate remeshing case.");
                    //}
                //}

                //for (int i = 0 ; i < newFaces.size() ; i++)
                //{
                    //auto newFace = newMesh.mesh.add_face(newFaces[i]);

                    //if (newFace == CGALMesh::null_face())
                    //{
                        //throw std::runtime_error("Failed to add a new triangle to the msh.");
                    //}

                    //newMesh.tetId[newFace] = newTetIds[i];
                //}

                //return constraintEdges;
            //}



            //void buildImpassableEdges(
                    //const SurfaceMesh &oldMesh, 
                    //SurfaceMesh &newMesh, 
                    //const std::map<CGALMesh::Vertex_index, CGALMesh::Vertex_index> &oldToNewVertexMap, 
                    //const std::map<CGALMesh::Edge_index, CGALMesh::Vertex_index> &edgeVertexMap, 
                    //const std::vector<std::pair<CGALMesh::Vertex_index, CGALMesh::Vertex_index>> &impassableEdges,
                    //CGALMesh::Property_map<CGALMesh::Edge_index, bool> &isImpassable)
            //{
                //for (auto e : this->mesh.edges())
                //{
                    //// Is the edge impassable
                    //const bool isEdgeImpassable = this->isImpassable[e];

                    //const CGALMesh::Halfedge_index h = mesh.halfedge(e);
                    //const CGALMesh::Vertex_index a = mesh.source(h);
                    //const CGALMesh::Vertex_index b = mesh.target(h);

                    //const CGALMesh::Vertex_index aNew = oldToNewVertexMap.at(a);
                    //const CGALMesh::Vertex_index bNew = oldToNewVertexMap.at(b);

                    //// If the edge has been subdivided
                    //if (edgeVertexMap.contains(e))
                    //{
                        //// The point of subdivision
                        //const CGALMesh::Vertex_index cNew = edgeVertexMap.at(e);

                        //const auto hNew0 = newMesh.mesh.halfedge(aNew, cNew);

                        //if (hNew0 == CGALMesh::null_halfedge())
                        //{
                            //throw std::runtime_error("Half-edge not found.");
                        //}

                        //const auto hNew1 = newMesh.mesh.halfedge(cNew, bNew);

                        //if (hNew1 == CGALMesh::null_halfedge())
                        //{
                            //throw std::runtime_error("Half-edge not found.");
                        //}

                        //const auto eNew0 = newMesh.mesh.edge(hNew0);
                        //const auto eNew1 = newMesh.mesh.edge(hNew1);

                        //newMesh.isImpassable[eNew0] = isEdgeImpassable;
                        //newMesh.isImpassable[eNew1] = isEdgeImpassable;
                    //}

                    //// If the edge has not been subdivided, just copy over the old value
                    //else
                    //{
                        //const auto hNew = newMesh.mesh.halfedge(aNew, bNew);

                        //if (hNew == CGALMesh::null_halfedge())
                        //{
                            //throw std::runtime_error("Half-edge not found.");
                        //}

                        //const auto eNew = newMesh.mesh.edge(hNew);

                        //newMesh.isImpassable[eNew] = isEdgeImpassable;
                    //}
                //}

                //// Fill in the new impassable edges
                ////
                //for (const auto &[aNew, bNew] : impassableEdges)
                //{
                    //const auto hNew = newMesh.mesh.halfedge(aNew, bNew);
                    //const auto eNew = newMesh.mesh.edge(hNew);
                    //newMesh.isImpassable[eNew] = true;
                //}
            //}





            ////CGALMesh::Property_map<CGALMesh::Vertex_index, int> getVertexColours(CGALMesh &cgalMesh, const double isovalue)
            ////{
                ////auto [vertexColour, created] = cgalMesh.add_property_map<CGALMesh::Vertex_index, int>("v:colour", -2);

                ////if (false == created)
                ////{
                    ////throw std::runtime_error("Could not make mesh colour array.");
                ////}

                ////for (auto v : cgalMesh.vertices())
                ////{
                    ////const double e = edgeParam[v];

                    ////if (CGAL::compare(e, isovalue) == CGAL::SMALLER)
                    ////{
                        ////vertexColour[v] = -1;
                    ////}
                    ////else if (CGAL::compare(e, isovalue) == CGAL::EQUAL)
                    ////{
                        //////std::cerr << "GRAY VERTEX!";
                        ////vertexColour[v] = 0;
                    ////}
                    ////else
                    ////{
                        ////vertexColour[v] = +1;
                    ////}
                ////}

                ////return vertexColour;
            ////}
            
};
