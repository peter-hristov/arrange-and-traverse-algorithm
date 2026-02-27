#pragma once

#include "./CGALTypedefs.h"

#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <CGAL/enum.h>
#include <iostream>
#include <stack>
#include <unordered_map>
#include <unordered_set>

#include "./FiberPoint.h"
#include "./TetMesh.h"
#include "./ReebSpace2.h"
#include "./Arrangement.h"
#include "./Fiber.h"

#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>

class SurfaceMesh
{
    public:
        CGALMesh mesh;
        CGALMesh::Property_map<CGALMesh::Vertex_index, double> edgeParam;
        CGALMesh::Property_map<CGALMesh::Face_index, int> tetId;
        CGALMesh::Property_map<CGALMesh::Face_index, int> sheetId;

        SurfaceMesh()
        {

        }

        // Create mesh from a from a triangle soup
        SurfaceMesh(const std::vector<std::array<double, 3>> &vertexCoordinates, std::vector<std::array<int, 3>> triangles, const std::vector<double> &vertexEdgePara, const std::vector<int> &tetId)
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


        }





        //std::vector<FiberPoint> getFiberPoints(const std::vector<int> triangleSheets = {})
        //{
            //std::vector<FiberPoint> allFiberPoints;

            //std::array<float, 3> triangleColour{1.0, 1.0, 0.0};
            //for (auto f : mesh.faces())
            //{
                //for (auto v : vertices_around_face(mesh.halfedge(f), mesh))
                //{
                    //std::array<double, 3> point = { mesh.point(v)[0], mesh.point(v)[1], mesh.point(v)[2] };

                    //allFiberPoints.push_back(FiberPoint(
                                //point,
                                //triangleColour, 
                                //1,
                                //-1
                                //));
                //}
            //}
            //return allFiberPoints;
        //}


        //std::vector<std::array<double, 3>> vertexCoordinates;
        //std::vector<std::array<int, 3>> triangles;
        ///
        //std::vector<double> edgeParam;
        //std::vector<bool> isVertexSingular;
        //std::vector<int> triangleTetId;
        //std::vector<int> triangleComponentId;




        std::vector<FiberPoint> getFiberPoints()
        {
            std::vector<FiberPoint> allFiberPoints;

            for (auto f : mesh.faces())
            {

                const int sheetId = this->sheetId[f];

                // Default triangle colour
                std::array<float, 3> triangleColour{1.0, 1.0, 0.0};

                if (sheetId != -1)
                {
                    triangleColour = fiber::fiberColours[sheetId % fiber::fiberColours.size()];
                }

                for (auto v : vertices_around_face(mesh.halfedge(f), mesh))
                {
                    std::array<double, 3> point = { mesh.point(v)[0], mesh.point(v)[1], mesh.point(v)[2] };

                    allFiberPoints.push_back(FiberPoint(
                                point,
                                triangleColour, 
                                1,
                                -1
                                ));
                }

            }

            return allFiberPoints;

        }

        //void print()
        //{
            //std::cout << "Number of vertices = " << this->vertexCoordinates.size() << "\n";
            //std::cout << "Number of triangles = " << this->triangles.size() << "\n";

            //for (int i = 0 ; i < this->vertexCoordinates.size() ; i++)
            //{
                //printf("Vertex %d with coordinates (%f, %f, %f) and edgePara %f and isSingular %d.\n", i, vertexCoordinates[i][0], vertexCoordinates[i][1], vertexCoordinates[i][2], edgeParam[i], (int)isVertexSingular[i]);
            //}

            //for (int i = 0 ; i < this->triangles.size() ; i++)
            //{
                //printf("Triangle %d with vertices (%d, %d, %d) and tetId %d\n", i, triangles[i][0], triangles[i][1], triangles[i][2], triangleTetId[i]);
            //}

        //}




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

        //std::array<double, 3> triangleMidpoint(
                //const std::array<double, 3>& v0,
                //const std::array<double, 3>& v1,
                //const std::array<double, 3>& v2)
        //{
            //return {
                //(v0[0] + v1[0] + v2[0]) / 3.0,
                //(v0[1] + v1[1] + v2[1]) / 3.0,
                //(v0[2] + v1[2] + v2[2]) / 3.0
            //};
        //}


        // Function to compute midpoint of a triangular face
        std::array<double, 3> triangleMidpoint(CGALMesh::Face_index f)
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

        int computeTriangleSheetId(TetMesh &tetMesh, Arrangement &singularArrangement, ReebSpace2 reebSpace, const CGALMesh::Face_index triangle)
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

            return -1;
        }



        void computeTriangleSheets(TetMesh &tetMesh, Arrangement &singularArrangement, ReebSpace2 reebSpace)
        {
            // Compute the connectivity of the mesh
            //
            std::vector<std::size_t> component(num_faces(mesh));
            std::size_t num = CGAL::Polygon_mesh_processing::connected_components(mesh, CGAL::make_property_map(component));

            std::cout << "There are " << num << " components\n";

            std::map<int, int> componentToSheetId;

            for (auto face : mesh.faces())
            {
                const int componentId = component[face];

                //if (false == componentToSheetId.contains(componentId))
                //{
                    //componentToSheetId[componentId] = this->computeTriangleSheetId(tetMesh, singularArrangement, reebSpace, face);
                //}

                //this->sheetId[face] =  componentToSheetId.at(componentId);

                this->sheetId[face]    = this->computeTriangleSheetId(tetMesh, singularArrangement, reebSpace, face);

                //std::cout << "The sheet id id " << this->sheetId[face] << std::endl;
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


        // Marching triangles
        SurfaceMesh splitSingularTriangles(const double isovalue)
        {
            // Colour triangles as below, at and above
            auto [vertexColour, created] = this->mesh.add_property_map<CGALMesh::Vertex_index, int>("v:colour", -2);
            
            if (false == created)
            {
                throw std::runtime_error("Could not make mesh array.");
            }

            for (auto v : this->mesh.vertices())
            {
                const double e = edgeParam[v];

                if (CGAL::compare(e, isovalue) == CGAL::SMALLER)
                {
                    vertexColour[v] = -1;
                }
                else if (CGAL::compare(e, isovalue) == CGAL::EQUAL)
                {
                    std::cerr << "GRAY VERTEX!";
                    vertexColour[v] = 0;
                }
                else
                {
                    vertexColour[v] = +1;
                }
            }


            SurfaceMesh newMesh;
            std::tie(newMesh.edgeParam, created) = newMesh.mesh.add_property_map<CGALMesh::Vertex_index,double>("v:edgeParam", -1.0);

            std::map<CGALMesh::Vertex_index, CGALMesh::Vertex_index> oldToNewVertexMap;

            // Copy over the vertices
            //
            for (auto v : this->mesh.vertices())
            {
                const auto &point = this->mesh.point(v);
                const double e = this->edgeParam[v];

                CGALMesh::Vertex_index newVertexId = newMesh.mesh.add_vertex(CartesianPoint_3(point[0], point[1], point[2]));
                newMesh.edgeParam[newVertexId] = e;

                oldToNewVertexMap[v] = newVertexId;
            }


            // Old Mesh ID and NEW mesh vertex
            std::map<CGALMesh::Edge_index, CGALMesh::Vertex_index> edgeVertexMap;

            // 1. Find all edges crossing the isovalue
            for (auto e : this->mesh.edges())
            {
                const auto h = this->mesh.halfedge(e);

                const auto v0 = this->mesh.source(h);
                const auto v1 = this->mesh.target(h);

                const double val0 = edgeParam[v0];
                const double val1 = edgeParam[v1];

                if (vertexColour[v0] * vertexColour[v1] == -1) // edge crosses isovalue
                {
                    const CartesianPoint_3 edgeVertex = interpolate_vertex(v0, v1, isovalue);

                    CGALMesh::Vertex_index edgeVertexIndex = newMesh.mesh.add_vertex(edgeVertex);
                    newMesh.edgeParam[edgeVertexIndex] = isovalue;
                    edgeVertexMap[e] = edgeVertexIndex;

                    //const CGALMesh::Vertex_index edgeVertexIndex = mesh.add_vertex(edgeVertex);
                    //edgeVertexMap[e] = edgeVertexIndex;
                    //edgeParam[edgeVertexIndex] = isovalue;
                }
            }


            std::tie(newMesh.tetId, created) = newMesh.mesh.add_property_map<CGALMesh::Face_index, int>("f:tetId", -1);
            std::tie(newMesh.sheetId, created) = newMesh.mesh.add_property_map<CGALMesh::Face_index, int>("f:sheetId", -1);


            std::vector<std::array<CGALMesh::Vertex_index, 3>> newFaces;

            std::vector<int> newTetIds;

            // 3. Collect faces that are affected by splits
            //
            for (auto f : this->mesh.faces())
            {
                const int tetId = this->tetId[f];

                std::vector<CGALMesh::Halfedge_index> activeHalfEdges;
                for (auto h : halfedges_around_face(mesh.halfedge(f), mesh))
                {
                    auto e = mesh.edge(h);

                    if (edgeVertexMap.contains(e))
                    {
                        activeHalfEdges.push_back(h);
                    }
                }

                if (activeHalfEdges.size() == 1)
                {
                    const auto h = activeHalfEdges[0];

                    const CGALMesh::Vertex_index a = oldToNewVertexMap[mesh.source(h)];
                    const CGALMesh::Vertex_index b = oldToNewVertexMap[mesh.target(h)];
                    const CGALMesh::Vertex_index c = oldToNewVertexMap[mesh.target(mesh.next(h))];
                    const CGALMesh::Vertex_index d = edgeVertexMap.at(mesh.edge(h));

                    newFaces.push_back({b, c, d});
                    newTetIds.push_back(tetId);

                    newFaces.push_back({a, d, c});
                    newTetIds.push_back(tetId);
                }

                if (activeHalfEdges.size() == 2)
                {
                    // Find the 3rd vertex

                    auto h1 = activeHalfEdges[0];
                    auto h2 = activeHalfEdges[1];

                    // Swap to make sure that next(h1) = h2
                    if (mesh.next(h2) == h1)
                    {
                        std::swap(h1, h2);
                    }

                    const CGALMesh::Vertex_index a = oldToNewVertexMap[mesh.target(h2)];
                    const CGALMesh::Vertex_index b = oldToNewVertexMap[mesh.source(h1)];
                    const CGALMesh::Vertex_index c = oldToNewVertexMap[mesh.target(h1)];
                    const CGALMesh::Vertex_index d1 = edgeVertexMap.at(mesh.edge(h1));
                    const CGALMesh::Vertex_index d2 = edgeVertexMap.at(mesh.edge(h2));

                    newFaces.push_back({d1, c, d2});
                    newTetIds.push_back(tetId);

                    newFaces.push_back({d2, a, d1});
                    newTetIds.push_back(tetId);

                    newFaces.push_back({a, b, d1});
                    newTetIds.push_back(tetId);
                }

                else
                {
                }
            }




            for (int i = 0 ; i < newFaces.size() ; i++)
            {
                const std::array<CGALMesh::Vertex_index, 3> faceVertices = newFaces[i];

                auto newFace = newMesh.mesh.add_face(faceVertices);

                if (newFace == CGALMesh::null_face())
                {
                    throw std::runtime_error("Failed to add a new triangle to the msh.");
                }

                newMesh.tetId[newFace] = newTetIds[i];
            }

            

            // 4. Retriangulate affected faces
            //CGAL::Polygon_mesh_processing::triangulate_faces(facesToTriangulate, newMesh.mesh);

            if (false == CGAL::is_valid_polygon_mesh(newMesh.mesh))
            {
                throw std::runtime_error("New mesh is not valid.");
            }

            std::cout << "Previous number of faces : " << mesh.num_faces() << std::endl;
            std::cout << "New number of faces : " << newMesh.mesh.num_faces() << std::endl;

            return newMesh;
        }


        // Interpolate a point along an edge for a given isovalue
        //std::array<double,3> interpolateEdge(int a, int b, double iso) const
        //{
            //const double valA = edgeParam[a];
            //const double valB = edgeParam[b];

            //// Avoid division by zero (flat edge)
            //const double alpha = (valB != valA) ? (iso - valA) / (valB - valA) : 0.5;

            //const auto &posA = vertexCoordinates[a];
            //const auto &posB = vertexCoordinates[b];

            //return {
                //posA[0] + alpha * (posB[0] - posA[0]),
                    //posA[1] + alpha * (posB[1] - posA[1]),
                    //posA[2] + alpha * (posB[2] - posA[2])
            //};
        //}


        //Mesh to_cgal_mesh()
        //{
            //// 1. Build polygon soup
            //std::vector<CartesianPoint_3> points;
            //std::vector<std::vector<std::size_t>> polygons;
            //for (auto& p : vertexCoordinates)
            //{
                //points.push_back(CartesianPoint_3(p[0],p[1],p[2]));
            //}
            //for (auto& tri : triangles)
            //{
                //polygons.push_back({(std::size_t)tri[0], (std::size_t)tri[1], (std::size_t)tri[2]});
            //}

            //// 2. Remove duplicate points / faces
            //CGAL::Polygon_mesh_processing::merge_duplicate_points_in_polygon_soup(points, polygons);
            //CGAL::Polygon_mesh_processing::merge_duplicate_polygons_in_polygon_soup(points, polygons);
            //CGAL::Polygon_mesh_processing::orient_polygon_soup(points, polygons);

            //// 3. Build Surface_mesh
            //Mesh mesh;
            //CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, polygons, mesh);


            //return mesh;
        //}


        //void computeConnectedComponents(Mesh mesh)
        //{
            //std::vector<std::size_t> component(num_faces(mesh));
            //std::size_t num = CGAL::Polygon_mesh_processing::connected_components(mesh, CGAL::make_property_map(component));

            //std::cerr << "The number of components is ---------------------" << num << std::endl;
            //namespace PMP = CGAL::Polygon_mesh_processing;

            //std::cerr << "Vertices: " << num_vertices(mesh) << "\n";
            //std::cerr << "Faces: " << num_faces(mesh) << "\n";


            //std::cerr << "Is valid: "
                //<< CGAL::is_valid_polygon_mesh(mesh)
                //<< "\n";

            //std::vector< boost::graph_traits<Mesh>::halfedge_descriptor > borders;
            //PMP::border_halfedges(faces(mesh), mesh, std::back_inserter(borders));

            //std::size_t boundary_edge_count = borders.size() / 2; // each edge appears twice
            //std::cerr << "Boundary edges: " << boundary_edge_count << "\n";

        //}

        



        // Marching triangles
        //SurfaceMesh splitSingularTriangles(const double isovalue)
        //{
            //std::map<std::set<int>, std::array<double, 3>> triangleIntersectionPoints;

            //// Find the intersected edges as well as the points of intersection
            ////
            //for (const auto &triangle : this->triangles)
            //{
                //const int v0 = triangle[0];
                //const int v1 = triangle[1];
                //const int v2 = triangle[2];

                //const bool v0Inside = this->edgeParam[v0] < isovalue;
                //const bool v1Inside = this->edgeParam[v1] < isovalue;
                //const bool v2Inside = this->edgeParam[v2] < isovalue;

                //if (v0Inside != v1Inside)
                //{
                    //triangleIntersectionPoints[{v0, v1}] = interpolateEdge(v0, v1, isovalue); 
                //}
                //if (v1Inside != v2Inside)
                //{
                    //triangleIntersectionPoints[{v1, v2}] = interpolateEdge(v1, v2, isovalue); 
                //}
                //if (v2Inside != v0Inside)
                //{
                    //triangleIntersectionPoints[{v2, v0}] = interpolateEdge(v2, v0, isovalue); 
                //}
            //}

            //// Set up the new vertices
            ////
            //std::map<std::set<int>, int> triangleIntersectionIndices;

            //SurfaceMesh newMesh;
            //newMesh.vertexCoordinates = this->vertexCoordinates;
            //newMesh.edgeParam = this->edgeParam;
            //newMesh.isVertexSingular = this->isVertexSingular;

            //for (const auto &[edge, point] : triangleIntersectionPoints)
            //{
                //newMesh.vertexCoordinates.push_back(point);
                //newMesh.edgeParam.push_back(isovalue);
                //newMesh.isVertexSingular.push_back(true);
                //triangleIntersectionIndices[edge] = newMesh.vertexCoordinates.size()-1;
            //}


            //// Set up the new triangles
            //for (int i = 0 ; i < this->triangles.size() ; i++)
            //{
                //const auto &triangle = this->triangles[i];
                //const int &tetId = this->triangleTetId[i];

                //const int v0 = triangle[0];
                //const int v1 = triangle[1];
                //const int v2 = triangle[2];

                //const bool v0v1Intersected = triangleIntersectionIndices.contains({v0, v1});
                //const bool v1v2Intersected = triangleIntersectionIndices.contains({v1, v2});
                //const bool v2v0Intersected = triangleIntersectionIndices.contains({v2, v0});


                //// No intersected, skip this case
                //if (v0v1Intersected + v1v2Intersected + v2v0Intersected == 0)
                //{
                    //newMesh.triangles.push_back({
                            //v0, 
                            //v1, 
                            //v2, 
                            //});

                    //newMesh.triangleTetId.push_back(tetId);

                //}
                //else if (v0v1Intersected + v1v2Intersected + v2v0Intersected == 1)
                //{
                    //// Rotate so that vB is the intersected vertex is vB and intersected edges is vAvC
                    ////
                    //int vA, vB, vC;

                    //if (triangleIntersectionIndices.contains({v0, v1}))
                    //{
                        //vA = v1; vB = v2; vC = v0;
                    //}

                    //else if (triangleIntersectionIndices.contains({v1, v2}))
                    //{
                        //vA = v2; vB = v0; vC = v1;
                    //}

                    //else if (triangleIntersectionIndices.contains({v2, v0}))
                    //{
                        //vA = v0; vB = v1; vC = v2;
                    //}
                    //else
                    //{
                        //throw std::runtime_error("Impossible else case.");
                    //}

                    //newMesh.triangles.push_back({
                            //vA, 
                            //vB, 
                            //triangleIntersectionIndices.at({vA, vC}), 
                            //});

                    //newMesh.triangles.push_back({
                            //vB, 
                            //vC, 
                            //triangleIntersectionIndices.at({vA, vC}), 
                            //});

                    //newMesh.triangleTetId.push_back(tetId);
                    //newMesh.triangleTetId.push_back(tetId);
                //}
                //else if (v0v1Intersected + v1v2Intersected + v2v0Intersected == 2)
                //{
                    //// Rotate so that vB is the odd one out (between the two intersected edges)
                    ////
                    //int vA, vB, vC;

                    //if (triangleIntersectionIndices.contains({v0, v1}) && triangleIntersectionIndices.contains({v1, v2}))
                    //{
                        //vA = v0; vB = v1; vC = v2;
                    //}

                    //else if (triangleIntersectionIndices.contains({v1, v2}) && triangleIntersectionIndices.contains({v0, v2}))
                    //{
                        //vA = v1; vB = v2; vC = v0;
                    //}

                    //else if (triangleIntersectionIndices.contains({v0, v1}) && triangleIntersectionIndices.contains({v0, v2}))
                    //{
                        //vA = v2; vB = v0; vC = v1;
                    //}
                    //else
                    //{
                        //throw std::runtime_error("Impossible else case.");
                    //}


                    ////         vB
                    ////         /\
                    ////        /  \
                    ////       /____\
                    ////      /      \
                    ////     /________\
                    ////    vC        vA
                    ////

                    //newMesh.triangles.push_back({
                            //vB, 
                            //triangleIntersectionIndices.at({vB, vC}), 
                            //triangleIntersectionIndices.at({vA, vB}), 
                            //});

                    //newMesh.triangles.push_back({
                            //vA, 
                            //triangleIntersectionIndices.at({vA, vB}), 
                            //triangleIntersectionIndices.at({vB, vC}), 
                            //});


                    //newMesh.triangles.push_back({
                            //vC,
                            //vA,
                            //triangleIntersectionIndices.at({vB, vC}), 
                            //});

                    //newMesh.triangleTetId.push_back(tetId);
                    //newMesh.triangleTetId.push_back(tetId);
                    //newMesh.triangleTetId.push_back(tetId);

                //}
                //else if (v0v1Intersected + v1v2Intersected + v2v0Intersected == 3)
                //{
                    //throw std::runtime_error("Degenerate triangle detected in remeshing.");
                //}
                //else
                //{
                    //throw std::runtime_error("Impossible else case.");

                //}
            //}


            //return newMesh;
        //}


};
