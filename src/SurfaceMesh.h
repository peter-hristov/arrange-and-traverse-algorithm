#pragma once

#include "./CGALTypedefs.h"

#include <gmp.h>
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

// TODO Copy and move constructors
class SurfaceMesh
{
    public:
        CGALMesh mesh;

        // Named keys for the maps defined on the simplices of the mesh

        // The edge param from TTK, based on "P. Klacansky, J. Tierny, H. Carr and Z. Geng, "Fast and Exact Fiber Surfaces for Tetrahedral Meshes," in IEEE Transactions on Visualization and Computer Graphics, vol. 23, no. 7, pp. 1782-1795, 1 July 2017, doi: 10.1109/TVCG.2016.2570215)"
        static constexpr const char* EDGE_PARAM_KEY    = "v:edgeParam";

        // Which tet this triangle is in
        static constexpr const char* TET_ID_KEY        = "f:tetId";

        // The sheet label of this triangle after segmentation
        static constexpr const char* SHEET_ID_KEY      = "f:sheetId";

        // Whether the edge is a singular fiber which segments the fiber surfaces
        static constexpr const char* IMPASSABLE_KEY    = "e:isImpassable";

        // Which connected components this triangle is in with respect to the impassable edges
        static constexpr const char* COMPONENT_ID_KEY  = "f:componentId";

        const double epsilon = 1e-10;


        SurfaceMesh() 
        { 
            mesh.add_property_map<CGALMesh::Vertex_index, double>(EDGE_PARAM_KEY,    -1.0);
            mesh.add_property_map<CGALMesh::Face_index,   int>   (TET_ID_KEY,        -1);
            mesh.add_property_map<CGALMesh::Face_index,   int>   (SHEET_ID_KEY,      -1);
            mesh.add_property_map<CGALMesh::Face_index,   int>   (COMPONENT_ID_KEY,  -1);
            mesh.add_property_map<CGALMesh::Edge_index,   bool>  (IMPASSABLE_KEY,    false);
        }

        // Create mesh from a from a triangle soup
        SurfaceMesh(const std::vector<std::array<double, 3>> &vertexCoordinates, const std::vector<std::array<int, 3>> &triangles, const std::vector<double> &vertexEdgePara, const std::vector<int> &tetId);
        


        // Helpers to get the maps of the mesh
        auto edgeParam()
        {
            auto r = mesh.property_map<CGALMesh::Vertex_index, double>(EDGE_PARAM_KEY);
            if (!r.has_value()) throw std::runtime_error("edgeParam not initialized");
            return r.value();
        }

        auto tetId()
        {
            auto r = mesh.property_map<CGALMesh::Face_index, int>(TET_ID_KEY);
            if (!r.has_value()) throw std::runtime_error("tetId not initialized");
            return r.value();
        }

        auto sheetId()
        {
            auto r = mesh.property_map<CGALMesh::Face_index, int>(SHEET_ID_KEY);
            if (!r.has_value()) throw std::runtime_error("sheetId not initialized");
            return r.value();
        }

        auto componentId()
        {
            auto r = mesh.property_map<CGALMesh::Face_index, int>(COMPONENT_ID_KEY);
            if (!r.has_value()) throw std::runtime_error("componentId not initialized");
            return r.value();
        }

        auto isImpassable()
        {
            auto r = mesh.property_map<CGALMesh::Edge_index, bool>(IMPASSABLE_KEY);
            if (!r.has_value()) throw std::runtime_error("isImpassable not initialized");
            return r.value();
        }

        auto edgeParam() const
        {
            auto r = mesh.property_map<CGALMesh::Vertex_index, double>(EDGE_PARAM_KEY);
            if (!r.has_value()) throw std::runtime_error("edgeParam not initialized");
            return r.value();
        }

        auto tetId() const
        {
            auto r = mesh.property_map<CGALMesh::Face_index, int>(TET_ID_KEY);
            if (!r.has_value()) throw std::runtime_error("tetId not initialized");
            return r.value();
        }

        auto sheetId() const
        {
            auto r = mesh.property_map<CGALMesh::Face_index, int>(SHEET_ID_KEY);
            if (!r.has_value()) throw std::runtime_error("sheetId not initialized");
            return r.value();
        }

        auto componentId() const
        {
            auto r = mesh.property_map<CGALMesh::Face_index, int>(COMPONENT_ID_KEY);
            if (!r.has_value()) throw std::runtime_error("componentId not initialized");
            return r.value();
        }

        auto isImpassable() const
        {
            auto r = mesh.property_map<CGALMesh::Edge_index, bool>(IMPASSABLE_KEY);
            if (!r.has_value()) throw std::runtime_error("isImpassable not initialized");
            return r.value();
        }

        void print();
        void printSheetHistogram(ReebSpace2 &reebSpace);

        // Subdivision
        //
        CartesianPoint_3 interpolate_vertex(CGALMesh::Vertex_index v0, CGALMesh::Vertex_index v1, double isovalue);
        std::pair<CGALMesh::Property_map<CGALMesh::Vertex_index, int>, std::vector<CGALMesh::Vertex_index>> getVertexColours(CGALMesh &cgalMesh, const double isovalue);
        void subdivideMesh(const std::vector<double> &isovalues);
        void subdivideMeshOnce(const double isovalue);
        void triangulateMesh();
        void repairMesh();
        void validateMesh();
        std::unordered_set<CGALMesh::Edge_index> getActiveEdges(const std::vector<CGALMesh::Vertex_index> &grayVertices, const CGALMesh::Property_map<CGALMesh::Vertex_index, int> &vertexColour);

        // Labeling
        //
        bool bfsComponentFromSeed(const TetMesh &tetMesh, const Arrangement &singularArrangement, const int seedTriangleId, const std::vector<int> &tetTriangleIds, const CartesianPoint &controlPoint, std::vector<bool> &visited);

        int findFiberPointComponent(const TetMesh &tetMesh, const Arrangement &singularArrangement, const std::vector<std::pair<int, int>> &fiberSeeds, const std::vector<int> &tetTriangleIds, const Segment_2 &controlSegment, const double pointAlpha);
        
        int computeTriangleSheetId(TetMesh &tetMesh, Arrangement &singularArrangement, ReebSpace2 &reebSpace, const CGALMesh::Face_index &triangle, const std::vector<std::tuple<K::FT, int, int>> &intersectedSegments, const Segment_2 &controlSegment);

        void computeTriangleSheets(TetMesh &tetMesh, Arrangement &singularArrangement, ReebSpace2 &reebSpace, const std::vector<std::tuple<K::FT, int, int>> &intersectedSegments, const Segment_2 &controlSegment);

};
