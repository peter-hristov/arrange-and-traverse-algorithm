#pragma once

#include "./CGALTypedefs.h"

#include <gmp.h>
#include <iostream>
#include <stack>
#include <string>
#include <unordered_map>
#include <unordered_set>

#include "./TetMesh.h"
#include "./ReebSpace2.h"
#include "./Arrangement.h"
#include "./FiberGraph.h"

// TODO Copy and move constructors
class FiberSurface
{
    public:
        CGALMesh mesh;

        // Used to decide which vertices are black/gray/white in the segmentation
        static constexpr double epsilon = 0.001;

        //
        // Named keys for the maps defined on the simplices of the mesh
        //

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

        // Factory method
        static FiberSurface constructSegmentedFiberSurface(TetMesh &, Arrangement &, ReebSpace2 &, const std::vector<std::array<double, 2>> &, const std::set<int> & = {});

        //
        // Helpers to get the maps of the mesh
        //
        CGALMesh::Property_map<CGALMesh::Vertex_index, double> edgeParam();
        CGALMesh::Property_map<CGALMesh::Vertex_index, double> edgeParam() const;
        CGALMesh::Property_map<CGALMesh::Face_index,   int>    tetId();
        CGALMesh::Property_map<CGALMesh::Face_index,   int>    tetId() const;
        CGALMesh::Property_map<CGALMesh::Face_index,   int>    sheetId();
        CGALMesh::Property_map<CGALMesh::Face_index,   int>    sheetId() const;
        CGALMesh::Property_map<CGALMesh::Face_index,   int>    componentId();
        CGALMesh::Property_map<CGALMesh::Face_index,   int>    componentId() const;
        CGALMesh::Property_map<CGALMesh::Edge_index,   bool>   isImpassable();
        CGALMesh::Property_map<CGALMesh::Edge_index,   bool>   isImpassable() const;


        //
        // Constructors
        //
        FiberSurface();
        FiberSurface(const std::vector<std::array<double, 3>> &vertexCoordinates, const std::vector<std::array<int, 3>> &triangles, const std::vector<double> &vertexEdgePara, const std::vector<int> &tetId);
        

        //
        // Print stuff
        //
        void print();
        void printSheetHistogram(ReebSpace2 &reebSpace);



        //
        // Subdivision
        //

        // Remesh along the contour (singular fiber) for a given isovalue using the edgeParam map
        void remeshOnce(const double isovalue);

        // Remesh for multiple isovalues
        void remesh(const std::vector<double> &isovalues);

        // After remeshing many segmented face will not be triangles any more so we trinagulate
        void triangulate();

        // (Optional) repair mesh after the remeshing, expensive
        void repairMesh();

        // (Optional) make sure the mesh is valid and not degenerate
        void validateMesh();

        // While remeshing, interpolate along an edge at the isovalue to get the coordinates of the newly added vertex
        CartesianPoint_3 interpolateVertex(CGALMesh::Vertex_index v0, CGALMesh::Vertex_index v1, double isovalue);

        // While remeshing,set a colour to each vertex (black, gray, white; <, =, >)
        std::pair<CGALMesh::Property_map<CGALMesh::Vertex_index, int>, std::vector<CGALMesh::Vertex_index>> getVertexColours(CGALMesh &cgalMesh, const double isovalue);

        // While remeshing, find all the active edges in remeshing with a BFS from the gray vertices (singular points)
        std::unordered_set<CGALMesh::Edge_index> getActiveEdges(const std::vector<CGALMesh::Vertex_index> &grayVertices, const CGALMesh::Property_map<CGALMesh::Vertex_index, int> &vertexColour);



        //
        // Labeling
        //

        // The bfs procedure used in findFiberPointComponent
        bool bfsComponentFromSeed(const TetMesh &tetMesh, const Arrangement &singularArrangement, const int seedTriangleId, const std::vector<int> &tetTriangleIds, const CartesianPoint &controlPoint, std::vector<bool> &visited);

        // Once we have the seed set for fiber point, grow the components to see which one contains a tet which contains the triangle
        int findFiberPointComponent(const TetMesh &tetMesh, const Arrangement &singularArrangement, const std::vector<std::pair<int, int>> &fiberSeeds, const CartesianPoint &controlPoint, const int &tetId);
        
        // Label an individual triangle of a fiber surface based on their corresponding Reeb space sheet (do only after remeshing)
        int labelTriangle(TetMesh &tetMesh, Arrangement &singularArrangement, ReebSpace2 &reebSpace, const CGALMesh::Face_index &triangle, const std::vector<std::tuple<K::FT, int, int>> &intersectedSegments, const Segment_2 &controlSegment);

        // Label the triangles of a fiber surface based on their corresponding Reeb space sheet (do only after remeshing)
        void labelFiberSurface(TetMesh &tetMesh, Arrangement &singularArrangement, ReebSpace2 &reebSpace, const std::vector<std::tuple<K::FT, int, int>> &intersectedSegments, const Segment_2 &controlSegment);

        void filterTriangles(const std::set<int>&);

};
