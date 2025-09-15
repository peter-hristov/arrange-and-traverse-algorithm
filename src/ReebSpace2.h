#pragma once

#include "./CGALTypedefs.h"

#include <queue>

#include "./TetMesh.h"
#include "./Arrangement.h"
#include "./ReebSpace.h"
#include "./PreimageGraph.h"

class ReebSpace2
{
    public:

        //
        // Geometric computation, which segments intersect which others, in which order and whether their upper/lower triangles are flipped
        std::map<Halfedge_const_handle, std::vector<std::pair<int, bool>>> edgeRegionSegments;
        std::map<Halfedge_const_handle, std::vector<int>> vertexRegionSegments;

        // 
        // Plus/Minus Triangles for each region
        //std::map<Halfedge_const_handle, std::vector<std::vector<int>>> edgeRegionMinusTriangles;
        //std::map<Halfedge_const_handle, std::vector<std::vector<int>>> edgeRegionPlusTriangles;

        std::map<Halfedge_const_handle, std::vector<int>> edgeCrossingMinusTriangles;
        std::map<Halfedge_const_handle, std::vector<int>> edgeCrossingPlusTriangles;

        std::map<Halfedge_const_handle, std::vector<std::vector<int>>> vertexRegionMinusTriangles;
        std::map<Halfedge_const_handle, std::vector<std::vector<int>>> vertexRegionPlusTriangles;

        //
        // Preimage graphs and correspondence graph
        std::map<Face_const_handle, std::vector<int>> correspondenceGraph;
        std::map<Halfedge_const_handle, std::pair<PreimageGraph, PreimageGraph>> preimageGraphs;
        std::map<Halfedge_const_handle, std::pair<PreimageGraph, PreimageGraph>> preimageGraphsCached;


        //
        // Geometric computation
        void computeEdgeRegionSegments(const TetMesh &tetMesh, Arrangement &singularArrangement);
        void computeVertexRegionSegments(const TetMesh &tetMesh, Arrangement &singularArrangement);

        // When sorting the segments around a vertex, compare two of them
        bool compareRegularSegments(const Halfedge_const_handle &halfEdge, const Point_2& b, const Point_2& c);

        bool doSegmentEndpointsOverlap(const Segment_2 &s1, const Segment_2 &s2);
        bool ifSegmentInHalfEdgeRegion(Arrangement_2::Halfedge_around_vertex_const_circulator &halfEdgeCirculator, const Segment_2 &segment);
        Halfedge_const_handle getSegmentRegion(Vertex_const_handle &vertexHandle, const Segment_2 &segment);

        //
        // Plus/Minus triangles for each region
        void computeEdgeRegionMinusPlusTriangles(const TetMesh &tetMesh, Arrangement &singularArrangement);
        void computeEdgeCrossingMinusPlusTriangles(const TetMesh &tetMesh, Arrangement &singularArrangement);
        void computeVertexRegionMinusPlusTriangles(const TetMesh &tetMesh, Arrangement &singularArrangement);


        //
        // Compute the actual REeb space
        void loopFace(TetMesh &tetMesh, const Halfedge_const_handle &seedHalfEdge, std::queue<Halfedge_const_handle> &traversalQueue, std::set<Face_const_handle> &visited, Arrangement &singularArrangement);
        void traverse(TetMesh &tetMesh, Arrangement &singularArrangement);

        // 
        // Unit Tests
        void unitTest(const TetMesh &tetMesh, Arrangement &singularArrangement, Arrangement &regularArrangement);
        void unitTestComparePreimageGraphs(const TetMesh &tetMesh, Arrangement &singularArrangement, Arrangement &regularArrangement, ReebSpace &rs);
        bool areHalfEdgeRegionMapsEqual(const std::map<Halfedge_const_handle, std::set<int>>& a, const std::map<Halfedge_const_handle, std::set<int>>& b);
};
