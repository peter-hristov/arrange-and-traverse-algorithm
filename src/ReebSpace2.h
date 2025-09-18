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
        //std::map<Halfedge_const_handle, std::pair<int, bool>> edgeCrossingSegments2;
        //std::map<Halfedge_const_handle, std::vector<std::pair<int, bool>>> edgeRegionSegments2;
        //std::map<Halfedge_const_handle, std::vector<std::pair<int, bool>>> vertexRegionSegments22;

        std::vector<std::pair<int, bool>> edgeCrossingSegments;
        std::vector<std::vector<std::pair<int, bool>>> vertexRegionSegments;
        std::vector<std::vector<std::pair<int, bool>>> edgeRegionSegments;

        std::vector<bool> isHalfEdgePseudoSingular;

        //
        // Preimage graphs and correspondence graph
        std::vector<std::vector<int>> correspondenceGraph;
        std::vector<PreimageGraph> preimageGraphs;

        //std::map<Face_const_handle, std::vector<int>> correspondenceGraph;
        //std::map<Face_const_handle, PreimageGraph> preimageGraphs;
        std::map<Halfedge_const_handle, std::pair<PreimageGraph, PreimageGraph>> preimageGraphsCached;


        //
        // Geometric computation
        void computeEdgeRegionSegments(const TetMesh &tetMesh, Arrangement &singularArrangement);
        void computeVertexRegionSegments(const TetMesh &tetMesh, Arrangement &singularArrangement);

        void computeEdgeRegionSegments2(const TetMesh &tetMesh, Arrangement &singularArrangement);

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

        // @TODO Move this to the TetMesh
        void assignHalfEdgePseudoSingular(const TetMesh &tetMesh, Arrangement &singularArrangement);


        //
        // Compute the actual REeb space
        void loopFace(TetMesh &tetMesh, const Halfedge_const_handle &seedHalfEdge, std::queue<Halfedge_const_handle> &traversalQueue, std::vector<bool> &visited, const bool);
        void traverse(TetMesh &tetMesh, Arrangement &singularArrangement, const bool);

        // 
        // Unit Tests
        void unitTest(const TetMesh &tetMesh, Arrangement &singularArrangement, Arrangement &regularArrangement);
        bool unitTestComparePreimageGraphs(const TetMesh &tetMesh, Arrangement &singularArrangement, Arrangement &regularArrangement, ReebSpace &rs);
        bool areHalfEdgeRegionMapsEqual(const std::map<Halfedge_const_handle, std::set<int>>& a, const std::map<Halfedge_const_handle, std::set<int>>& b);
};
