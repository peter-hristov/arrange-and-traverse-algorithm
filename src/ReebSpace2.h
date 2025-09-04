#pragma once

#include <queue>

#include "./TetMesh.h"
#include "./Arrangement.h"
#include "./ReebSpace.h"
#include "./PreimageGraph.h"

class ReebSpace2
{
    public:

        //
        // Geometric computation
        //
        std::map<Halfedge_const_handle, std::set<int>> edgeRegionSegments;
        std::map<Halfedge_const_handle, std::set<int>> vertexRegionSegments;

        void computeEdgeRegionSegments(const TetMesh &tetMesh, Arrangement &singularArrangement);
        void computeVertexRegionSegments(const TetMesh &tetMesh, Arrangement &singularArrangement);

        bool doSegmentEndpointsOverlap(const Segment_2 &s1, const Segment_2 &s2);
        bool ifSegmentInHalfEdgeRegion(Arrangement_2::Halfedge_around_vertex_const_circulator &halfEdgeCirculator, const Segment_2 &segment);
        Halfedge_const_handle getSegmentRegion(Vertex_const_handle &vertexHandle, const Segment_2 &segment);


        //void loopFace(const TetMesh &tetMesh, const Halfedge_const_handle &halfEdgeSeed);
        void loopFace(const TetMesh &tetMesh, const Halfedge_const_handle &seedHalfEdge, std::queue<Halfedge_const_handle> &traversalQueue, std::set<Face_const_handle> &visited);
        void traverse(const TetMesh &tetMesh, Arrangement &singularArrangement);


        std::map<Halfedge_const_handle, std::pair<PreimageGraph, PreimageGraph>> preimageGraphs;

        // Used for uint tests
        std::map<Halfedge_const_handle, std::pair<PreimageGraph, PreimageGraph>> preimageGraphsCached;


        // Make sure the image of the Jacobi set is simple (connected, no weird things going on).
        void checkInitialAssumptions(const TetMesh &tetMesh, Arrangement &singularArrangement);

        void unitTest(const TetMesh &tetMesh, Arrangement &singularArrangement, Arrangement &regularArrangement);
        void unitTestComparePreimageGraphs(const TetMesh &tetMesh, Arrangement &singularArrangement, Arrangement &regularArrangement, ReebSpace &rs);

        std::map<Halfedge_const_handle, std::vector<int>> edgeRegionMinusTriangles;
        std::map<Halfedge_const_handle, std::vector<int>> edgeRegionPlusTriangles;

        std::map<Halfedge_const_handle, std::vector<int>> edgeCrossingMinusTriangles;
        std::map<Halfedge_const_handle, std::vector<int>> edgeCrossingPlusTriangles;

        std::map<Halfedge_const_handle, std::vector<int>> vertexRegionMinusTriangles;
        std::map<Halfedge_const_handle, std::vector<int>> vertexRegionPlusTriangles;


        void computeEdgeRegionMinusPlusTriangles(const TetMesh &tetMesh, Arrangement &singularArrangement);
        void computeEdgeCrossingMinusPlusTriangles(const TetMesh &tetMesh, Arrangement &singularArrangement);
        void computeVertexRegionMinusPlusTriangles(const TetMesh &tetMesh, Arrangement &singularArrangement);

        bool areHalfEdgeRegionMapsEqual(const std::map<Halfedge_const_handle, std::set<int>>& a, const std::map<Halfedge_const_handle, std::set<int>>& b);


};
