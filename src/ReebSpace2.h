#pragma once

#include "./CGALTypedefs.h"

#include <queue>

#include "./TetMesh.h"
#include "./Arrangement.h"
#include "./ReebSpace.h"
#include "./PreimageGraph.h"

#include "DisjointSetSimple.h"

class ReebSpace2
{
    public:

        // Index by half-edge ID, 
        // For each half-edge $e$, get its originating segment ID and whether it is left-to-right and right-to-left (direction of travel)
        std::vector<std::pair<int, bool>> edgeCrossingSegments;

        // Indexed by half-edge ID
        // For each helf-edge $e$, get the ordered list of regular segment IDs that intersect the "edge region" and whether they are left-to-right (direction of travel)
        // The "edge region" is defined as the interior of the half-edge and the order is linear along the half-edge
        std::vector<std::vector<std::pair<int, bool>>> edgeRegionSegments;

        // Indexed by half-edge ID
        // For each helf-edge $e$, get the ordered list of regular segment IDs that intersect the vertex region (have the same source as $e$) and whether they are left-to-right (direction of travel)
        // The vertex region is defined as the area between $e$ and $e->next$ around $e->source$, the order is CCW  as per CGAL convention
        std::vector<std::vector<std::pair<int, bool>>> vertexRegionSegments;

        // Index by singular face Id in the singular arrangement
        // For each singular face get the list of fiber graph class component (FGCC) Ids
        // If two FGCCs have the same ID, they are in the same sheet, if they have different IDs check in correspondenceGraphDS whether they are in the same sheets
        std::vector<std::vector<int>> correspondenceGraph;

        // Indexed by fiber graph class ID
        // For each fiber graph class component (FGCC), get the Reeb space sheet it belonds to
        DisjointSetSimple correspondenceGraphDS;


        //  \    .    .   .   .    .    /
        //   \   .   .    .    .   .   /
        //    \  .  .     .     .  .  /
        //     \ . .  o1  .  o2  . . /
        //      \.________.________./
        //       a        .         b 
        //                .
        //
        // Indexed by half-edge ID
        // For a half-edge $ab$, where "\, _, /" represent singular segments and ... represent regular ones,
        // The first preimage graph is o1: at the source $a$, after all regular segments
        // The second preimage graph is o2: at the target $b$, before all regular segments
        std::vector<std::pair<PreimageGraph, PreimageGraph>> preimageGraphs;

        std::vector<PreimageGraph> preimageGraphPerFace;
        
        std::set<std::set<int>> areSheetsConnected;


        std::vector<std::vector<std::array<float, 2>>> sheetBoundary;
        int orderIndex = 0;

        std::map<int, std::vector<int>> trianglesPerSheet;
        int numberOfSheets;

        std::map<int, double> sheetArea;
        std::map<int, double> sheetAreaProportion;


        // <Geometric computation>
        //
        void computeEdgeRegionSegments(const TetMesh &tetMesh, Arrangement &singularArrangement);
        void computeVertexRegionSegments(const TetMesh &tetMesh, Arrangement &singularArrangement);

        void computeEdgeRegionSegments2(const TetMesh &tetMesh, Arrangement &singularArrangement);

        // When sorting the segments around a vertex, compare two of them
        bool compareRegularSegments(const Halfedge_const_handle &halfEdge, const Point_2& b, const Point_2& c);
        bool doSegmentEndpointsOverlap(const Segment_2 &s1, const Segment_2 &s2);
        bool ifSegmentInHalfEdgeRegion(Arrangement_2::Halfedge_around_vertex_const_circulator &halfEdgeCirculator, const Segment_2 &segment);
        Halfedge_const_handle getSegmentRegion(Vertex_const_handle &vertexHandle, const Segment_2 &segment);

        // Plus/Minus triangles for each region
        void computeEdgeRegionMinusPlusTriangles(const TetMesh &tetMesh, Arrangement &singularArrangement);
        void computeEdgeCrossingMinusPlusTriangles(const TetMesh &tetMesh, Arrangement &singularArrangement);
        void computeVertexRegionMinusPlusTriangles(const TetMesh &tetMesh, Arrangement &singularArrangement);

        //
        // </Geometric computation>




        // Compute the actual Reeb space with "singular arrange and traverse"
        void seedFace(TetMesh &tetMesh, const Halfedge_const_handle &seedHalfEdge);
        void loopFace(TetMesh &tetMesh, const Halfedge_const_handle &seedHalfEdge);
        void traverse(TetMesh &tetMesh, Arrangement &singularArrangement, const bool);

        // Unit Tests
        void unitTest(const TetMesh &tetMesh, Arrangement &singularArrangement, Arrangement &regularArrangement);
        bool unitTestComparePreimageGraphs(const TetMesh &tetMesh, Arrangement &singularArrangement, Arrangement &regularArrangement, ReebSpace &rs);
        bool areHalfEdgeRegionMapsEqual(const std::map<Halfedge_const_handle, std::set<int>>& a, const std::map<Halfedge_const_handle, std::set<int>>& b);

        void computeSheets(Arrangement &singularArrangement);

};










