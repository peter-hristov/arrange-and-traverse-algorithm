#pragma once

#include "./CGALTypedefs.h"

#include <queue>

#include "./TetMesh.h"
#include "./Arrangement.h"
#include "./ReebSpace.h"
#include "./FiberGraph.h"
#include "./DisjointSetSimple.h"

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
        //       a        .        b 
        //                .
        //
        // Indexed by half-edge ID
        // For a half-edge $ab$, where "\, _, /" represent singular segments and ... represent regular ones,
        // The first preimage graph is o1: at the source $a$, after all regular segments
        // The second preimage graph is o2: at the target $b$, before all regular segments
        std::vector<std::pair<FiberGraph, FiberGraph>> fiberGraphs;

        // Save one fiber graph per face, for visualisation after
        std::vector<FiberGraph> representativeFiberGraphs;

        // Save one fiber graph seed per face, for visualisation after
        // Each seed is <triangleId, componentId>
        std::vector<std::vector<std::pair<int, int>>> representativeFiberGraphSeeds;



        
        // This tells you whether two sheets share a boundary (a half-edge or a vertex).
        // These are the edges of the sheet graph
        std::set<std::set<int>> areSheetsConnected;


        std::vector<std::vector<std::array<float, 2>>> sheetBoundary;
        //int orderIndex = 0;

        std::map<int, std::vector<int>> trianglesPerSheet;
        int numberOfSheets;

        std::map<int, double> sheetArea;
        std::map<int, double> sheetAreaProportion;


        std::map<int, std::vector<std::vector<Halfedge_const_handle>>> sheetBoundaries;





        // <Geometric computation>
        //
        void computeEdgeRegionSegments(const TetMesh &tetMesh, Arrangement &singularArrangement);
        void computeEdgeRegionSegments2(const TetMesh &tetMesh, Arrangement &singularArrangement);
        void computeEdgeRegionSegments3(const TetMesh &tetMesh, Arrangement &singularArrangement);
        void computeVertexRegionSegments(const TetMesh &tetMesh, Arrangement &singularArrangement);

        // Geometry helper functions
        bool compareRegularSegments(const Halfedge_const_handle &halfEdge, const Point_2& b, const Point_2& c);
        bool doSegmentEndpointsOverlap(const Segment_2 &s1, const Segment_2 &s2);
        bool ifSegmentInHalfEdgeRegion(const Halfedge_const_handle &halfEdge, const Segment_2 &segment);
        Halfedge_const_handle getSegmentRegion(Vertex_const_handle &vertexHandle, const Segment_2 &segment);

        // Plus/Minus triangles for each region
        void determineEdgeRegionSegmentsOrientation(const TetMesh &tetMesh, Arrangement &singularArrangement);
        void determineEdgeCrossingSegmentsOriantation(const TetMesh &tetMesh, Arrangement &singularArrangement);
        void determineVertexRegionSegmentsOrientation(const TetMesh &tetMesh, Arrangement &singularArrangement);

        //
        // </Geometric computation>







        // <Combinatorial computation>
        //

        // Compute the first fiber graph of a half-edge $e$, assuming the fiber graphs of its twin $e->twin()$ are already computed.
        void seedFace(TetMesh &, const Halfedge_const_handle &);

        // Compute the fiber graphs of all half-edge of a given face
        void loopFace(TetMesh &, const Halfedge_const_handle &);

        // Traverse the faces of the singular arrangement to compute the Reeb space
        void traverse(TetMesh &, Arrangement &, const bool);


        //
        // </Combinatorial computation>


        // Additional functions
        //
        void computeSheets(Arrangement &singularArrangement);

        void computeSheetBoundaries(Arrangement &singularArrangement);
        std::vector<std::vector<std::array<double, 2>>> computeSheetControlPolygons(const int);


        FiberGraph computeFiberGraph(TetMesh &tetMesh, Arrangement &singularArrangement, std::array<double, 2> controlPoint);
        int computeFiberGraphReverse(TetMesh &tetMesh, Arrangement &singularArrangement, std::array<double, 2> controlPoint, std::set<int> initialTriangles);

        FiberGraph computeFiberGraph2(TetMesh &tetMesh, Arrangement &singularArrangement, std::array<double, 2> controlPoint);
        FiberGraph computeFiberGraph3(TetMesh &tetMesh, Arrangement &singularArrangement, const Segment_2 &controlSegment, const K::FT &pointAlpha, const std::vector<std::tuple<K::FT, int, int>> &);


        std::vector<std::pair<int, int>> computeSeedFibers(TetMesh &tetMesh, Arrangement &singularArrangement, std::array<double, 2> controlPoint, const std::set<int> &);
        std::vector<std::pair<int, int>> computeSeedFibersGivenLine(TetMesh &tetMesh, Arrangement &singularArrangement, const Segment_2 &controlSegment, const K::FT &pointAlpha, const std::vector<std::tuple<K::FT, int, int>> &);

        // Unit Tests, mostly depricated
        void unitTest(const TetMesh &tetMesh, Arrangement &singularArrangement, Arrangement &regularArrangement);
        bool unitTestCompareFiberGraphs(const TetMesh &tetMesh, Arrangement &singularArrangement, Arrangement &regularArrangement, ReebSpace &rs);
        bool areHalfEdgeRegionMapsEqual(const std::map<Halfedge_const_handle, std::set<int>>& a, const std::map<Halfedge_const_handle, std::set<int>>& b);


};










