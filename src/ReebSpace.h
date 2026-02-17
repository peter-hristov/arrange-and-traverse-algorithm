#pragma once

#include <vector>

#include "./TetMesh.h"
#include "./Arrangement.h"

class ReebSpace
{
    public:
        // The connected components of the preimage graph for each face in the arrangement
        std::vector<DisjointSet<int>> preimageGraphs;

        // Each vertex in the correspondence graph is a pair of {faceId, componentId}
        // The components of the correspondence graph are the sheets of the Reeb space
        DisjointSet<std::pair<int, int>> correspondenceGraph;

        // For each faceID we have a number of triangle seeds, which are given in fiberComponentId
        std::vector<std::vector<std::pair<int, int>>> fiberSeeds;



        std::unordered_map<int, CartesianPolygon_2> sheetPolygon;
        std::map<int, double> sheetArea;
        // Which sheets are incomplete, or degenerate
        std::unordered_set<int> incompleteSheets;
        // Maps the IDs of the Reeb space sheets to consequitive integers, useful for colouring things
        std::unordered_map<int, int> sheetConsequitiveIndices;
        // Jacobi type computing from the Reeb space algorithm, for comparison with the Jacobi type computed from the mesh
        std::unordered_map<std::pair<int, int>, int, MyHash<std::pair<int, int>>> jacobiType;




        // Computes the correspondence between two faces given a half edge between them
        void computeTwinFacePreimageGraph(const TetMesh &tetMesh, const Arrangement &arrangement, const Arrangement_2::Halfedge_const_handle &);

        // Computes the correspondence between fiber components of two incided faces
        void determineCorrespondence(const TetMesh &tetMesh, const Arrangement &arrangement, const Arrangement_2::Halfedge_const_handle &);

        // Compute the preimage graphs Gi for each cell in the arrangement and their correspondence
        void computeTraversal(const TetMesh &tetMesh, const Arrangement &arrangement, const bool, const bool);

        // Compute the geometry of the sheets of the Reeb space
        void computeReebSpacePostprocess(const TetMesh &tetMesh, const Arrangement &arrangement);

        void computeSheetGeometry(const TetMesh &tetMesh, const Arrangement &arrangement);

        void computeSheetArea(const TetMesh &tetMesh, const Arrangement &arrangement);

        void printTopSheets(const TetMesh &tetMesh, const Arrangement &arrangement, const int &numberToPrint);

        //void countIntersectionsTypes(TetMesh &tetMesh, Arrangement &arrangement);


        //
        // Helper functions
        //
        std::pair<std::vector<int>, std::vector<int>> getMinusPlusTrianglesIndex(const TetMesh &tetMesh, const Arrangement &arrangement, const Arrangement_2::Halfedge_const_handle currentHalfEdge);
        void testTraverseArrangement(const Arrangement &arrangement);
};
