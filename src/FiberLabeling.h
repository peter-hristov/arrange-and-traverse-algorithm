#pragma once

#include "./TetMesh.h"
#include "./Arrangement.h"
#include "./ReebSpace2.h"

namespace fiber::labeling
{
    std::vector<std::pair<int, int>> computeSeedFibers(TetMesh &tetMesh, Arrangement &singularArrangement, ReebSpace2 &reebSpace, const std::array<double, 2> &controlPoint, const std::set<int> &selectedSheetIds);
    std::vector<std::pair<int, int>> computeSeedFibersGivenLine(TetMesh &tetMesh, Arrangement &singularArrangement, ReebSpace2 &reebSpace, const Segment_2 &controlSegment, const K::FT &pointAlpha, const std::vector<std::tuple<K::FT, int, int>> &intersectedSegments);

    // This is mostly depricated because it is too expensive to store fiber graphs per face on \bar{A}
    FiberGraph computeFiberGraph(TetMesh &tetMesh, Arrangement &singularArrangement, ReebSpace2 &reebSpace, std::array<double, 2> controlPoint);
}
