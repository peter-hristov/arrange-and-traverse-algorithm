#pragma once

#include <vector>

#include "./TetMesh.h"
#include "./ReebSpace2.h"
#include "./Arrangement.h"
#include "./FiberSurface.h"


namespace seeds
{
    typedef std::pair<int, int> SeedPair;
    typedef std::vector<SeedPair> SeedSet;

    // Computing Labeled Seet sets
    SeedSet computeFiberSeedSet(TetMesh &tetMesh, Arrangement &singularArrangement, ReebSpace2 &reebSpace, const std::array<double, 2> &controlPoint, const std::set<int> &selectedSheetIds);
    SeedSet computeFiberSeedSetGivenLine(TetMesh &tetMesh, Arrangement &singularArrangement, ReebSpace2 &reebSpace, const Segment_2 &controlSegment, const K::FT &pointAlpha, const std::vector<std::tuple<K::FT, int, int>> &intersectedSegments);
};

