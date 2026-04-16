#pragma once

#include <vector>

#include "./TetMesh.h"
#include "./ReebSpace2.h"
#include "./Arrangement.h"
#include "./FiberSurface.h"

typedef std::pair<int, int> SeedPair;

class SeedSet
{
    public:

        std::vector<SeedPair> seedPairs;

        // Computing Labeled Seet sets
        static std::vector<std::pair<int, int>> computeFiberSeedSet(TetMesh &tetMesh, Arrangement &singularArrangement, ReebSpace2 &reebSpace, const std::array<double, 2> &controlPoint, const std::set<int> &selectedSheetIds);

        static std::vector<std::pair<int, int>> computeFiberSeedSetGivenLine(TetMesh &tetMesh, Arrangement &singularArrangement, ReebSpace2 &reebSpace, const Segment_2 &controlSegment, const K::FT &pointAlpha, const std::vector<std::tuple<K::FT, int, int>> &intersectedSegments);

};

