#pragma once

#include <vector>

#include "./TetMesh.h"
#include "./FiberPoint.h"
#include "./ReebSpace2.h"
#include "./Arrangement.h"
#include "./FiberSurface.h"

namespace fiber
{
    std::vector<FiberPoint> computeLabeledFiber(TetMesh &, Arrangement &, ReebSpace2 &, std::array<double, 2>, const std::set<int> &);

    std::vector<FiberPoint> growSeedSet(const TetMesh &tetMesh, Arrangement &arrangement, ReebSpace2 &reebSpace, const std::array<double, 2> &fiberPoint, const std::vector<std::pair<int, int>> &fiberSeeds);

};
