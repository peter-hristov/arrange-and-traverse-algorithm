#pragma once

#include <vector>

#include <vtkPolyData.h>
#include <vtkSmartPointer.h>

#include "./TetMesh.h"
#include "./Arrangement.h"
#include "./ReebSpace.h"
#include "./ReebSpace2.h"
#include "./FiberSurface.h"
#include "./Fiber.h"


class Data
{
  public:

    // Ideally, I would like a move constructor, but there's some issues in computing the search structure for the arrangement, something is not moved.
    // For not, just pass by reference

    TetMesh &tetMesh;

    Arrangement &arrangement;
    ReebSpace &reebSpace;

    Arrangement &singularArrangement;
    ReebSpace2 &reebSpace2;

    Data(TetMesh& tm, Arrangement& a, Arrangement& sa, ReebSpace& rs, ReebSpace2 &rs2)
        : tetMesh(tm),
        arrangement(a),
        singularArrangement(sa),
        reebSpace(rs),
        reebSpace2(rs2)
    {}

    std::vector<Fiber> fibers;
    std::vector<FiberSurface> fiberSurfaces;
    std::vector<FiberSurface> featureSurfaces;


    vtkSmartPointer<vtkPolyData> molecule;

    std::set<int> selectedSheetIds;

    bool zeroAxis;
};
