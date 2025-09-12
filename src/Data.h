#pragma once

#include <GL/gl.h>
#include <unordered_map>
#ifdef __APPLE__
#include <OpenGL/glu.h>
#include <OpenGL/gl.h>
#else
#include <GL/glu.h>
#include <GL/gl.h>
#endif

#include <cassert>
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <qpoint.h>
#include <qvector.h>

#include "./FiberPoint.h"
#include "./TetMesh.h"
#include "./Arrangement.h"
#include "./ReebSpace.h"
#include "./ReebSpace2.h"

class Data
{
  public:

    // Ideally, I would like a move constructor, but there's some issues in computing the search structure for the arrangement, something is not moved.
    // For not, just pass by reference

    TetMesh &tetMesh;
    Arrangement &arrangement;
    Arrangement &singularArrangement;
    ReebSpace &reebSpace;
    ReebSpace2 &reebSpace2;

    Data(TetMesh& tm, Arrangement& a, Arrangement& sa, ReebSpace& rs, ReebSpace2 &rs2)
        : tetMesh(tm),
        arrangement(a),
        singularArrangement(sa),
        reebSpace(rs),
        reebSpace2(rs2)
    {}

    std::string outputFibersFile = "./fibers.vtp";
};
