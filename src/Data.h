#pragma once

#include <GL/gl.h>
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

#include "./FaceFiber.h"

#include <qpoint.h>
#include <qvector.h>

class Data
{
  public:
    Data() {}

    // Min/max range coordinates
    GLfloat minF, maxF;
    GLfloat minG, maxG;

    // Min/max domain coordinates
    GLfloat minX, maxX;
    GLfloat minY, maxY;
    GLfloat minZ, maxZ;

    int xDim, yDim, zDim;

    // Compute the min/max F, G and X, Y, Z coordinates
    void computeMinMaxRangeDomainCoordinates();

    // Which of the tets in the tetrahedra array contain a fiber
    std::vector<bool> tetsWithFibers;

    std::string longnameF, longnameG, units;

    QVector<QPointF> mousePoints;

    // Tetraheda given in vertex IDs
    std::vector<std::vector<size_t>> tetrahedra;

    // The coordinates of vertices in the domain
    std::vector<std::vector<GLfloat>> vertexDomainCoordinates;

    // The coordinates of vertices in the domain, split into two arrays
    std::vector<GLfloat> vertexCoordinatesF;
    std::vector<GLfloat> vertexCoordinatesG;

    // These are the points on the faces of actives tets that have fibers (two per tet)
    std::vector<FaceFiberPoint> faceFibers;

    void computeTetExitPoints(const float, const float, const std::vector<float> = {1,1,1});

    void readData(std::string);

    // Read data file from a txt file
    // Construct a statndard hexahedral grid (split into tets)
    void readDataGrid(const std::string);

    size_t trippleToIndex(const size_t, const size_t, const size_t);
    void addTetsForCube(const size_t, const size_t, const size_t);

};
