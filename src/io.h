#pragma once

#include "./CGALTypedefs.h"

#include <string>

#include "./TetMesh.h"
#include "./SurfaceMesh.h"
#include "./Arrangement.h"
#include "./ReebSpace.h"
#include "./ReebSpace2.h"

#include "./FiberPoint.h"

#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>


namespace io
{
    TetMesh readData(const std::string&);
    TetMesh readDataTxt(const std::string&);
    TetMesh readDataVtu(const std::string&);
    //SurfaceMesh readDataVtuTTK(const std::string &filename);
    
    SurfaceMesh readDataVtuTTK(const std::string &filename, double u1, double v1, double u2, double v2);
    SurfaceMesh computeFiberSurface(vtkSmartPointer<vtkUnstructuredGrid> mesh, double u1, double v1, double u2, double v2);

    vtkSmartPointer<vtkPolyData> readMolecule(const std::string& filename);

    SurfaceMesh readDataVtp(const std::string&);
    CGALMesh readCGALMesh(const std::string& filename);

    void saveReebSpace(const ReebSpace2 &, const std::string&); 
    ReebSpace2 loadReebSpace(const std::string& filename);

    void saveSheets(const TetMesh &tetMesh, const Arrangement &arrangement, const ReebSpace &reebSpace, const std::string &outputSheetPolygonsFilename);
    void saveSheets2(const TetMesh &tetMesh, const Arrangement &arrangement, ReebSpace2 &reebSpace, const std::string &outputSheetPolygonsFilename);
    void saveFibers(const std::vector<FiberPoint>&, const std::string&);

    void saveSheetGraph(ReebSpace2 &reebSpace, const std::string&);

    void saveFiberSurface(std::vector<SurfaceMesh> &mesh, const std::string& filename);

    vtkSmartPointer<vtkPolyData> buildFiberSurfacePolyData(SurfaceMesh& surfMesh);

    void writeImpassableEdgesToVTK(const SurfaceMesh& surfMesh, const std::string& filename);
    void saveSheetsFeatures(const TetMesh &tetMesh, const Arrangement &arrangement, ReebSpace2 &reebSpace, const std::string &outputSheetPolygonsFilename);

    void saveOriginalMesh(const std::string, vtkSmartPointer<vtkUnstructuredGrid>);

    std::vector<FiberPoint> generatefFaceFibersForSheet(const TetMesh &tetMesh, Arrangement &arrangement, ReebSpace &reebSpace, const int sheetId, const int numberOfFiberPoints);
    void generatefFaceFibersForSheets(const TetMesh &tetMesh, Arrangement &arrangement, ReebSpace &reebSpace, const int sheetOutputCount, const int numberOfFiberPoints, const std::string);

    void printTriangle(const TetMesh &tetMesh, const int &triangleId);
}
