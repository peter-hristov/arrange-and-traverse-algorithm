#include <filesystem>

#include <vtkSmartPointer.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCell.h>
#include <vtkDataArray.h>
#include <vtkPointData.h>
#include <random>
#include <ranges>

#include <vtkCell.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyLine.h>
#include <vtkDataArray.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPolyLine.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>

//#include <vtkPolygon.h>
#include <vtkIntArray.h>
#include <vtkCellData.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkTriangleFilter.h>
#include <vtkTriangle.h>




#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLUnstructuredGridReader.h>

#include "./io.h"
#include "./TetMesh.h"
#include "./Fiber.h"

TetMesh io::readData(const std::string &filename)
{
    std::filesystem::path filePath(filename);
    
    if (!std::filesystem::exists(filePath)) 
    {
        throw std::runtime_error("File does not exist: " + filename);
    }

    std::string extension = filePath.extension().string();
    if (extension == ".vtu") 
    {
        return io::readDataVtu(filename);
    } 
    else if (extension == ".txt") 
    {
        return io::readDataTxt(filename);
    } 

    throw std::runtime_error("Unsupported file type: " + extension);
}

TetMesh io::readDataVtu(const std::string &filename)
{
    TetMesh tetMesh;

    // Read the VTU file
    vtkSmartPointer<vtkXMLUnstructuredGridReader> reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    reader->SetFileName(filename.c_str());

    reader->Update();

    vtkSmartPointer<vtkUnstructuredGrid> mesh = reader->GetOutput();
    if (!mesh)
    {
        throw std::runtime_error("Failed to get mesh output from the file: " + filename);
    }

    if (mesh->GetNumberOfPoints() == 0)
    {
        throw std::runtime_error("Mesh contains no points: " + filename);
    }

    if (mesh->GetNumberOfCells() == 0)
    {
        throw std::runtime_error("Mesh contains no cells: " + filename);
    }

    // Set deault names for the range axis
    tetMesh.longnameF = "f";
    tetMesh.longnameG = "g";

    int numVertices = mesh->GetPoints()->GetNumberOfPoints(); 
    int numTets = mesh->GetNumberOfCells();

    // Initialize all the data arrays
    tetMesh.vertexCoordinatesF = std::vector<double>(numVertices);
    tetMesh.vertexCoordinatesG = std::vector<double>(numVertices);
    tetMesh.tetrahedra = std::vector<std::array<int, 4>>(numTets);
    tetMesh.vertexDomainCoordinates = std::vector<std::array<float, 3>>(numVertices);

    // Print vertices
    vtkSmartPointer<vtkPoints> points = mesh->GetPoints();
    //std::cout << "Vertices:\n";
    for (vtkIdType i = 0; i < points->GetNumberOfPoints(); i++) {
        double p[3];
        points->GetPoint(i, p);
        //std::cout << "Vertex " << i << ": (" << p[0] << ", " << p[1] << ", " << p[2] << ")\n";

        tetMesh.vertexDomainCoordinates[i][0] = p[0];
        tetMesh.vertexDomainCoordinates[i][1] = p[1];
        tetMesh.vertexDomainCoordinates[i][2] = p[2];
    }

    // Print tetrahedra
    //std::cout << "\nTetrahedra:\n";
    for (vtkIdType i = 0; i < mesh->GetNumberOfCells(); i++) {
        vtkCell* cell = mesh->GetCell(i);
        if (cell->GetNumberOfPoints() == 4) { // Tetrahedron check
            //std::cout << "Tetrahedron " << i << ": ";
            for (vtkIdType j = 0; j < 4; j++) {
                //std::cout << cell->GetPointId(j) << " ";
                tetMesh.tetrahedra[i][j] = cell->GetPointId(j);
            }
            //std::cout << "\n";
        }
    }

    // Print vertex data arrays
    //std::cout << "\nVertex Data Arrays:\n";
    vtkPointData* pointData = mesh->GetPointData();

    assert(pointData->GetNumberOfArrays() >= 2);

    vtkDataArray* fDataArray = pointData->GetArray(1);
    vtkDataArray* gDataArray = pointData->GetArray(0);

    assert(fDataArray->GetNumberOfTuples() == numVertices);
    assert(gDataArray->GetNumberOfTuples() == numVertices);

    for (vtkIdType i = 0; i < fDataArray->GetNumberOfTuples(); i++) 
    {
        tetMesh.vertexCoordinatesF[i] = fDataArray->GetTuple1(i);
    }

    for (vtkIdType i = 0; i < gDataArray->GetNumberOfTuples(); i++) 
    {
        tetMesh.vertexCoordinatesG[i] = gDataArray->GetTuple1(i);
    }

    return tetMesh;
}


TetMesh io::readDataTxt(const std::string &filename)
{
    TetMesh tetMesh;
    
    // Set deault names for the range axis
    tetMesh.longnameF = "f";
    tetMesh.longnameG = "g";

    // Open data file
    std::ifstream dataFile (filename);
    if (false == dataFile.is_open()) 
    { 
        throw std::runtime_error("Could not open file: " + filename);
    }


    // Read in data in a string and skip the comments
    std::string rawStringData;
    std::string myline;
    while (dataFile) {
        std::getline (dataFile, myline);
        if (myline[0] == '#')
        {
            //std::cout << myline << '\n';
        }
        else
        {
            rawStringData += " " + myline;
        }
    }

    // Set up the inputstream from the string
    std::istringstream dataStream(rawStringData);

    // Read in the number of vertices and tets
    int numVertices, numTets;
    dataStream >> numVertices >> numTets;

    // Initialize all the data arrays
    tetMesh.vertexCoordinatesF = std::vector<double>(numVertices);
    tetMesh.vertexCoordinatesG = std::vector<double>(numVertices);
    tetMesh.tetrahedra = std::vector<std::array<int, 4>>(numTets);
    tetMesh.vertexDomainCoordinates = std::vector<std::array<float, 3>>(numVertices);

    // Read in the domain coordinates
    for  (int i = 0 ; i < numVertices ; i++)
    {
        dataStream >> tetMesh.vertexDomainCoordinates[i][0];
        dataStream >> tetMesh.vertexDomainCoordinates[i][1];
        dataStream >> tetMesh.vertexDomainCoordinates[i][2];
    }

    // Read in the range coordinates
    for  (int i = 0 ; i < numVertices ; i++)
    {
        dataStream >> tetMesh.vertexCoordinatesF[i];
        dataStream >> tetMesh.vertexCoordinatesG[i];
    }
    
    // Read in the tetrahedron configuration
    for  (int i = 0 ; i < numTets ; i++)
    {
        dataStream >> tetMesh.tetrahedra[i][0];
        dataStream >> tetMesh.tetrahedra[i][1];
        dataStream >> tetMesh.tetrahedra[i][2];
        dataStream >> tetMesh.tetrahedra[i][3];
    }

    return tetMesh;
}


void io::saveSheetsFeatures(const TetMesh &tetMesh,
                        const Arrangement &arrangement,
                        ReebSpace2 &reebSpace2,
                        const std::string &filename)
{

    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> triangles = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkIntArray> sheetIdsArray = vtkSmartPointer<vtkIntArray>::New();

    sheetIdsArray->SetName("SheetID");
    sheetIdsArray->SetNumberOfComponents(1);

    // Add all points (assume all unique, in order)
    for (const auto& v : tetMesh.vertexDomainCoordinates)
    {
        points->InsertNextPoint(v[0], v[1], v[2]);
    }

    for (const auto& [sheetId, triangleIds] : reebSpace2.trianglesPerSheet)
    {
        for (int triangleId : triangleIds)
        {
            const std::set<int>& triangle = tetMesh.triangles[triangleId];

            vtkIdType ids[3];
            int idx = 0;
            for (int vertexId : triangle)
                ids[idx++] = vertexId;

            vtkSmartPointer<vtkTriangle> tri = vtkSmartPointer<vtkTriangle>::New();
            tri->GetPointIds()->SetId(0, ids[0]);
            tri->GetPointIds()->SetId(1, ids[1]);
            tri->GetPointIds()->SetId(2, ids[2]);

            triangles->InsertNextCell(tri);
            sheetIdsArray->InsertNextValue(sheetId);
        }
    }

    // Build polydata
    polyData->SetPoints(points);
    polyData->SetPolys(triangles);

    // Attach sheet ID array to cell data
    polyData->GetCellData()->AddArray(sheetIdsArray);

    // Write to file
    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName(filename.c_str());
    writer->SetInputData(polyData);
    writer->Write();
}

void io::saveSheets2(const TetMesh &tetMesh,
                        const Arrangement &arrangement,
                        ReebSpace2 &reebSpace,
                        const std::string &outputSheetPolygonsFilename)
{
    auto points = vtkSmartPointer<vtkPoints>::New();
    auto polys = vtkSmartPointer<vtkCellArray>::New();
    auto sheetIds = vtkSmartPointer<vtkIntArray>::New();
    auto faceIds = vtkSmartPointer<vtkIntArray>::New();

    sheetIds->SetName("SheetId");
    faceIds->SetName("FaceId");


    std::map<Vertex_const_handle, int> vertexIdMap;

    for (auto v = arrangement.arr.vertices_begin(); v != arrangement.arr.vertices_end(); ++v)
    {
        const Point_2 &p = v->point();
        vtkIdType id = points->InsertNextPoint(CGAL::to_double(p.x()), CGAL::to_double(p.y()), 0.0);

        vertexIdMap[v] = vertexIdMap.size();
    }


    std::cout << "---------------------------------------- Outputting sheets\n";

    for (auto fit = arrangement.arr.faces_begin(); fit != arrangement.arr.faces_end(); ++fit)
    {
        if (fit->is_unbounded())
            continue;

        const int faceId = fit->data();

        std::vector<vtkIdType> ptIds;

        auto circ = fit->outer_ccb();
        auto start = circ;
        do
        {
            const int pointId = vertexIdMap.at(circ->target());
            ptIds.push_back(pointId);
            ++circ;
        } while (circ != start);

        for (const int &componentId : reebSpace.correspondenceGraph[faceId])
        {
            const int sheetId = reebSpace.correspondenceGraphDS.find(componentId);

            // Insert the polygon directly into vtkCellArray
            polys->InsertNextCell(ptIds.size(), ptIds.data());
            sheetIds->InsertNextValue(sheetId);
            faceIds->InsertNextValue(faceId);
        }
    }

    // Create polydata
    auto polyData = vtkSmartPointer<vtkPolyData>::New();
    polyData->SetPoints(points);
    polyData->SetPolys(polys);
    polyData->GetCellData()->AddArray(sheetIds);
    polyData->GetCellData()->AddArray(faceIds);

    // --- Triangulate polygons --- // very important, otherwise paraview will triangulate and can sometimes fill in missing polygons
    auto triangleFilter = vtkSmartPointer<vtkTriangleFilter>::New();
    triangleFilter->SetInputData(polyData);
    triangleFilter->Update();
    auto triangulatedPolyData = triangleFilter->GetOutput();

    // Write triangulated polydata
    auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName(outputSheetPolygonsFilename.c_str());
    writer->SetInputData(triangulatedPolyData);
    writer->SetDataModeToAscii(); // optional for debugging
    writer->Write();

    std::cout << "Saved " << polys->GetNumberOfCells() 
        << " polygons to " << outputSheetPolygonsFilename << std::endl;
}

void io::saveSheets(const TetMesh &tetMesh, const Arrangement &arrangement, const ReebSpace &reebSpace, const std::string &outputSheetPolygonsFilename)
{
    auto points = vtkSmartPointer<vtkPoints>::New();
    auto polys = vtkSmartPointer<vtkCellArray>::New();
    auto sheetIds = vtkSmartPointer<vtkIntArray>::New();

    sheetIds->SetName("SheetId");


    for (const auto &[sheetId, polygon] : reebSpace.sheetPolygon)
    {
        std::vector<vtkIdType> ptIds;

        for (int i = 0 ; i < polygon.size() ; i++)
        {
            const CartesianPoint &point = polygon[i];
            double u = point.x();
            double v = point.y();

            vtkIdType id = points->InsertNextPoint(u, v, 0.0);
            ptIds.push_back(id);
        }

        polys->InsertNextCell(ptIds.size(), ptIds.data());
        sheetIds->InsertNextValue(sheetId);
    }

    // Create polydata
    auto polyData = vtkSmartPointer<vtkPolyData>::New();
    polyData->SetPoints(points);
    polyData->SetPolys(polys);
    polyData->GetCellData()->AddArray(sheetIds);

    // --- Triangulate polygons --- // very important, otherwise paraview will triangulate and can sometimes fill in missing polygons
    auto triangleFilter = vtkSmartPointer<vtkTriangleFilter>::New();
    triangleFilter->SetInputData(polyData);
    triangleFilter->Update();
    auto triangulatedPolyData = triangleFilter->GetOutput();

    // Write triangulated polydata
    auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName(outputSheetPolygonsFilename.c_str());
    writer->SetInputData(triangulatedPolyData);
    writer->SetDataModeToAscii(); // optional for debugging
    writer->Write();

    std::cout << "Saved " << polys->GetNumberOfCells() 
        << " polygons to " << outputSheetPolygonsFilename << std::endl;

}



void io::saveFibers(const std::string &outputFile, const std::vector<FiberPoint> &fiberPoints)
{
    std::cout << "Saving fibers in " << outputFile << std::endl;
    //std::cout << "The fiber has size " << this->faceFibers.size() << std::endl;  

    // 1. Create the points
    auto points = vtkSmartPointer<vtkPoints>::New();
    auto idArray = vtkSmartPointer<vtkIntArray>::New();
    auto colourArray = vtkSmartPointer<vtkDoubleArray>::New();

    idArray->SetName("SheetId");
    idArray->SetNumberOfComponents(1);

    colourArray->SetName("Colour");
    colourArray->SetNumberOfComponents(3);

    for (const FiberPoint &p : fiberPoints)
    {
        points->InsertNextPoint(p.point.data());
        idArray->InsertNextValue(p.sheetId);
        colourArray->InsertNextTuple(p.colour.data());
    }

    // 3. Create the cells (wrap polyline in cell array)
    auto cells = vtkSmartPointer<vtkCellArray>::New();
    for (int i = 1 ; i < fiberPoints.size() ; i+=2)
    {
        if (fiberPoints[i-1].sheetId == fiberPoints[i].sheetId)
        {
            // One edge segment
            auto polyLine = vtkSmartPointer<vtkPolyLine>::New();
            polyLine->GetPointIds()->SetNumberOfIds(2);
            polyLine->GetPointIds()->SetId(0, i-1);
            polyLine->GetPointIds()->SetId(1, i);

            cells->InsertNextCell(polyLine);
        }
    }

    // 4. Create the polydata object
    auto polyData = vtkSmartPointer<vtkPolyData>::New();
    polyData->SetPoints(points);
    polyData->SetLines(cells);

    // 5. Attach the VertexID array to the point data
    polyData->GetPointData()->AddArray(idArray);
    polyData->GetPointData()->AddArray(colourArray);
    polyData->GetPointData()->SetScalars(colourArray);  // optional: for coloring

    // 6. Write to .vtp file (XML format)
    auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName(outputFile.c_str());
    writer->SetInputData(polyData);
    writer->Write();
}


std::vector<FiberPoint> io::generatefFaceFibersForSheet(const TetMesh &tetMesh, Arrangement &arrangement, ReebSpace &reebSpace, const int sheetId, const int numberOfFiberPoints)
{
    CartesianPolygon_2 &polygon = reebSpace.sheetPolygon.at(sheetId);

    if (polygon.size() == 0)
    {
        return {};
    }

    // Compute the controid so that we can pull all verties towards it
    CartesianPoint centroid = CGAL::centroid(polygon.vertices_begin(), polygon.vertices_end());

    // If need only one, get it at the center
    if (numberOfFiberPoints == 1)
    {
        const std::array<double, 2> fiberPoint = {(double)centroid.x(), (double)centroid.y()};
        const std::vector<FiberPoint> fiber = fiber::computeFiber(tetMesh, arrangement, reebSpace, fiberPoint, sheetId);
        printf("The fiber size is %d\n", fiber.size());
        return fiber;
    }

    std::vector<std::array<double, 2>> fiberPoints;


    // If we need more, sample along the boundary
    for (const CartesianPoint &point : polygon) 
    {
        // Get point from CGAL (and convert to double )
        double u = point.x();
        double v = point.y();

        // Interpolate closer to the centroid to make sure we are in the sheet ( if the sheet is "convex enough")
        const double alpha = 0.2;
        u = (1 - alpha) * u + alpha * centroid.x();
        v = (1 - alpha) * v + alpha * centroid.y();

        fiberPoints.push_back({u, v});
    }

    std::vector<FiberPoint> sheetFibers;

    // Calculate step size we only want some of the fiber points, not all
    double step = static_cast<double>(fiberPoints.size() - 1) / (numberOfFiberPoints - 1);

    for (int i = 0; i < numberOfFiberPoints; ++i) 
    {
        int index = static_cast<int>(i * step);

        const std::array<double, 2> fiberPoint = {fiberPoints[index][0], fiberPoints[index][1]};
        const std::vector<FiberPoint> fiber = fiber::computeFiber(tetMesh, arrangement, reebSpace, fiberPoint, sheetId);

        printf("The fiber size is %d\n", fiber.size());
        sheetFibers.insert(sheetFibers.end(), fiber.begin(), fiber.end());
    }

    return sheetFibers;
}

void io::generatefFaceFibersForSheets(const TetMesh &tetMesh, Arrangement &arrangement, ReebSpace &reebSpace, const int sheetOutputCount, const int numberOfFiberPoints, const std::string folderPath)
{
    namespace fs = std::filesystem;

    fs::path folderPathFs(folderPath);
    if (!fs::exists(folderPathFs)) 
    {
        fs::create_directory(folderPathFs);
    }

    for (const auto &[sheetId, colourId] : reebSpace.sheetConsequitiveIndices)
    {
        if (reebSpace.incompleteSheets.contains(sheetId))
        {
            printf("Skipping fiber %d, it's incomplete.",  sheetId);
        }

        if (colourId > sheetOutputCount || reebSpace.incompleteSheets.contains(sheetId))
        {
            continue;
        }

        std::cout << "-------------------------------------------------------------------------------------------- Generating fibers for sheet " << sheetId << "..." << std::endl;
        const std::vector<FiberPoint> sheetFibers = io::generatefFaceFibersForSheet(tetMesh, arrangement, reebSpace, sheetId, numberOfFiberPoints);

        //std::cout << "Saving fibers..." << std::endl;
        std::string outputFile = folderPathFs.string() + "/fibers_" + std::to_string(sheetId) + ".vtp";
        io::saveFibers(outputFile, sheetFibers);
    }
}

void io::printTriangle(const TetMesh &tetMesh, const int &triangleId)
{
    const std::set<int> triangle = tetMesh.triangles[triangleId];

    for (const int &vertexId : triangle)
    {
        printf("%d ", vertexId);

    }
    printf("\n");
    for (const int &vertexId : triangle)
    {
        printf("[%.1f, %.1f]\n", tetMesh.vertexCoordinatesF[vertexId], tetMesh.vertexCoordinatesG[vertexId]);
    }
    printf("--------\n");
}
