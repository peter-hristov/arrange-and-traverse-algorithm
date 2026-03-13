#include <filesystem>

#include <vtkUnstructuredGrid.h>
#include <vtkCell.h>
#include <vtkDataArray.h>
#include <vtkPointData.h>
#include <random>
#include <ranges>

#include <vtkCell.h>
#include <vtkCellArray.h>
#include <vtkLine.h>
#include <vtkPoints.h>
#include <vtkPolyLine.h>
#include <vtkDataArray.h>
#include <vtkPointData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPolyLine.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>

//#include <vtkPolygon.h>
#include <vtkIntArray.h>
#include <vtkCellData.h>
#include <vtkTriangleFilter.h>
#include <vtkTriangle.h>


#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkStaticCleanPolyData.h>

#include<ttkFiberSurface.h>
#include<ttkRangePolygon.h>



#include "./io.h"
#include "./TetMesh.h"
#include "./Fiber.h"
#include "./SurfaceMesh.h"
#include "src/CGALTypedefs.h"

SurfaceMesh getSurfaceMesh(vtkPolyData* polyData)
{
    if (!polyData)
    {
        std::cerr << "Polydata is not valid.\n";
    }

    if (polyData->GetPoints()->GetNumberOfPoints() == 0)
    {
        return {};
    }

    // Clean up if it's a triangle soup, merge duplicated triangles
    vtkSmartPointer<vtkStaticCleanPolyData> cleaner = vtkSmartPointer<vtkStaticCleanPolyData>::New();
    cleaner->SetInputData(polyData);
    //cleaner->SetTolerance(1e-6);
    cleaner->SetTolerance(0);
    cleaner->Update();

    vtkSmartPointer<vtkPolyData> cleanedPolyData = cleaner->GetOutput();
    
    vtkPoints* points = cleanedPolyData->GetPoints();
    vtkCellArray* cells = cleanedPolyData->GetPolys();

    // Get first scalar array (vertex data)
    vtkDataArray* scalarEdgeParam = cleanedPolyData->GetPointData()->GetArray("EdgeParameterization");  
    vtkDataArray* scalarTetId = cleanedPolyData->GetCellData()->GetArray("TetIds");  

    if (!scalarEdgeParam)
    {
        std::cerr << "No point scalar data found.\n";
        return {};
    }
    if (!scalarTetId)
    {
        std::cerr << "No cell scalar data found.\n";
        return {};
    }

    //std::cout << "Before Number of points: " << polyData->GetPoints()->GetNumberOfPoints() << "\n";
    //std::cout << "Before Number of cells: " << polyData->GetNumberOfCells() << "\n";

    //std::cout << "Number of points: " << points->GetNumberOfPoints() << "\n";
    //std::cout << "Number of cells: " << cleanedPolyData->GetNumberOfCells() << "\n";

    std::vector<std::array<double, 3>> vertexCoordinates(points->GetNumberOfPoints()); 
    std::vector<double> edgeParam(points->GetNumberOfPoints()); 

    
    for (vtkIdType i = 0; i < points->GetNumberOfPoints(); ++i)
    {
        double p[3];
        points->GetPoint(i, p);
        vertexCoordinates[i] = {static_cast<double>(p[0]),
                                     static_cast<double>(p[1]),
                                     static_cast<double>(p[2])};

        edgeParam[i] = scalarEdgeParam->GetTuple1(i);
    }

    std::vector<std::array<int, 3>> triangles; 
    std::vector<int> triangleTetId;

    // --- read triangles ---
    vtkIdType npts = 0;
    const vtkIdType* pts = nullptr;

    cells->InitTraversal();
    while (cells->GetNextCell(npts, pts))
    {
        if (npts != 3)
            continue; // skip non-triangle cells

        triangles.push_back({
                static_cast<int>(pts[0]),
                static_cast<int>(pts[1]),
                static_cast<int>(pts[2])}
                );

        // Get the TetId for this cell
        vtkIdType cellId = triangles.size() - 1; 
        int tetId = static_cast<int>(scalarTetId->GetTuple1(cellId));
        triangleTetId.push_back(tetId);

    }


    return SurfaceMesh(vertexCoordinates, triangles, edgeParam, triangleTetId);


    // Manual merge
    //vtkPoints* points = merged->GetPoints();
    //vtkCellArray* cells = merged->GetPolys();

    //// Get first scalar array (vertex data)
    //vtkDataArray* scalarEdgeParam = merged->GetPointData()->GetArray("EdgeParameterization");  
    //vtkDataArray* scalarTetId = merged->GetCellData()->GetArray("TetIds");  

    //if (!scalarEdgeParam)
    //{
        //std::cerr << "No point scalar data found.\n";
        //return {};
    //}
    //if (!scalarTetId)
    //{
        //std::cerr << "No cell scalar data found.\n";
        //return {};
    //}


    //SurfaceMesh mesh;

    //// Point -> newIndex
    //std::map<std::array<double, 3>, int> uniquePoints;

    //// Old Index -> new Index
    //std::vector<int> oldIndices(points->GetNumberOfPoints());

    //for (vtkIdType i = 0; i < points->GetNumberOfPoints(); ++i)
    //{
        //double p[3];
        //points->GetPoint(i, p);
        //std::array<double, 3> point({p[0], p[1], p[2]});

        //if (false == uniquePoints.contains(point))
        //{
            //const int newIndex = uniquePoints.size();
            //uniquePoints[point] = newIndex;
            //oldIndices[i] = newIndex;
        //}
        //else
        //{

            //oldIndices[i] = uniquePoints.at(point);

        //}

        ////mesh.edgeParam[i] = scalarEdgeParam->GetTuple1(i);
    //}


    //mesh.edgeParam.resize(uniquePoints.size());
    //for (vtkIdType i = 0; i < points->GetNumberOfPoints(); ++i)
    //{
        //const int newIndex = oldIndices[i];
        //mesh.edgeParam[newIndex] = scalarEdgeParam->GetTuple1(i);
    //}

    //mesh.vertexCoordinates.resize(uniquePoints.size());
    //for (auto const& [pt, idx] : uniquePoints)
    //{
        //mesh.vertexCoordinates[idx] =
        //{
            //static_cast<float>(pt[0]),
            //static_cast<float>(pt[1]),
            //static_cast<float>(pt[2])
        //};
    //}

    //// --- read triangles ---
    //vtkIdType npts = 0;
    //const vtkIdType* pts = nullptr;

    //cells->InitTraversal();
    //while (cells->GetNextCell(npts, pts))
    //{
        //if (npts != 3)
        //{
            //throw std::runtime_error("Cell in a vtp is not a triangle: ");
        //}

        //mesh.triangles.push_back({
                //static_cast<int>(oldIndices[pts[0]]),
                //static_cast<int>(oldIndices[pts[1]]),
                //static_cast<int>(oldIndices[pts[2]])}
                //);

        //// Get the TetId for this cell
        //vtkIdType cellId = mesh.triangles.size() - 1; 
        //int tetId = static_cast<int>(scalarTetId->GetTuple1(cellId));
        //mesh.triangleTetId.push_back(tetId);

    //}

    //mesh.isVertexSingular = std::vector<bool>(mesh.vertexCoordinates.size(), false);

    //std::cout << "Number of points: " << points->GetNumberOfPoints() << "\n";
    //std::cout << "new Number of points: " << mesh.vertexCoordinates.size() << "\n";
    //std::cout << "Number of cells: " << polyData->GetNumberOfCells() << "\n";

    //return mesh;


    //const vtkIdType* pts = nullptr;

    //cells->InitTraversal();

    //while (cells->GetNextCell(npts, pts))
    //{
        //if (npts != 3)
        //{
            //std::cerr << "Non-triangle cell encountered.\n";
            //continue;
        //}

        //std::cout << "Triangle:\n";

        //for (int i = 0; i < 3; ++i)
        //{
            //double p[3];
            //points->GetPoint(pts[i], p);

            //double value = scalars->GetTuple1(pts[i]);

            //std::cout << "  Vertex " << i
                      //<< " | ID: " << pts[i]
                      //<< " | Pos: (" << p[0] << ", "
                                     //<< p[1] << ", "
                                     //<< p[2] << ")"
                      //<< " | Scalar: " << value
                      //<< "\n";

            //std::array<float, 3> point = {(float)p[0], (float)p[1], (float)p[2]};


        //}
    //}

}


CGALMesh io::readCGALMesh(const std::string& filename)
{
    // --- 1. Read VTP using VTK ---
    auto reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
    reader->SetFileName(filename.c_str());
    reader->Update();

    vtkPolyData* polyData = reader->GetOutput();

    // --- 2. Extract vertices ---
    std::vector<CartesianPoint_3> points;

    vtkPoints* vtkPts = polyData->GetPoints();
    points.reserve(vtkPts->GetNumberOfPoints());
    for (vtkIdType i = 0; i < vtkPts->GetNumberOfPoints(); ++i)
    {
        double p[3];
        vtkPts->GetPoint(i, p);
        points.emplace_back(p[0], p[1], p[2]);
    }

    // --- 3. Extract polygons (triangles) ---
    std::vector<std::vector<std::size_t>> polygons;

    vtkCellArray* cells = polyData->GetPolys();
    vtkIdType npts;
    const vtkIdType* ptIds;
    cells->InitTraversal();
    while (cells->GetNextCell(npts, ptIds))
    {
        if (npts != 3) continue; // skip non-triangles
        polygons.push_back({static_cast<std::size_t>(ptIds[0]),
                static_cast<std::size_t>(ptIds[1]),
                static_cast<std::size_t>(ptIds[2])});
    }

    std::cout << "Read VTP: " << points.size() << " vertices, "
        << polygons.size() << " triangles.\n";

    // --- 4. Remove duplicates ---
    //CGAL::Polygon_mesh_processing::merge_duplicate_points_in_polygon_soup(points, polygons);
    //CGAL::Polygon_mesh_processing::merge_duplicate_polygons_in_polygon_soup(points, polygons);

    //// --- 5. Convert polygon soup to Surface_mesh ---
    //Mesh mesh;
    //CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, polygons, mesh);



    CGAL::Polygon_mesh_processing::merge_duplicate_points_in_polygon_soup(points, polygons);
    CGAL::Polygon_mesh_processing::merge_duplicate_polygons_in_polygon_soup(points, polygons);

    // 3. Build Surface_mesh
    CGALMesh mesh;
    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, polygons, mesh);


    return mesh;
}

void io::saveFiberSurface(SurfaceMesh& surfMesh, const std::string& filename)
{
    const CGALMesh& mesh = surfMesh.mesh;

    // -------------------------------------------------------------------------
    // 1. Points
    // -------------------------------------------------------------------------
    auto vtkPts = vtkSmartPointer<vtkPoints>::New();
    vtkPts->SetNumberOfPoints(mesh.number_of_vertices());

    // Map CGAL vertex index -> VTK point id
    std::unordered_map<CGALMesh::Vertex_index, vtkIdType> vertexMap;
    vertexMap.reserve(mesh.number_of_vertices());

    // edgeParam lives on vertices
    auto vtkEdgeParam = vtkSmartPointer<vtkDoubleArray>::New();
    vtkEdgeParam->SetName("edgeParam");
    vtkEdgeParam->SetNumberOfComponents(1);
    vtkEdgeParam->SetNumberOfTuples(mesh.number_of_vertices());

    vtkIdType pid = 0;
    for (auto v : mesh.vertices())
    {
        const auto& pt = mesh.point(v);
        vtkPts->SetPoint(pid, CGAL::to_double(pt.x()),
                              CGAL::to_double(pt.y()),
                              CGAL::to_double(pt.z()));

        double ep = surfMesh.edgeParam[v];
        vtkEdgeParam->SetValue(pid, ep);

        vertexMap[v] = pid++;
    }

    // -------------------------------------------------------------------------
    // 2. Faces (triangles assumed; adapt for polygons if needed)
    // -------------------------------------------------------------------------
    auto vtkCells = vtkSmartPointer<vtkCellArray>::New();

    auto vtkTetId   = vtkSmartPointer<vtkIntArray>::New();
    vtkTetId->SetName("tetId");
    vtkTetId->SetNumberOfComponents(1);

    auto vtkSheetId = vtkSmartPointer<vtkIntArray>::New();
    vtkSheetId->SetName("sheetId");
    vtkSheetId->SetNumberOfComponents(1);

    auto vtkComponentId = vtkSmartPointer<vtkIntArray>::New();
    vtkComponentId->SetName("componentId");
    vtkComponentId->SetNumberOfComponents(1);

    auto vtkTriangleId = vtkSmartPointer<vtkIntArray>::New();
    vtkTriangleId->SetName("triangleId");
    vtkTriangleId->SetNumberOfComponents(1);

    for (auto f : mesh.faces())
    {
        // Collect vertices of this face
        std::vector<vtkIdType> ids;
        for (auto v : CGAL::vertices_around_face(mesh.halfedge(f), mesh))
            ids.push_back(vertexMap.at(v));

        vtkCells->InsertNextCell(static_cast<vtkIdType>(ids.size()), ids.data());

        int tid = surfMesh.tetId[f];
        int sid = surfMesh.sheetId[f];
        int cid = surfMesh.componentId[f];

        vtkTetId->InsertNextValue(tid);
        vtkSheetId->InsertNextValue(sid);
        vtkComponentId->InsertNextValue(cid);
        vtkTriangleId->InsertNextValue(f.idx());
    }

    // -------------------------------------------------------------------------
    // 3. Edge attribute: isImpassable
    //    VTK PolyData has no native edge arrays, so we store it as a field-data
    //    array (one value per edge, ordered by CGAL edge index).
    //    Alternatively you could split faces into lines; field data is simpler.
    // -------------------------------------------------------------------------
    auto vtkImpassable = vtkSmartPointer<vtkUnsignedCharArray>::New();
    vtkImpassable->SetName("isImpassable");
    vtkImpassable->SetNumberOfComponents(1);
    vtkImpassable->SetNumberOfTuples(mesh.number_of_edges());

    vtkIdType eid = 0;
    for (auto e : mesh.edges())
    {
        bool val = surfMesh.isImpassable[e];
        vtkImpassable->SetValue(eid++, val);
    }

    // -------------------------------------------------------------------------
    // 4. Assemble vtkPolyData
    // -------------------------------------------------------------------------
    auto polyData = vtkSmartPointer<vtkPolyData>::New();
    polyData->SetPoints(vtkPts);
    polyData->SetPolys(vtkCells);

    polyData->GetPointData()->AddArray(vtkEdgeParam);   // vertex attribute
    polyData->GetCellData()->AddArray(vtkTetId);        // face attributes
    polyData->GetCellData()->AddArray(vtkSheetId);
    polyData->GetCellData()->AddArray(vtkComponentId);
    polyData->GetCellData()->AddArray(vtkTriangleId);
    polyData->GetFieldData()->AddArray(vtkImpassable);  // edge attribute

    // -------------------------------------------------------------------------
    // 5. Write .vtp (XML PolyData)
    // -------------------------------------------------------------------------
    auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName(filename.c_str());
    writer->SetInputData(polyData);
    writer->SetDataModeToBinary();   // or SetDataModeToAscii()
    writer->Write();
}

void io::writeImpassableEdgesToVTK(const SurfaceMesh& surfMesh, const std::string& filename)
{
    const CGALMesh& mesh = surfMesh.mesh;

    // Points
    auto pts = vtkSmartPointer<vtkPoints>::New();
    pts->SetNumberOfPoints(mesh.number_of_vertices());

    std::unordered_map<CGALMesh::Vertex_index, vtkIdType> vertexMap;
    vertexMap.reserve(mesh.number_of_vertices());

    vtkIdType pid = 0;
    for (auto v : mesh.vertices())
    {
        const auto& pt = mesh.point(v);
        pts->SetPoint(pid, CGAL::to_double(pt.x()),
                          CGAL::to_double(pt.y()),
                          CGAL::to_double(pt.z()));
        vertexMap[v] = pid++;
    }

    // Impassable edges as line cells
    auto cells = vtkSmartPointer<vtkCellArray>::New();

    for (auto e : mesh.edges())
    {
        if (!surfMesh.isImpassable[e]) continue;

        auto h = mesh.halfedge(e);
        vtkIdType v0 = vertexMap.at(mesh.source(h));
        vtkIdType v1 = vertexMap.at(mesh.target(h));
        vtkIdType line[2] = {v0, v1};
        cells->InsertNextCell(2, line);
    }

    // Assemble
    auto polyData = vtkSmartPointer<vtkPolyData>::New();
    polyData->SetPoints(pts);
    polyData->SetLines(cells);

    // Write
    auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName(filename.c_str());
    writer->SetInputData(polyData);
    writer->SetDataModeToBinary();
    writer->Write();
}


SurfaceMesh io::readDataVtp(const std::string &filename)
{
    // Read VTP file
    vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
    reader->SetFileName(filename.c_str());
    reader->Update();

    vtkPolyData* polyData = reader->GetOutput();

    return getSurfaceMesh(polyData);
}


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

SurfaceMesh io::computeFiberSurface(vtkSmartPointer<vtkUnstructuredGrid> mesh, double u1, double v1, double u2, double v2)
{
    // This is correct I tested now
    std::string field1Name = mesh->GetPointData()->GetArrayName(0);
    std::string field2Name = mesh->GetPointData()->GetArrayName(1);

    // Create a polyline for the range polygon
    //
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    double start[3] = {v1, u1, 0.0};
    double end[3]   = {v2, u2, 0.0};
    vtkIdType idStart = points->InsertNextPoint(start);
    vtkIdType idEnd   = points->InsertNextPoint(end);
    vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
    line->GetPointIds()->SetId(0, idStart);
    line->GetPointIds()->SetId(1, idEnd);
    vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
    cells->InsertNextCell(line);
    vtkSmartPointer<vtkUnstructuredGrid> polyline = vtkSmartPointer<vtkUnstructuredGrid>::New();
    polyline->SetPoints(points);
    polyline->SetCells(VTK_LINE, cells);






// --- 4. Create first scalar array ---
    vtkSmartPointer<vtkDoubleArray> scalars1 = vtkSmartPointer<vtkDoubleArray>::New();
    scalars1->SetName(field1Name.c_str());
    scalars1->SetNumberOfComponents(1);
    scalars1->SetNumberOfTuples(points->GetNumberOfPoints());
    scalars1->SetValue(idStart, v1);
    scalars1->SetValue(idEnd, v2);

    polyline->GetPointData()->AddArray(scalars1);

    // --- 5. Create second scalar array ---
    vtkSmartPointer<vtkDoubleArray> scalars2 = vtkSmartPointer<vtkDoubleArray>::New();
    scalars2->SetName(field2Name.c_str());
    scalars2->SetNumberOfComponents(1);
    scalars2->SetNumberOfTuples(points->GetNumberOfPoints());
    scalars2->SetValue(idStart, u1);
    scalars2->SetValue(idEnd, u2);

    polyline->GetPointData()->AddArray(scalars2);

    // Optional: choose which is active for default visualization
    polyline->GetPointData()->SetScalars(scalars1); // or scalars2


    //std::string field1Name = mesh->GetPointData()->GetArrayName(0);
    //std::string field2Name = mesh->GetPointData()->GetArrayName(1);
    //int totalPoints = mesh->GetPointData()->GetArray(0)->GetSize();

    vtkSmartPointer<ttkFiberSurface> fiberSurface = vtkSmartPointer<ttkFiberSurface>::New();
    // --- Set multiple inputs ---
    fiberSurface->SetInputData(0, mesh);     // Input0
    fiberSurface->SetInputData(1, polyline);   // Input1
    fiberSurface->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, field1Name.c_str()); // scalar1
    fiberSurface->SetInputArrayToProcess(1, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, field2Name.c_str()); // scalar2
    fiberSurface->SetInputArrayToProcess(2, 1, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, field1Name.c_str()); // scalar1
    fiberSurface->SetInputArrayToProcess(3, 1, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, field2Name.c_str()); // scalar2
    fiberSurface->Update();

    //std::cout << "Field names are " << field1Name.c_str() << " and " << field2Name.c_str() << std::endl;

    vtkPolyData* fiberSurfMesh = vtkPolyData::SafeDownCast(fiberSurface->GetOutput());


    //std::cout << "The fiber surface has " << fiberSurfMesh->GetNumberOfCells() << " cells.\n";


    return getSurfaceMesh(fiberSurfMesh);
}

SurfaceMesh io::readDataVtuTTK(const std::string &filename, double u1, double v1, double u2, double v2)
{
    // Read the VTU file
    vtkSmartPointer<vtkXMLUnstructuredGridReader> reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    reader->SetFileName(filename.c_str());

    reader->Update();

    vtkSmartPointer<vtkUnstructuredGrid> mesh = reader->GetOutput();
    if (!mesh)
    {
        throw std::runtime_error("Failed to get mesh output from the file: " + filename);
    }

    return io::computeFiberSurface(mesh, u1, v1, u2, v2);
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

    tetMesh.originalMesh = mesh;

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





void io::saveSheetGraph(ReebSpace2 &reebSpace, const std::string &filename)
{
    const double minSize = 1.0;  // inches
    const double maxSize = 10.0;

    std::ofstream out(filename);
    if (!out.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }

    out << "graph G {\n";
    out << "    overlap=false;\n";
    out << "    splines=true;\n";
    out << "    node [shape=circle, style=filled, fillcolor=lightblue, color=blue];\n";

    for (const auto &[node, prop] : reebSpace.sheetAreaProportion)
    {
        double size = minSize + prop * (maxSize - minSize);

        out << "    " << node
            << " [width=" << size
            << ", height=" << size
            << ", fixedsize=true];\n";
    }

    // Write edges
    for (const auto &edgeSet : reebSpace.areSheetsConnected) {
        auto it = edgeSet.begin();
        int a = *it++;
        int b = *it;
        out << "    " << a << " -- " << b << ";\n";
    }

    out << "}\n";
    out.close();
}

vtkSmartPointer<vtkPolyData> io::readMolecule(const std::string& filename)
{
    vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
    reader->SetFileName(filename.c_str());
    reader->Update();

    vtkSmartPointer<vtkPolyData> polyData = reader->GetOutput();

    vtkPoints* points = polyData->GetPoints();
    vtkCellArray* lines = polyData->GetLines();
    vtkCellArray* verts = polyData->GetVerts();

    //printf("Read: %s\n", filename.c_str());
    //printf("  Points    : %lld\n", points->GetNumberOfPoints());
    //printf("  Lines     : %lld\n", lines->GetNumberOfCells());
    //printf("  Vertices  : %lld\n", verts->GetNumberOfCells());

    return polyData;
}

void io::saveReebSpace(const ReebSpace2 &reebSpace, const std::string& filename)
{
    std::ofstream out(filename, std::ios::binary);

    auto writeInt  = [&](int v)  { out.write(reinterpret_cast<const char*>(&v), sizeof(v)); };
    auto writeBool = [&](bool v) { out.write(reinterpret_cast<const char*>(&v), sizeof(v)); };

    // edgeCrossingSegments
    writeInt(reebSpace.edgeCrossingSegments.size());
    for (const auto& [id, dir] : reebSpace.edgeCrossingSegments) { writeInt(id); writeBool(dir); }

    // edgeRegionSegments
    writeInt(reebSpace.edgeRegionSegments.size());
    for (const auto& vec : reebSpace.edgeRegionSegments) {
        writeInt(vec.size());
        for (const auto& [id, dir] : vec) { writeInt(id); writeBool(dir); }
    }

    // vertexRegionSegments
    writeInt(reebSpace.vertexRegionSegments.size());
    for (const auto& vec : reebSpace.vertexRegionSegments) {
        writeInt(vec.size());
        for (const auto& [id, dir] : vec) { writeInt(id); writeBool(dir); }
    }

    // correspondenceGraph
    writeInt(reebSpace.correspondenceGraph.size());
    for (const auto& vec : reebSpace.correspondenceGraph) {
        writeInt(vec.size());
        for (int id : vec) writeInt(id);
    }

    // correspondenceGraphDS
    writeInt(reebSpace.correspondenceGraphDS.parent.size());
    for (int v : reebSpace.correspondenceGraphDS.parent) writeInt(v);
    for (int v : reebSpace.correspondenceGraphDS.rank)   writeInt(v);

    auto saveFiberGraph = [&](const FiberGraph& fg) {
        writeInt(fg.componentRoot.size());
        for (const auto& [k, v] : fg.componentRoot) { writeInt(k); writeInt(v); }
        writeInt(fg.componentRepresentative.size());
        for (const auto& [k, v] : fg.componentRepresentative) { writeInt(k); writeInt(v); }
    };

    // fiberGraphs
    writeInt(reebSpace.fiberGraphs.size());
    for (const auto& [fg1, fg2] : reebSpace.fiberGraphs) { saveFiberGraph(fg1); saveFiberGraph(fg2); }

    // representativeFiberGraphs
    writeInt(reebSpace.representativeFiberGraphs.size());
    for (const auto& fg : reebSpace.representativeFiberGraphs) saveFiberGraph(fg);
}

ReebSpace2 io::loadReebSpace(const std::string& filename)
{
    ReebSpace2 reebSpace;
    std::ifstream in(filename, std::ios::binary);

    auto readInt  = [&]() { int  v; in.read(reinterpret_cast<char*>(&v), sizeof(v)); return v; };
    auto readBool = [&]() { bool v; in.read(reinterpret_cast<char*>(&v), sizeof(v)); return v; };

    // edgeCrossingSegments
    int n = readInt();
    reebSpace.edgeCrossingSegments.resize(n);
    for (auto& [id, dir] : reebSpace.edgeCrossingSegments) { id = readInt(); dir = readBool(); }

    // edgeRegionSegments
    n = readInt();
    reebSpace.edgeRegionSegments.resize(n);
    for (auto& vec : reebSpace.edgeRegionSegments) {
        int m = readInt(); vec.resize(m);
        for (auto& [id, dir] : vec) { id = readInt(); dir = readBool(); }
    }

    // vertexRegionSegments
    n = readInt();
    reebSpace.vertexRegionSegments.resize(n);
    for (auto& vec : reebSpace.vertexRegionSegments) {
        int m = readInt(); vec.resize(m);
        for (auto& [id, dir] : vec) { id = readInt(); dir = readBool(); }
    }

    // correspondenceGraph
    n = readInt();
    reebSpace.correspondenceGraph.resize(n);
    for (auto& vec : reebSpace.correspondenceGraph) {
        int m = readInt(); vec.resize(m);
        for (int& id : vec) id = readInt();
    }

    // correspondenceGraphDS
    n = readInt();
    reebSpace.correspondenceGraphDS.parent.resize(n);
    reebSpace.correspondenceGraphDS.rank.resize(n);
    for (int& v : reebSpace.correspondenceGraphDS.parent) v = readInt();
    for (int& v : reebSpace.correspondenceGraphDS.rank)   v = readInt();


    auto loadFiberGraph = [&](FiberGraph& fg) {
        int n = readInt();
        fg.componentRoot.clear(); fg.componentRoot.reserve(n);
        for (int i = 0; i < n; i++) { int k = readInt(), v = readInt(); fg.componentRoot[k] = v; }
        n = readInt();
        fg.componentRepresentative.clear(); fg.componentRepresentative.reserve(n);
        for (int i = 0; i < n; i++) { int k = readInt(), v = readInt(); fg.componentRepresentative[k] = v; }
    };


    // fiberGraphs
    n = readInt();
    reebSpace.fiberGraphs.resize(n);
    for (auto& [fg1, fg2] : reebSpace.fiberGraphs) { loadFiberGraph(fg1); loadFiberGraph(fg2); }

    // representativeFiberGraphs
    n = readInt();
    reebSpace.representativeFiberGraphs.resize(n);
    for (auto& fg : reebSpace.representativeFiberGraphs) loadFiberGraph(fg);

    return reebSpace;
}

void io::saveOriginalMesh(const std::string filename, vtkSmartPointer<vtkUnstructuredGrid> originalMesh)
{
    auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName(filename.c_str());
    writer->SetInputData(originalMesh);
    writer->Write();
}
