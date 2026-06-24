#include "./src/CGALTypedefs.h"


#include <CGAL/number_utils.h>
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkLine.h>

#include "./src/TIMT.h"
#include "./src/Timer.h"
#include "./src/Fiber.h"


// Helper function to interpolate domain coordinates based on range coordiantes
std::array<float, 3> interpolateDomainCoordinate(
    const Segment_2 &segment,
    const Point_2 &pProj,
    const std::array<float, 3> &domA,
    const std::array<float, 3> &domB)
{
    const Point_2 &a = segment.source();
    const Point_2 &b = segment.target();

    const auto abx = b.x() - a.x();
    const auto aby = b.y() - a.y();
    const auto apx = pProj.x() - a.x();
    const auto apy = pProj.y() - a.y();
    const auto abLenSq = abx * abx + aby * aby;

    // t such that pProj = a + t * (b - a)
    const double t = CGAL::to_double((apx * abx + apy * aby) / abLenSq);

    std::array<float, 3> interpolated;
    for (int d = 0; d < 3; ++d)
    {
        interpolated[d] = static_cast<float>(
            (1.0 - t) * static_cast<double>(domA[d]) + t * static_cast<double>(domB[d]));
    }
    return interpolated;
}





timt::TopologyGraph timt::computeInexactTopologyGraph(TetMesh &tetMesh, const std::array<double, 2> &p)
{
    TopologyGraph tg;

    Timer::start();

    // 1. Set up range and domain coordinates
    //
    tg.vertexDomainCoordinates = tetMesh.vertexDomainCoordinates;
    tg.vertexRangeCoordinates = std::vector<Point_2>(tetMesh.vertexCoordinatesF.size());

    for (int i = 0 ; i < tetMesh.vertexCoordinatesF.size() ; i++)
    {
        const double u = tetMesh.vertexCoordinatesF[i];
        const double v = tetMesh.vertexCoordinatesG[i];

        tg.vertexRangeCoordinates[i] = Point_2(u, v);
    };
    Timer::stop("Converting range coordinates to exact  :");


    // 2. Compute all the distances
    //
    Timer::start();

    Point_2 pExact(p[0], p[1]);
    tg.vertexRangeDistances = std::vector<double>(tg.vertexRangeCoordinates.size());

    for (int vId = 0 ; vId < tg.vertexRangeCoordinates.size() ; vId++)
    {
        tg.vertexRangeDistances[vId] = CGAL::to_double(CGAL::squared_distance(pExact, tg.vertexRangeCoordinates[vId]));
    }
    Timer::stop("Computing vertex range distances       :");


    // 3. Construct the edges of the topology graph (using all edges)
    //
    Timer::start();

    for (int edgeId = 0 ; edgeId < tetMesh.edges.size() ; edgeId++)
    {
        const std::array<int, 2> &edge = tetMesh.edges[edgeId];

        tg.edges.insert({edge[0], edge[1]});
    }

    Timer::stop("Computing inexact topology graph       :");

    std::cout << "The inexact topology graph has " << tg.vertexRangeCoordinates.size() << " points and " << tg.edges.size() << " edges.\n";

    return tg;

}

timt::TopologyGraph timt::computeExactTopologyGraph(TetMesh &tetMesh, Arrangement &singularArrangement, ReebSpace2 &reebSpace, const std::array<double, 2> &p)
{
    TopologyGraph tg;

    Timer::start();
    tg.vertexDomainCoordinates = tetMesh.vertexDomainCoordinates;
    tg.vertexRangeCoordinates = std::vector<Point_2>(tetMesh.vertexCoordinatesF.size());

    // 1. Set up exact points
    for (int i = 0 ; i < tetMesh.vertexCoordinatesF.size() ; i++)
    {
        const double u = tetMesh.vertexCoordinatesF[i];
        const double v = tetMesh.vertexCoordinatesG[i];
        const Point_2 point(u, v);

        tg.vertexRangeCoordinates[i] = point;
    };
    Timer::stop("Setting up exact vertex points         :");


    // 2. Turn to segments
    std::vector<Segment_2> allSegments(tetMesh.edges.size());

    Timer::start();
    for (int edgeId = 0 ; edgeId < tetMesh.edges.size() ; edgeId++)
    {
        const std::array<int, 2> &edge = tetMesh.edges[edgeId];

        allSegments[edgeId] = Segment_2(tg.vertexRangeCoordinates[edge[0]], tg.vertexRangeCoordinates[edge[1]]);
    }
    Timer::stop("Setting up segments                    :");

    // 3. Turn p to an exact point
    Point_2 pExact(p[0], p[1]);


    std::map<std::array<int, 2>, int> edgeSubdivisionPointIndex;

    // 4. Subdivide edges
    Timer::start();
    for (int edgeId = 0 ; edgeId < tetMesh.edges.size() ; edgeId++)
    {
        const std::array<int, 2> &edge = tetMesh.edges[edgeId];

        const Segment_2 &segment = allSegments[edgeId];
        const Point_2 pProj = segment.supporting_line().projection(pExact);

        if (CGAL::collinear_are_ordered_along_line(segment.source(), pProj, segment.target()))
        {
            tg.vertexRangeCoordinates.push_back(pProj);
            tg.vertexDomainCoordinates.push_back(interpolateDomainCoordinate(segment, pProj, tg.vertexDomainCoordinates[edge[0]], tg.vertexDomainCoordinates[edge[1]]));

            const int newVertexId = static_cast<int>(tg.vertexRangeCoordinates.size()) - 1;
            edgeSubdivisionPointIndex[edge] = newVertexId;
        }
    }
    Timer::stop("Subdividing edges                      :");



    // 5. Compute all the distances
    //
    Timer::start();
    tg.vertexRangeDistances = std::vector<double>(tg.vertexRangeCoordinates.size());

    for (int vId = 0 ; vId < tg.vertexRangeCoordinates.size() ; vId++)
    {
        tg.vertexRangeDistances[vId] = CGAL::to_double(CGAL::squared_distance(pExact, tg.vertexRangeCoordinates[vId]));
    }
    Timer::stop("Computing vertex range distances       :");



    // Compute active tets
    //
    //
    //std::vector<std::set<int>> activeTrianglesPerComponent = Fiber::computeActiveTrianglesPerComponent(tetMesh, singularArrangement, reebSpace, p);
    std::vector<std::set<int>> activeTrianglesPerComponent;

    // This new component will be a new vertex in the topolgy graph, this tells us its index
    std::vector<int> componentTgraphIndex;

    for (int i = 0 ; i < activeTrianglesPerComponent.size() ; i++)
    {
        tg.vertexRangeCoordinates.push_back(pExact);
        tg.vertexRangeDistances.push_back(0.0);
        tg.vertexDomainCoordinates.push_back({-1, -1, -1});

        componentTgraphIndex.push_back(tg.vertexDomainCoordinates.size() - 1);
    }


    // 6. Construct the edges of the topology graph
    //
    Timer::start();

    for (int i = 0 ; i <  tetMesh.tetrahedra.size() ; i++)
    {
        const auto &tet = tetMesh.tetrahedra[i];

        // Save all the vertices and the distance to sort later
        std::vector<std::pair<double, int>> tetVertices;

        // Add all the tet vertices
        for (int a = 0 ; a < 4 ; a++)
        {
            tetVertices.push_back({tg.vertexRangeDistances[tet[a]], tet[a]});
        }

        // For all edges, all any (if there are) internal points
        for (int a = 0 ; a < 4 ; a++)
        {
            for (int b = a + 1 ; b < 4 ; b++)
            {
                int aIndex = tet[a];
                int bIndex = tet[b];

                if (aIndex > bIndex)
                {
                    std::swap(aIndex, bIndex);
                }

                std::array<int, 2> edge = {aIndex, bIndex};

                if (edgeSubdivisionPointIndex.contains(edge))
                {
                    const int vertexId = edgeSubdivisionPointIndex.at(edge);
                    tetVertices.push_back({tg.vertexRangeDistances[vertexId], vertexId});
                }
            }
        }

        // For all triangles, if any one of them is active, the tet is active, so add the zero-fiber component vertex
        for (int a = 0 ; a < 4 ; a++)
        {
            for (int b = a + 1 ; b < 4 ; b++)
            {
                for (int c = b + 1 ; c < 4 ; c++)
                {
                    const std::set<int> triangleSet{tet[a], tet[b], tet[c]};
                    const int triangleIndex = tetMesh.triangleIndices.at(triangleSet);

                    for (int id = 0 ; id < activeTrianglesPerComponent.size() ; id++)
                    {
                        if (activeTrianglesPerComponent[id].contains(triangleIndex))
                        {
                            tetVertices.push_back({0.0, componentTgraphIndex[id]});

                            break;
                        }
                    }
                }
            }
        }

        std::sort(tetVertices.begin(), tetVertices.end());

        for (int i = 0 ; i < tetVertices.size() - 1 ; i++)
        {
            tg.edges.insert({tetVertices[i].second, tetVertices[i+1].second});
        }
    }
    Timer::stop("Computing topology graph               :");


    std::cout << "The exact topology graph has " << tg.vertexRangeCoordinates.size() << " points and " << tg.edges.size() << " edges.\n";


    return tg;
}




void timt::writeTopologyGraphToVTP(
    const TopologyGraph &tg,
    const std::string &filename)
{

    // Points
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    points->SetNumberOfPoints(tg.vertexDomainCoordinates.size());

    //for (size_t i = 0; i < domainCoordinates.size(); ++i)
    //{
        //const std::array<float, 3> &coord = domainCoordinates[i];
        //points->SetPoint(i, coord[0], coord[1], coord[2]);
    //}

    for (size_t i = 0; i < tg.vertexRangeCoordinates.size(); ++i)
    {
        const float u = CGAL::to_double(tg.vertexRangeCoordinates[i].x());
        const float v = CGAL::to_double(tg.vertexRangeCoordinates[i].y());
        points->SetPoint(i, u, v, 0);
    }


    // Lines (edges)
    vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
    for (const auto &edge : tg.edges)
    {
        vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
        line->GetPointIds()->SetId(0, edge.first);
        line->GetPointIds()->SetId(1, edge.second);
        lines->InsertNextCell(line);
    }
    // Scalar array: distances (per-point)
    vtkSmartPointer<vtkDoubleArray> distanceArray = vtkSmartPointer<vtkDoubleArray>::New();
    distanceArray->SetName("Distance");
    distanceArray->SetNumberOfComponents(1);
    distanceArray->SetNumberOfTuples(tg.vertexRangeDistances.size());
    for (size_t i = 0; i < tg.vertexRangeDistances.size(); ++i)
    {
        distanceArray->SetValue(i, tg.vertexRangeDistances[i]);
    }
    // Assemble PolyData
    vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
    polyData->SetPoints(points);
    polyData->SetLines(lines);
    polyData->GetPointData()->AddArray(distanceArray);
    polyData->GetPointData()->SetActiveScalars("Distance");
    // Write
    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName(filename.c_str());
    writer->SetInputData(polyData);
    writer->Write();
}
