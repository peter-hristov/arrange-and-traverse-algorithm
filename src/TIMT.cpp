#include "./src/CGALTypedefs.h"


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

void timt::generateTopologyGraph(TetMesh &tetMesh, const std::array<double, 2> &p)
{

    Timer::start();
    std::vector<Point_2> allVertexPoints(tetMesh.vertexCoordinatesF.size());
    std::vector<std::array<float, 3>> allDomainCoordinates(tetMesh.vertexCoordinatesF.size());

    // 1. Set up exact points
    for (int i = 0 ; i < tetMesh.vertexCoordinatesF.size() ; i++)
    {
        const double u = tetMesh.vertexCoordinatesF[i];
        const double v = tetMesh.vertexCoordinatesG[i];
        const Point_2 point(u, v);

        allVertexPoints[i] = point;
        allDomainCoordinates[i] = tetMesh.vertexDomainCoordinates[i];
    };
    Timer::stop("Setting up exact vertex points         :");


    // 2. Turn to segments
    std::vector<Segment_2> allSegments(tetMesh.edges.size());

    Timer::start();
    for (int edgeId = 0 ; edgeId < tetMesh.edges.size() ; edgeId++)
    {
        const std::array<int, 2> &edge = tetMesh.edges[edgeId];

        allSegments[edgeId] = Segment_2(allVertexPoints[edge[0]], allVertexPoints[edge[1]]);
    }
    Timer::stop("Setting up segments                    :");

    // 3. Turn p to an exact point
    Point_2 pExact(p[0], p[1]);

    std::vector<std::pair<int, int>> topologyGraphEdges;

    // 4. Make the topology graph
    Timer::start();
    for (int edgeId = 0 ; edgeId < tetMesh.edges.size() ; edgeId++)
    {
        const std::array<int, 2> &edge = tetMesh.edges[edgeId];

        const Segment_2 &segment = allSegments[edgeId];
        const Point_2 pProj = segment.supporting_line().projection(pExact);


        if (CGAL::collinear_are_ordered_along_line(segment.source(), pProj, segment.target()))
        {
            allVertexPoints.push_back(pProj);
            allDomainCoordinates.push_back(interpolateDomainCoordinate(segment, pProj, allDomainCoordinates[edge[0]], allDomainCoordinates[edge[1]]));

            const int newVertexId = static_cast<int>(allVertexPoints.size()) - 1;
            topologyGraphEdges.push_back({edge[0], newVertexId});
            topologyGraphEdges.push_back({newVertexId, edge[1]});

            //std::cout << "Adding split edge " << edge[0] << " -> " << allVertexPoints.size() << std::endl;
            //std::cout << "Adding split edge " << allVertexPoints.size() << " -> " << edge[1] << std::endl;
        }
        else
        {
            topologyGraphEdges.push_back({edge[0], edge[1]});

            //std::cout << "Adding edge " << edge[0] << " -> " << edge[1] << std::endl;

        }
    }
    Timer::stop("Computing topology graph               :");

    // 5. Compute all the distances
    Timer::start();
    std::vector<double> allVertexPointsDistances(allVertexPoints.size());

    for (int vId = 0 ; vId < allVertexPoints.size() ; vId++)
    {
        allVertexPointsDistances[vId] = CGAL::to_double(CGAL::squared_distance(pExact, allVertexPoints[vId]));
    }



    Timer::stop("Computing distances                    :");

    writeTopologyGraphToVTP(allDomainCoordinates, topologyGraphEdges, allVertexPointsDistances, "topologyGraph.vtp");
}

void timt::writeTopologyGraphToVTP(
    const std::vector<std::array<float, 3>> &domainCoordinates,
    const std::vector<std::pair<int, int>> &edges,
    const std::vector<double> &vertexDistances,
    const std::string &filename)
{
    // Points
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    points->SetNumberOfPoints(domainCoordinates.size());
    for (size_t i = 0; i < domainCoordinates.size(); ++i)
    {
        const std::array<float, 3> &coord = domainCoordinates[i];
        points->SetPoint(i, coord[0], coord[1], coord[2]);
    }
    // Lines (edges)
    vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
    for (const auto &edge : edges)
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
    distanceArray->SetNumberOfTuples(vertexDistances.size());
    for (size_t i = 0; i < vertexDistances.size(); ++i)
    {
        distanceArray->SetValue(i, vertexDistances[i]);
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
