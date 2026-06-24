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
            allVertexPoints.push_back(pProj);
            allDomainCoordinates.push_back(interpolateDomainCoordinate(segment, pProj, allDomainCoordinates[edge[0]], allDomainCoordinates[edge[1]]));

            const int newVertexId = static_cast<int>(allVertexPoints.size()) - 1;
            //topologyGraphEdges.push_back({edge[0], newVertexId});
            //topologyGraphEdges.push_back({newVertexId, edge[1]});

            edgeSubdivisionPointIndex[edge] = newVertexId;

            //std::cout << "Adding split edge " << edge[0] << " -> " << allVertexPoints.size() << std::endl;
            //std::cout << "Adding split edge " << allVertexPoints.size() << " -> " << edge[1] << std::endl;
        }
        else
        {
            //topologyGraphEdges.push_back({edge[0], edge[1]});

            //std::cout << "Adding edge " << edge[0] << " -> " << edge[1] << std::endl;

        }
    }
    Timer::stop("Subdividing edges                      :");



    // 5. Compute all the distances
    Timer::start();
    std::vector<double> allVertexPointsDistances(allVertexPoints.size());

    for (int vId = 0 ; vId < allVertexPoints.size() ; vId++)
    {
        allVertexPointsDistances[vId] = CGAL::to_double(CGAL::squared_distance(pExact, allVertexPoints[vId]));
    }
    Timer::stop("Computing distances                    :");

    // 6. Construct the edges of the topology graph
    Timer::start();
    for (int i = 0 ; i <  tetMesh.tetrahedra.size() ; i++)
    {
        const auto &tet = tetMesh.tetrahedra[i];

        // Save all the vertices and the distance to sort later
        std::vector<std::pair<double, int>> tetVertices;

        // Add all the tet vertices
        for (int a = 0 ; a < 4 ; a++)
        {
            tetVertices.push_back({allVertexPointsDistances[tet[a]], tet[a]});
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
                    tetVertices.push_back({allVertexPointsDistances[vertexId], vertexId});
                }
            }
        }

        std::sort(tetVertices.begin(), tetVertices.end());

        for (int i = 0 ; i < tetVertices.size() - 1 ; i++)
        {
            topologyGraphEdges.push_back({tetVertices[i].second, tetVertices[i+1].second});
        }
    }

    Timer::stop("Computing topology graph               :");

    writeTopologyGraphToVTP(allDomainCoordinates, allVertexPoints, topologyGraphEdges, allVertexPointsDistances, "topologyGraph.vtp");
}

void timt::writeTopologyGraphToVTP(
    const std::vector<std::array<float, 3>> &domainCoordinates,
    const std::vector<Point_2> &rangeCoordinates,
    const std::vector<std::pair<int, int>> &edges,
    const std::vector<double> &vertexDistances,
    const std::string &filename)
{
    // Points
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    points->SetNumberOfPoints(domainCoordinates.size());
    //for (size_t i = 0; i < domainCoordinates.size(); ++i)
    //{
        //const std::array<float, 3> &coord = domainCoordinates[i];
        //points->SetPoint(i, coord[0], coord[1], coord[2]);
    //}

    for (size_t i = 0; i < rangeCoordinates.size(); ++i)
    {
        const float u = CGAL::to_double(rangeCoordinates[i].x());
        const float v = CGAL::to_double(rangeCoordinates[i].y());
        points->SetPoint(i, u, v, 0);
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
