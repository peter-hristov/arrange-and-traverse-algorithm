#include "./CGALTypedefs.h"

#include <cstddef>
#include <filesystem>

#include <GL/glut.h>
#include <QApplication>

#include "./io.h"
#include "./Timer.h"
#include "./ReebSpace.h"
#include "./ReebSpace2.h"
#include "./Data.h"
#include "./Arrangement.h"
#include "./utility/CLI11.hpp"
#include "./TracerVisualiserWindow.h"
#include "./ReebSpace2.h"

using namespace std;

int main(int argc, char* argv[])
{
    CLI::App cliApp("Reeb Space Fiber Visualiser");

    string filename;
    cliApp.add_option("--file, -f", filename, "Input data filename. Has to be either .txt of .vti.")->required();

    bool performanceRun = false;
    cliApp.add_flag("--performanceRun, -p", performanceRun, "Only compute the Reeb space, no graphics..");

    bool discardFiberSeeds = false;
    cliApp.add_flag("--discardPreimageGraphs, -d", discardFiberSeeds, "Discard the seeds for generating fibers based on sheets, discard to save a bit of memory (not too much).");

    float perturbationEpsilon = 1e-2;
    cliApp.add_option("--epsilon, -e", perturbationEpsilon, "Strength of the numerial perturbation in the range [-e, e].");

    string outputSheetPolygonsFilename;
    cliApp.add_option("--outputSheetPolygons, -o", outputSheetPolygonsFilename, "Filename where to output the coordinates of the polygons that represent each sheet.");

    int fiberSampling = 1;
    cliApp.add_option("--fiberSampling, -s", fiberSampling, "When saving fibers per component, how many do we save. Default is to save the centroid, otherwise sample along the boundary.");

    int sheetOutputCount = 10;
    cliApp.add_option("--sheetOutputCount", sheetOutputCount, "How many sheets to sample for automatic feature extraction.");

    string outputSheetFibersFolder;
    cliApp.add_option("--outputSheetFibersFolder", outputSheetFibersFolder, "Folder in which to ouput fiber for each sheet.");

    //string outputFibersFilename = "./fibers.vtp";
    //cliApp.add_option("--outputFibers", outputSheetPolygonsFilename, "Filename where to save the visible fiber components. Must be .vtp");

    CLI11_PARSE(cliApp, argc, argv);

    // For convenience
    if (performanceRun == true)
    {
        discardFiberSeeds = true;
    }

    // Read, perturb and sort the indices of the vertices lexicographically (by their range position).
    TetMesh tetMesh;
    try
    {
        Timer::start();
        tetMesh = io::readData(filename);
        Timer::stop("Reading input data                     :");
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error: " << e.what() << '\n';
        return 1;
    }
    Timer::start();
    tetMesh.perturbRangeValues(perturbationEpsilon);
    Timer::stop("Perturbing range values                :");

    Timer::start();
    tetMesh.sortVertices();
    Timer::stop("Sorting range points                   :");

    Timer::start();
    tetMesh.computeBoundingBoxes();
    Timer::stop("Computing bounding boxes               :");

    Timer::start();
    tetMesh.computeCombinatorialStructure();
    Timer::stop("Computing edges, triangles and tets    :");

    Timer::start();
    tetMesh.computeUpperLowerLinkAndStar();
    Timer::stop("Computing upper/lower links and stars  :");

    Timer::start();
    tetMesh.computeSingularEdgeTypes();
    Timer::stop("Computing singular edges               :");

    //Timer::start();

    Timer::start();
    Arrangement singularArrangement;
    singularArrangement.computeArrangement(tetMesh, Arrangement::SegmentMode::UseSingularSegments);
    Timer::stop("Singular Arrangement                   :");





    Timer::start();
    singularArrangement.checkInitialAssumptions(tetMesh);
    Timer::stop("Making sure all faces are simple       :");

    std::cout << "\n------------------------------------------------------------- A FRESH START.\n";

    std::vector<Segment_2> segments;
    for (const auto &[edge, type] : tetMesh.edgeSingularTypes) 
    {
        if (type == -1)
        {
            segments.push_back(Segment_2(singularArrangement.arrangementPoints[edge[0]], singularArrangement.arrangementPoints[edge[1]]));
        }

    }

    Timer::start();
    CGAL::insert(singularArrangement.arr, segments.begin(), segments.end());

    std::cout << "We have added " << segments.size() << " new segments to the arrangement.\n";
    std::cout << "The arrangement size:"
        << "   |V| = " << singularArrangement.arr.number_of_vertices()
        << ",  |E| = " << singularArrangement.arr.number_of_edges()
        << ",  |F| = " << singularArrangement.arr.number_of_faces() << std::endl << std::endl;
    Timer::stop("Singular Arrangement 2                 :");

    tetMesh.regularEdgesNumber -= segments.size();
    tetMesh.singularEdgesNumber += segments.size();


    Timer::start();
    singularArrangement.checkInitialAssumptions(tetMesh);
    Timer::stop("Making sure all faces are simple       :");









    // New computation
    ReebSpace2 reebSpace2;

    Timer::start();
    reebSpace2.computeEdgeRegionSegments(tetMesh, singularArrangement);
    Timer::stop("Computed red/blud intersetions         :");

    Timer::start();
    reebSpace2.computeVertexRegionSegments(tetMesh, singularArrangement);
    Timer::stop("Computed vertex regions                :");

    Timer::start();
    reebSpace2.computeEdgeRegionMinusPlusTriangles(tetMesh, singularArrangement);
    Timer::stop("Edge regions plus/minus triangles      :");

    Timer::start();
    reebSpace2.computeEdgeCrossingMinusPlusTriangles(tetMesh, singularArrangement);
    Timer::stop("Edge crossing plus/minus triangles     :");

    Timer::start();
    reebSpace2.computeVertexRegionMinusPlusTriangles(tetMesh, singularArrangement);
    Timer::stop("Vertex regions plus/minus triangles    :");

    //Timer::start();
    //reebSpace2.unitTest(tetMesh, singularArrangement, arrangement);
    //Timer::stop("Geometric computation unit tests       :");


    //std::map<std::array<int, 2>, std::vector<int>> upperStarTriangles;
    //std::map<std::array<int, 2>, std::vector<int>> lowerStarTriangles;
    
    //for (const auto& [edge, triangles] : tetMesh.upperStarTriangles)
    //{
        //printf("\n----------------------------------------------\n");
        ////std::cout << "Edge edge [" << he->source()->point() << " -> " << he->target()->point() << "]";
        //printf("Edge [%f, %f] -> [%f, %f]\n", tetMesh.vertexCoordinatesF[edge[0]], tetMesh.vertexCoordinatesG[edge[0]], tetMesh.vertexCoordinatesF[edge[1]], tetMesh.vertexCoordinatesG[edge[1]]);
        //printf("Edge %d -> %d\n", edge[0], edge[1]);
        //printf("\n----------------------------------------------\n\n");

        //std::cout << "Upper star triangles:\n";
        //for (const int &triangleId : tetMesh.upperStarTriangles[edge])
        //{
            //io::printTriangle(tetMesh, triangleId);
        //}
        //std::cout << "Lower star triangles:\n";
        //for (const int &triangleId : tetMesh.lowerStarTriangles[edge])
        //{
            //io::printTriangle(tetMesh, triangleId);
        //}

    //}

    //for (auto he = singularArrangement.arr.halfedges_begin(); he != singularArrangement.arr.halfedges_end(); ++he)
    //{
        //printf("\n----------------------------------------------\n");
        //std::cout << "Half edge [" << he->source()->point() << " -> " << he->target()->point() << "]";
        //printf("\n----------------------------------------------\n\n");

        //std::cout << "Edge region minus triangles:\n";
        //for (const int &triangleId : reebSpace2.edgeRegionMinusTriangles[he])
        //{
            //io::printTriangle(tetMesh, triangleId);
        //}

        //std::cout << "\n\nEdge region plus triangles:\n";
        //for (const int &triangleId : reebSpace2.edgeRegionPlusTriangles[he])
        //{
            //io::printTriangle(tetMesh, triangleId);
        //}

        //std::cout << "\n\nEdge crossing minus triangles:\n";
        //for (const int &triangleId : reebSpace2.edgeCrossingMinusTriangles[he])
        //{
            //io::printTriangle(tetMesh, triangleId);
        //}

        //std::cout << "\n\nEdge crossing plus triangles:\n";
        //for (const int &triangleId : reebSpace2.edgeCrossingPlusTriangles[he])
        //{
            //io::printTriangle(tetMesh, triangleId);
        //}

        //std::cout << "\n\nVertex region minus triangles:\n";
        //for (const int &triangleId : reebSpace2.vertexRegionMinusTriangles[he])
        //{
            //io::printTriangle(tetMesh, triangleId);
        //}

        //std::cout << "\n\nVertex region plus triangles:\n";
        //for (const int &triangleId : reebSpace2.vertexRegionPlusTriangles[he])
        //{
            //io::printTriangle(tetMesh, triangleId);
        //}

    //}


    //std::cout << "Press Enter to continue...";
    //std::cin.get();

    Timer::start();
    reebSpace2.traverse(tetMesh, singularArrangement);
    Timer::stop("Computed singular traversal            :");

    int correspondenceGraphSize = 0;
    for (const auto &[faceHandle, correspondenceGraph] : reebSpace2.correspondenceGraph)
    {
        correspondenceGraphSize += correspondenceGraph.size();
    }

    std::cout << "The size of the corresponde graph is " << correspondenceGraphSize << std::endl;

    //return 0;



    Timer::start();
    Arrangement arrangement;
    //arrangement.computeArrangement(tetMesh, Arrangement::SegmentMode::UseAllSegments);
    Timer::stop("Arrangement                            :");

    Timer::start();
    //arrangement.computePointLocationDataStructure();
    Timer::stop("Arrangement search structure           :");


    Timer::start();
    ReebSpace reebSpace;
    //reebSpace.computeTraversal(tetMesh, arrangement, discardFiberSeeds);
    Timer::stop("Computed {G_F} and H                   :");


    Timer::start();
    //reebSpace2.unitTestComparePreimageGraphs(tetMesh, singularArrangement, arrangement, reebSpace);
    Timer::stop("Comparing preimage graphs              :");


    //return 0;

    

    std::cout << "Postprocessing..." << std::endl;
    Timer::start();
    //reebSpace.computeSheetGeometry(tetMesh, arrangement);
    //reebSpace.computeSheetArea(tetMesh, arrangement);
    //reebSpace.printTopSheets(tetMesh, arrangement, 20);
    Timer::stop("Computed RS(f) Postprocess             :");

    if (performanceRun == true)
    {
        return 0;
    }

    if (false == outputSheetPolygonsFilename.empty())
    {
        try
        {
            io::saveSheets(tetMesh, arrangement, reebSpace, outputSheetPolygonsFilename);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Error: " << e.what() << '\n';
            return 1;
        }
    }


    if (false == outputSheetFibersFolder.empty())
    {
        try
        {
            io::generatefFaceFibersForSheets(tetMesh, arrangement, reebSpace, sheetOutputCount, fiberSampling, outputSheetFibersFolder);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Error: " << e.what() << '\n';
            return 1;
        }
    }


    // Set up QT Application
    QApplication app(argc, argv);
    glutInit(&argc, argv);

    // Package all my data for visualisation
    Data data(tetMesh, arrangement, singularArrangement, reebSpace, reebSpace2);

    // Create the widget
    TracerVisualiserWindow* window = new TracerVisualiserWindow(NULL, data);
    window->setWindowTitle("Fiber Visualiser");

    // Make the window full screen by default
    window->showMaximized();

    // Show the label
    window->show();

    // start it running
    app.exec();

    // clean up
    delete window;

    // return to caller
    return 0;
} // main()
