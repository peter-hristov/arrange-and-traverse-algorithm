#include "./CGALTypedefs.h"

#include <cstddef>
#include <filesystem>
#include <fstream>

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
#include "./UnitTests.h"


using namespace std;

int main(int argc, char* argv[])
{
    CLI::App cliApp("Reeb Space Fiber Visualiser");

    string filename;
    cliApp.add_option("--file, -f", filename, "Input data filename. Has to be either .txt of .vti.")->required();

    bool performanceRun = false;
    cliApp.add_flag("--performanceRun, -p", performanceRun, "Only compute the Reeb space, no graphics..");

    bool unitTestPreimageGraphs = false;
    cliApp.add_flag("--unitTestPreimageGraphs, -u", unitTestPreimageGraphs, "Compute the preimage graphs with arrange and traverse and check if they are the same as with singular arrange and traverse.");

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

    //
    // TetMesh computation
    //

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



    //
    // Arrangement Computation
    //

    Timer::start();
    Arrangement singularArrangement;
    singularArrangement.computeArrangement(tetMesh, Arrangement::SegmentMode::UseSingularSegments);
    Timer::stop("Initial Singular Arrangement           :");

    // The initial computation may have nested faces, we need to connect them
    Timer::start();
    singularArrangement.connectNestedFaces(tetMesh);
    Timer::stop("Making sure all faces are simple       :");

    // Recompute arrangement if needed
    if (tetMesh.pseudoSingularEdgesNumber > 0)
    {
        singularArrangement = Arrangement();

        // Recompute the arrangement with enough new segments to avoid connect all nested faces
        Timer::start();
        singularArrangement.computeArrangement(tetMesh, Arrangement::SegmentMode::UseSingularSegments);
        Timer::stop("Singular Arrangement                   :");
    }

    Timer::start();
    singularArrangement.assignIndices();
    Timer::stop("Assigning indices to the arrangement   :");


    //
    // Reeb space computation.
    //
    ReebSpace2 reebSpace2;

    Timer::start();
    reebSpace2.computeEdgeRegionSegments(tetMesh, singularArrangement);
    Timer::stop("Computed red/blud intersetions         :");

    Timer::start();
    reebSpace2.computeEdgeRegionMinusPlusTriangles(tetMesh, singularArrangement);
    Timer::stop("Edge regions plus/minus triangles      :");

    Timer::start();
    reebSpace2.computeVertexRegionSegments(tetMesh, singularArrangement);
    Timer::stop("Computed vertex regions                :");

    Timer::start();
    reebSpace2.computeVertexRegionMinusPlusTriangles(tetMesh, singularArrangement);
    Timer::stop("Vertex regions plus/minus triangles    :");

    Timer::start();
    reebSpace2.computeEdgeCrossingMinusPlusTriangles(tetMesh, singularArrangement);
    Timer::stop("Edge crossing plus/minus triangles     :");

    Timer::start();
    reebSpace2.assignHalfEdgePseudoSingular(tetMesh, singularArrangement);
    Timer::stop("Assigning pseudosingular edges         :");

    Timer::start();
    reebSpace2.traverse(tetMesh, singularArrangement, false);
    Timer::stop("Computed singular traversal            :");

    //Timer::start();
    //reebSpace2.computeSheets(singularArrangement);
    //Timer::stop("Postprocessing                         :");

    int correspondenceGraphSize = 0;
    for (const auto &correspondenceGraph : reebSpace2.correspondenceGraph)
    {
        correspondenceGraphSize += correspondenceGraph.size();
    }


    // This is the old computation, keep these empty unless we want to unit test
    //
    ReebSpace reebSpace;
    Arrangement arrangement;

    if (unitTestPreimageGraphs)
    {
        Timer::start();
        arrangement.computeArrangement(tetMesh, Arrangement::SegmentMode::UseAllSegments);
        Timer::stop("Arrangement                            :");

        Timer::start();
        arrangement.computePointLocationDataStructure();
        Timer::stop("Arrangement search structure           :");


        Timer::start();
        reebSpace.computeTraversal(tetMesh, arrangement, discardFiberSeeds);
        Timer::stop("Computed {G_F} and H                   :");

        if (reebSpace2.numberOfSheets != reebSpace.correspondenceGraph.getComponentRepresentatives().size())
        {
            std::cerr << "----------------------------------------------------------------------------------------------------------------\n";
            std::cerr << "--------------------------------- NUMBER OF SHEETS IS NOT EQUAL!!!-----------------------------------------------\n";
            std::cerr << "---------------------The NEW number of sheets is " << reebSpace2.numberOfSheets << std::endl;
            std::cerr << "---------------------The OLD number of sheets is " << reebSpace.correspondenceGraph.getComponentRepresentatives().size() << std::endl;
            std::cerr << "----------------------------------------------------------------------------------------------------------------\n";
            return 1;
        }

        Timer::start();
        bool areSheetsEqual = unitTests::testAreSheetsIdentical(tetMesh, arrangement, singularArrangement, reebSpace, reebSpace2);
        Timer::stop("Determinig whether the sheet are equal :");

        if (false == areSheetsEqual)
        {
            std::cerr << "----------------------------------------------------------------------------------------------------------------\n";
            std::cerr << "--------------------------------- THE SHEETS ARE NOT EQUAL!!!--------------------------------------------------\n";
            std::cerr << "----------------------------------------------------------------------------------------------------------------\n";
            return 1;
        }


        //Timer::start();
        //bool arePreimageGraphsEqual = reebSpace2.unitTestComparePreimageGraphs(tetMesh, singularArrangement, arrangement, reebSpace);
        //Timer::stop("Comparing preimage graphs              :");

        //if (false == arePreimageGraphsEqual)
        //{
            //std::cerr << "----------------------------------------------------------------------------------------------------------------\n";
            //std::cerr << "--------------------------------- PREIMAGE GRAPHS NOT EQUAL!!!--------------------------------------------------\n";
            //std::cerr << "----------------------------------------------------------------------------------------------------------------\n";
            //return 1;
        //}


        //std::cout << "Postprocessing..." << std::endl;
        //Timer::start();
        //reebSpace.computeSheetGeometry(tetMesh, arrangement);
        //reebSpace.computeSheetArea(tetMesh, arrangement);
        //reebSpace.printTopSheets(tetMesh, arrangement, 20);
        //Timer::stop("Computed RS(f) Postprocess             :");

        return 0;
    }


    //return 0;

    if (performanceRun == true)
    {
        return 0;
    }

    if (false == outputSheetPolygonsFilename.empty())
    {
        try
        {
            printf("SAVING SHEETS--------------------");
            io::saveSheets(tetMesh, arrangement, reebSpace, outputSheetPolygonsFilename + ".old.vtp");
            io::saveSheets2(tetMesh, singularArrangement, reebSpace2, outputSheetPolygonsFilename);
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
