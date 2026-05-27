#include "./CGALTypedefs.h"

#include <cstddef>
#include <filesystem>
#include <fstream>

#include <GL/glut.h>
#include <QApplication>
#include <string>

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
#include "./Performance.h"
#include "./Fiber.h"
#include "./LoadingBar.hpp"
#include "./TetMesh.h"
#include "./Augmentation.h"

using namespace std;

int main(int argc, char* argv[])
{
    CLI::App cliApp("Reeb Space Fiber Visualiser");

    string filename;
    cliApp.add_option("--file, -f", filename, "Input data filename. Has to be either .txt of .vti.")->required();

    string moleculeFilename;
    cliApp.add_option("--molecule, -m", moleculeFilename, "Input data filename. Has to be either .txt of .vti.");

    bool performanceRun = false;
    cliApp.add_flag("--performanceRun, -p", performanceRun, "Only compute the Reeb space, no graphics..");

    bool unitTestSheets = false;
    cliApp.add_flag("--unitTestSheets, -u", unitTestSheets, "Compute the fiber graphs with arrange and traverse and check if they are the same as with singular arrange and traverse.");

    bool unitTestFiberGraphs = false;
    cliApp.add_flag("--unitTestFiberGraphs, -U", unitTestFiberGraphs, "Compute the sheets with arrange and traverse and check if they are the same as with singular arrange and traverse.");

    bool discardFiberSeeds = false;
    cliApp.add_flag("--discardPreimageGraphs, -d", discardFiberSeeds, "Discard the seeds for generating fibers based on sheets, discard to save a bit of memory (not too much).");

    bool headless = false;
    cliApp.add_flag("--headless", headless, "Run without the graphical interface.");

    bool addZeroAxis = false;
    cliApp.add_flag("--zeroAxis", addZeroAxis, "Add the 0-0 axis to the fiber surface.");

    std::optional<double> perturbationEpsilon;
    cliApp.add_option("--epsilon, -e", perturbationEpsilon, "Strength of the numerial perturbation in the range [-e, e].");

    string outputSheetPolygonsFilename;
    cliApp.add_option("--outputSheetPolygons, -o", outputSheetPolygonsFilename, "Filename where to output the coordinates of the polygons that represent each sheet.");

    string saveReebSpaceFile;
    cliApp.add_option("--saveReebSpace, -s", saveReebSpaceFile, "Save the Reeb space to disk with this filename.");

    string readReebSpaceFile;
    cliApp.add_option("--readReebSpace, -l", readReebSpaceFile, "Load the Reeb space from disk with this filename.");

    string saveReebSpaceSheetsInfoFile;
    cliApp.add_option("--saveSheetInfo, -i", saveReebSpaceSheetsInfoFile, "Save info about the sheets.");

    string fiberBenchmarkFile;
    cliApp.add_option("--fiberPerformanceTimingsFile, -b", fiberBenchmarkFile, "Benchmakr for timings.");

    int sheetOutputCount = 10;
    cliApp.add_option("--sheetOutputCount", sheetOutputCount, "How many sheets to sample for automatic feature extraction.");

    string outputSheetFibersFolder;
    cliApp.add_option("--outputSheetFibersFolder", outputSheetFibersFolder, "Folder in which to ouput fiber for each sheet.");

    std::optional<float> fMin;
    cliApp.add_option("--fMin", fMin, "Set the min value for the f scalar field.");

    std::optional<float> gMin;
    cliApp.add_option("--gMin", gMin, "Set the min value for the g scalar field.");

    std::optional<float> fMax;
    cliApp.add_option("--fMax", fMax, "Set the max value for the f scalar field.");

    std::optional<float> gMax;
    cliApp.add_option("--gMax", gMax, "Set the max value for the g scalar field.");

    std::optional<float> fieldFValueFS;
    cliApp.add_option("--fieldFValueFS", fieldFValueFS, "Set value to compute an FS for the f field.");

    std::optional<float> fieldGValueFS;
    cliApp.add_option("--fieldGValueFS", fieldGValueFS, "Set value to compute an FS for the g field.");

    string fName = "";
    cliApp.add_option("--fName", fName, "The name of the f field to read from the input data.");

    string gName = "";
    cliApp.add_option("--gName", gName, "The name of the g field to read from the input data.");

    int sheetsToProcess = 20;
    cliApp.add_option("--sheetsToProcess", sheetsToProcess, "How many of the top sheets would you like to process?. Default is 20.");

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
        tetMesh = io::readData(filename, fName, gName);
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

    if (perturbationEpsilon.has_value())
    {
        Timer::start();
        tetMesh.perturbRangeValues(perturbationEpsilon.value(), fName, gName);
        Timer::stop("Perturbing range values                :");
    }

    Timer::start();
    tetMesh.sortVertices();
    Timer::stop("Sorting range points                   :");

    Timer::start();
    tetMesh.computeBoundingBoxes(fMin, fMax, gMin, gMax);
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

    Timer::start();
    tetMesh.computeSingularVertices();
    Timer::stop("Computing singular vertices            :");



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

    Timer::start();
    singularArrangement.assignHalfEdgePseudoSingular(tetMesh, singularArrangement);
    Timer::stop("Assigning pseudosingular edges         :");

    Timer::start();
    singularArrangement.computePointLocationDataStructure();
    Timer::stop("Arrangement search structure           :");





    //
    // Reeb space computation.
    //
    ReebSpace2 reebSpace2;

    if (readReebSpaceFile.empty())
    {
        Timer::start();
        reebSpace2.computeEdgeRegionSegments3(tetMesh, singularArrangement);
        Timer::stop("Computed red/blud intersetions         :");

        Timer::start();
        reebSpace2.determineEdgeRegionSegmentsOrientation(tetMesh, singularArrangement);
        Timer::stop("Edge regions plus/minus triangles      :");


        Timer::start();
        reebSpace2.computeVertexRegionSegments(tetMesh, singularArrangement);
        Timer::stop("Computed vertex regions                :");

        Timer::start();
        reebSpace2.determineVertexRegionSegmentsOrientation(tetMesh, singularArrangement);
        Timer::stop("Vertex regions plus/minus triangles    :");

        Timer::start();
        reebSpace2.determineEdgeCrossingSegmentsOriantation(tetMesh, singularArrangement);
        Timer::stop("Edge crossing plus/minus triangles     :");

        Timer::start();
        reebSpace2.traverse(tetMesh, singularArrangement, unitTestFiberGraphs);
        Timer::stop("Computed singular traversal            :");

        if (false == saveReebSpaceFile.empty())
        {
            io::saveReebSpace(reebSpace2, saveReebSpaceFile);
        }
    }
    else
    {
        Timer::start();
        reebSpace2 = io::loadReebSpace(readReebSpaceFile);
        Timer::stop("Read reeb space                        :");
    }








    Timer::start();
    reebSpace2.computeSheetBoundaries(singularArrangement);
    Timer::stop("Computing sheet boundaries             :");


    //for (const auto &[sheetId, boundary] : reebSpace2.sheetBoundaries)
    //{
        //printf("--------------------------------------------------- sheet with ID %d has this boundary: \n", sheetId);

        //for (const Halfedge_const_handle &he : boundary)
        //{
            //std::cout << he->source()->point() << std::endl;

        //}


    //}


    Timer::start();
    reebSpace2.computeSheets(singularArrangement);
    Timer::stop("Postprocessing                         :");


    // Another part of postprocessing to make sure we can do openmp for the fiber labeling of 
    for (auto he = singularArrangement.arr.halfedges_begin(); he != singularArrangement.arr.halfedges_end(); ++he)
    {
        auto& curve = he->curve();   // gets the Arr_segment_2 or whatever curve type
        curve.is_vertical();          // forces lazy _is_vertical initialization
    }


    Timer::start();
    singularArrangement.buildAABBtree(tetMesh);
    Timer::stop("Arrangement AABB Tree                  :");


    std::cout << "\nThe Reeb space has this many sheets " << reebSpace2.correspondenceGraphDS.countComponents() << std::endl;

    int correspondenceGraphSize = 0;
    for (const auto &correspondenceGraph : reebSpace2.correspondenceGraph)
    {
        correspondenceGraphSize += correspondenceGraph.size();
    }


    // This is the old computation, keep these empty unless we want to unit test
    //
    ReebSpace reebSpace;
    Arrangement arrangement;

    if (unitTestSheets || unitTestFiberGraphs)
    {
        Timer::start();
        arrangement.computeArrangement(tetMesh, Arrangement::SegmentMode::UseAllSegments);
        Timer::stop("Arrangement                            :");

        Timer::start();
        arrangement.computePointLocationDataStructure();
        Timer::stop("Arrangement search structure           :");


        Timer::start();
        reebSpace.computeTraversal(tetMesh, arrangement, discardFiberSeeds, unitTestFiberGraphs);
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


        if (true == unitTestFiberGraphs)
        {
            Timer::start();
            bool arePreimageGraphsEqual = reebSpace2.unitTestCompareFiberGraphs(tetMesh, singularArrangement, arrangement, reebSpace);
            Timer::stop("Comparing preimage graphs              :");

            if (false == arePreimageGraphsEqual)
            {
                std::cerr << "----------------------------------------------------------------------------------------------------------------\n";
                std::cerr << "--------------------------------- PREIMAGE GRAPHS NOT EQUAL!!!--------------------------------------------------\n";
                std::cerr << "----------------------------------------------------------------------------------------------------------------\n";
                return 1;
            }
        }


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
            //io::saveSheets(tetMesh, arrangement, reebSpace, outputSheetPolygonsFilename + ".old.vtp");
            io::saveSheets2(tetMesh, singularArrangement, reebSpace2, outputSheetPolygonsFilename);
            io::saveSheetsFeatures(tetMesh, singularArrangement, reebSpace2, outputSheetPolygonsFilename + ".features.vtp");
            io::saveSheetGraph(reebSpace2, outputSheetPolygonsFilename + ".graph.dot");
        }
        catch (const std::exception &e)
        {
            std::cerr << "Error: " << e.what() << '\n';
            return 1;
        }
    }

    if (fieldFValueFS.has_value())
    {
        for (int i = 0 ; i < sheetsToProcess ; i++)
        {
            const int sheetId = reebSpace2.orderSheet[i];

            const std::vector<std::array<double, 2>> controlPoints{
                    {fieldFValueFS.value(), -1e10}, 
                    {fieldFValueFS.value(), +1e10}
            };

            FiberSurface fs = FiberSurface::constructSegmentedFiberSurface(tetMesh, singularArrangement, reebSpace2, controlPoints, {sheetId});

            std::string fsFilename = "./output/fs.f." + std::to_string(sheetId) + ".vtp";
            io::saveFiberSurface({fs}, fsFilename);

            std::cout << "Saved f-field labeled FS in " << fsFilename << std::endl;
        }
    }

    if (fieldGValueFS.has_value())
    {
        for (int i = 0 ; i < sheetsToProcess ; i++)
        {
            const int sheetId = reebSpace2.orderSheet[i];

            const std::vector<std::array<double, 2>> controlPoints{
                    {-1e10, fieldGValueFS.value()}, 
                    {+1e10, fieldGValueFS.value()}
            };

            FiberSurface fs = FiberSurface::constructSegmentedFiberSurface(tetMesh, singularArrangement, reebSpace2, controlPoints, {sheetId});

            std::string fsFilename = "./output/fs.g." + std::to_string(sheetId) + ".vtp";
            io::saveFiberSurface({fs}, fsFilename);

            std::cout << "Saved g-field labeled FS in " << fsFilename << std::endl;
        }
    }


    if (false == fiberBenchmarkFile.empty())
    {
        //performance::testInteractiveFiberPerformance(tetMesh, singularArrangement, reebSpace2, 1000, fiberBenchmarkFile);
        performance::testInteractiveFiberSurfacePerformance(tetMesh, singularArrangement, reebSpace2, 100, fiberBenchmarkFile);
        return 0;
    }

    //io::saveOriginalMesh("og.vtu", tetMesh.originalMesh);
    //io::readDataVtp("/home/peter/Projects/data/reeb-space-test-data/nana/trajectories/State_2/fiberSurfaceExample.vtp");


    if (headless)
    {
        return 0;
    }

    // Set up QT Application
    QApplication app(argc, argv);
    glutInit(&argc, argv);

    // Package all my data for visualisation
    Data data(tetMesh, arrangement, singularArrangement, reebSpace, reebSpace2);
    data.zeroAxis = addZeroAxis;


    // New Stuff
    if (false == saveReebSpaceSheetsInfoFile.empty())
    {
        data.sheetRegularVertices = augmentation::computeRegularVertexSheets(tetMesh, singularArrangement, reebSpace2); 

        std::cerr << "\nSaved regular vertex sheets to file " << saveReebSpaceSheetsInfoFile;
        io::writeSheetData(saveReebSpaceSheetsInfoFile, tetMesh.isVertexSingular, reebSpace2.sheetArea, data.sheetRegularVertices);
        //return 0;
    }

    if (false == moleculeFilename.empty())
    {
        data.molecule = io::readMolecule(moleculeFilename);
    }

    // Create the widget
    TracerVisualiserWindow* window = new TracerVisualiserWindow(NULL, data);
    window->setWindowTitle("RS Explorer");

    // Make the window full screen by default
    //window->showMaximized();

    window->setWindowState(Qt::WindowNoState);
    //window->setMinimumSize(1800, 1200);
    window->setMinimumSize(1200, 800);
    window->showNormal();
    window->move(0, 0);

    // Show the label
    window->show();

    // start it running
    app.exec();

    // clean up
    delete window;

    // return to caller
    return 0;
} // main()
