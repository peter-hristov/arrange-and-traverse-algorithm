#include <cstddef>
#include <filesystem>

#include <GL/glut.h>
#include <QApplication>

#include "./io.h"
#include "./Timer.h"
#include "./ReebSpace.h"
#include "./Data.h"
#include "./CGALTypedefs.h"
#include "./Arrangement.h"
#include "./utility/CLI11.hpp"
#include "./TracerVisualiserWindow.h"

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
        tetMesh = io::readData(filename);
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error: " << e.what() << '\n';
        return 1;
    }
    Timer::start();
    tetMesh.perturbRangeValues(perturbationEpsilon);
    tetMesh.sortVertices();
    tetMesh.computeBoundingBoxes();
    tetMesh.computeCombinatorialStructure();
    tetMesh.computeUpperLowerLinkAndStar();
    tetMesh.computeSingularEdgeTypes();
    Timer::stop("Input mesh postprocessing              :");

    Timer::start();
    Arrangement arrangement;
    arrangement.computeArrangement(tetMesh);
    Timer::stop("Arrangement                            :");


    for (auto vit = arrangement.arr.vertices_begin(); vit != arrangement.arr.vertices_end(); ++vit)
    {
        const auto& v = *vit;

        if (false == arrangement.arrangementPointIndices.contains(vit->point()))
        {
            assert(v.degree() == 4);
            //std::cout << "Vertex at " << v.point() << " has degree " << v.degree() << "\n";
        } 
    }

    for (auto heit = arrangement.arr.halfedges_begin(); heit != arrangement.arr.halfedges_end(); ++heit)
    {
        // Get iterators to originating curves via the arrangement object
        auto oc_begin = arrangement.arr.originating_curves_begin(heit);
        auto oc_end   = arrangement.arr.originating_curves_end(heit);

        std::size_t count = std::distance(oc_begin, oc_end);

        assert(count == 1);

        //std::cout << "Halfedge from "
            //<< heit->source()->point() << " to "
            //<< heit->target()->point()
            //<< " has " << count << " originating curve(s)\n";
    }



    //cout << tetMesh.vertexCoordinatesF.size() << " " << tetMesh.edges.size() << "\n";

    //// Print vertices
    //for (int i = 0 ; i < tetMesh.vertexCoordinatesF.size() ; i++)
    //{
        //cout << tetMesh.vertexCoordinatesF[i] << " " << tetMesh.vertexCoordinatesG[i] << "\n";
    //}

    ////// Print edges
    //for (const auto &edge : tetMesh.edges)
    //{
        //cout << edge[0] << " " << edge[1] << "\n";

        ////cout << tetMesh.vertexCoordinatesF[edge[0]] << " " << tetMesh.vertexCoordinatesG[edge[0]] << " " << tetMesh.vertexCoordinatesF[edge[1]] << " " << tetMesh.vertexCoordinatesG[edge[1]] << "\n";
    //}



    return 0;

    Timer::start();
    arrangement.computePointLocationDataStructure();
    Timer::stop("Arrangement search structure           :");

    Timer::start();
    ReebSpace reebSpace;
    reebSpace.computeTraversal(tetMesh, arrangement, discardFiberSeeds);
    Timer::stop("Computed {G_F} and H                   :");

    std::cout << "Postprocessing..." << std::endl;
    Timer::start();
    reebSpace.computeSheetGeometry(tetMesh, arrangement);
    reebSpace.computeSheetArea(tetMesh, arrangement);
    reebSpace.printTopSheets(tetMesh, arrangement, 20);
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
    Data data(tetMesh, arrangement, reebSpace);

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
