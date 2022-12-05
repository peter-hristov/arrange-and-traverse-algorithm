#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include <QApplication>
#include "./TracerVisualiserWindow.h"

using namespace std;


int main(int argc, char* argv[])
{
    if (argc != 2)
    {
        cout << "The usage is ./tv9k <input_file>." << endl;
        return 0;
    }

    std::string filename = argv[1];

    // Set up the data class
    Data* data = new Data();

    // Read in data file
    //data->readData(filename);
    data->readDataGrid(filename);

    return 0;

    // Compute domain/range data bounds
    data->computeMinMaxRangeDomainCoordinates();

    // Set up QT Application
    QApplication app(argc, argv);
    glutInit(&argc, argv);

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
    delete data;

    // return to caller
    return 0;
} // main()
