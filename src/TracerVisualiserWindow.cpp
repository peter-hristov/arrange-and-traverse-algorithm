#include <iostream>
#include <qnamespace.h>
#include <string>
#include <utility>
#include <QCheckBox>
#include <QHBoxLayout>


#include "./TracerVisualiserWindow.h"
#include "./Data.h"
#include "./io.h"
#include "src/ReebSpace2.h"

using namespace std;
void
TracerVisualiserWindow::keyPressEvent(QKeyEvent* event)
{
    const int moveSpeed = 2;

    if (event->key() == Qt::Key_U) {
        QApplication::quit();
    }

    if (event->key() == Qt::Key_W &&
            event->modifiers() & Qt::ControlModifier)
    {
        close();   // closes this window
        return;
    }

    if (event->key() == Qt::Key_Backspace) {
        if (this->plotWidget->controlPoints.size() > 0)
        {
            this->plotWidget->controlPoints.pop_back();
        }
        this->plotWidget->recomputeFiber = true;
        this->plotWidget->update();
        this->update();
    }

    if (event->key() == Qt::Key_Return || event->key() == Qt::Key_Enter) {
        this->plotWidget->recomputeFiberSurface = true;
        this->plotWidget->update();
        this->update();
    }

    if (event->key() == Qt::Key_I) {
        this->plotWidget->mousePoint.setY(this->plotWidget->mousePoint.y() - moveSpeed);
        this->plotWidget->recomputeFiber = true;
        this->update();
    }
    if (event->key() == Qt::Key_J) {
        this->plotWidget->mousePoint.setX(this->plotWidget->mousePoint.x() - moveSpeed);
        this->plotWidget->recomputeFiber = true;
        this->update();
    }
    if (event->key() == Qt::Key_K) {
        this->plotWidget->mousePoint.setY(this->plotWidget->mousePoint.y() + moveSpeed);
        this->plotWidget->recomputeFiber = true;
        this->update();
    }
    if (event->key() == Qt::Key_L) {
        this->plotWidget->mousePoint.setX(this->plotWidget->mousePoint.x() + moveSpeed);
        this->plotWidget->recomputeFiber = true;
        this->update();
    }


    if (event->key() == Qt::Key_7) {
        this->tracerVisualiserWidget->fiberColour = 0;
        this->tracerVisualiserWidget->update();
    }
    if (event->key() == Qt::Key_8) {
        this->tracerVisualiserWidget->fiberColour = 1;
        this->tracerVisualiserWidget->update();
    }
    if (event->key() == Qt::Key_9) {
        this->tracerVisualiserWidget->fiberColour = 2;
        this->tracerVisualiserWidget->update();
    }

    if (event->key() == Qt::Key_C) {
        checkboxShowTraces->setChecked(!checkboxShowTraces->isChecked());
    }

    if (event->key() == Qt::Key_D) {
        this->computeTracedFiberSurfaceButton->click();
    }

    if (event->key() == Qt::Key_H) {
        //this->data.printSheetHistogram();
    }
    if (event->key() == Qt::Key_S) {
        // Save the fibers
        std::string filename = "./output/fibers.vtp";
        std::cout << "Saving fibers to " << filename << std::endl;
        io::saveFibers(this->tracerVisualiserWidget->faceFibers, filename);

        filename = "./output/fiber-surface.vtp";
        io::saveFiberSurface(this->data.surfaceMeshes, filename);
        std::cout << "Saving surfaces to to " << filename << std::endl;

        filename = "./output/reeb-space.png";
        this->plotWidget->saveToFile(filename);
        std::cout << "Saving surfaces to to " << filename << std::endl;
    }

}

TracerVisualiserWindow::TracerVisualiserWindow(QWidget* parent, Data &_data)
    : QWidget(parent)
      , data(_data)    // <-- initialize reference here
{

    // Initialise Widgets
    plotWidget = new PlotWidget(this, data);
    tracerVisualiserWidget = new TracerVisualiserWidget(this, data);

    plotWidget->sibling = tracerVisualiserWidget;
    tracerVisualiserWidget->sibling = plotWidget;

    checkboxShowFibers = new QCheckBox("Show Fibers");
    checkboxShowFibers->setChecked(true);

    checkboxShowFiberSurfaces = new QCheckBox("Show Fiber Surface");
    checkboxShowFiberSurfaces->setChecked(true);

    //checkboxShowFaces = new QCheckBox("Show Faces");
    //checkboxShowFaces->setChecked(true);


    vertexOpacitySlider = new QSlider(Qt::Horizontal);
    vertexOpacitySlider->setValue(this->tracerVisualiserWidget->vertexOpacity * 100);

    edgeOpacitySlider = new QSlider(Qt::Horizontal);
    edgeOpacitySlider->setValue(this->tracerVisualiserWidget->edgeOpacity * 100);

    faceOpacitySlider = new QSlider(Qt::Horizontal);
    faceOpacitySlider->setValue(this->tracerVisualiserWidget->faceOpacity * 100);

    fakeSlider = new QSlider(Qt::Horizontal);
    fakeSlider->setTracking(false);

    checkboxShowTraces = new QCheckBox("Show fiber point trace.");
    this->computeTracedFiberSurfaceButton = new QPushButton("Compute Traced Fiber Surface", this);
    this->clearAllButton = new QPushButton("Clear All", this);

    this->clearFibersButton = new QPushButton("Clear Fibers", this);
    this->clearFiberSurfaceButton = new QPushButton("Clear FS", this);
    this->clearSelectedSheetsButton = new QPushButton("Clear Sheets", this);


    // Create widgets
    this->spinBoxAddSheet = new QSpinBox(this);
    this->buttonAddSheet = new QPushButton("Add sheet", this);

    spinBoxAddSheet->setRange(1, FiberGraph::componentCount);


    //
    // Layouts
    //
    optionsLayout = new QGridLayout();
    optionsLayout2 = new QGridLayout();

    auto rowOneLayout = new QGridLayout();
    rowOneLayout->addWidget(checkboxShowFibers,0, 0);
    //rowOneLayout->addWidget(vertexOpacitySlider, 0, 1);
    rowOneLayout->addWidget(checkboxShowFiberSurfaces,1, 0);
    //rowOneLayout->addWidget(edgeOpacitySlider, 1, 1);
    //rowOneLayout->addWidget(checkboxShowFaces,2, 0);
    //rowOneLayout->addWidget(faceOpacitySlider, 2, 1);

    optionsLayout->addLayout(rowOneLayout, 0, 0);


    optionsLayout2->addWidget(clearFibersButton, 0, 0);
    optionsLayout2->addWidget(clearFiberSurfaceButton, 0, 1);
    optionsLayout2->addWidget(clearSelectedSheetsButton, 0, 2);
    optionsLayout2->addWidget(clearAllButton, 0, 3);



    optionsLayout2->addWidget(checkboxShowTraces, 1, 0);
    optionsLayout2->addWidget(computeTracedFiberSurfaceButton, 1, 1);
    //optionsLayout2->addWidget(fakeSlider, 1, 1);

    optionsLayout2->addWidget(spinBoxAddSheet, 2, 0);
    optionsLayout2->addWidget(buttonAddSheet, 2, 1);


    // Set up layout
    windowLayout = new QGridLayout(this);
    windowLayout->addWidget(tracerVisualiserWidget, 0, 0);
    windowLayout->addWidget(plotWidget, 0, 1);

    windowLayout->addLayout(optionsLayout, 1, 0);
    windowLayout->addLayout(optionsLayout2, 1, 1);

    connect(spinBoxAddSheet, &QSpinBox::editingFinished, buttonAddSheet, &QPushButton::click);

    connect(buttonAddSheet, &QPushButton::clicked, this, [this]() {
            const int sheetId = spinBoxAddSheet->value();

            if (this->data.reebSpace2.sheetArea.contains(sheetId))
            {
                this->tracerVisualiserWidget->selectedSheetIds.insert(sheetId);

                this->plotWidget->staticReebSpaceCache = nullptr;
                this->plotWidget->update();
                this->tracerVisualiserWidget->update();
            }
            });


    connect(this->clearFibersButton, &QPushButton::clicked, this, [this]() {
            this->plotWidget->fiberPointsTraces.clear();
            this->plotWidget->fiberPointsTraces.shrink_to_fit();

            this->tracerVisualiserWidget->clearFiber();

            this->plotWidget->update();
            this->tracerVisualiserWidget->update();
            });


    connect(this->clearFiberSurfaceButton, &QPushButton::clicked, this, [this]() {
            this->plotWidget->controlPoints.clear();
            this->plotWidget->controlPoints.shrink_to_fit();
            this->tracerVisualiserWidget->clearFiberSurface();

            this->plotWidget->update();
            this->tracerVisualiserWidget->update();
            });
    connect(this->clearSelectedSheetsButton, &QPushButton::clicked, this, [this]() {
            this->plotWidget->staticReebSpaceCache = nullptr;
            this->tracerVisualiserWidget->selectedSheetIds = {};

            this->plotWidget->update();
            this->tracerVisualiserWidget->update();
            });

    connect(this->clearAllButton, &QPushButton::clicked, this, [this]() {
            this->plotWidget->fiberPointsTraces.clear();
            this->plotWidget->fiberPointsTraces.shrink_to_fit();
            this->plotWidget->controlPoints.clear();
            this->plotWidget->controlPoints.shrink_to_fit();
            this->plotWidget->staticReebSpaceCache = nullptr;

            this->tracerVisualiserWidget->selectedSheetIds = {};
            this->tracerVisualiserWidget->clearFiber();
            this->tracerVisualiserWidget->clearFiberSurface();

            this->plotWidget->update();
            this->tracerVisualiserWidget->update();
            });

    connect(this->computeTracedFiberSurfaceButton, &QPushButton::clicked, this, [this]() {
            if (this->plotWidget->fiberPointsTraces.size() == 1)
            {
            this->plotWidget->controlPoints = std::move(this->plotWidget->fiberPointsTraces[0]);

            this->plotWidget->recomputeFiberSurface = true;
            this->plotWidget->update();
            this->tracerVisualiserWidget->update();
            }
            });

    connect(checkboxShowFibers, &QCheckBox::toggled, [=](bool checked) {
            this->tracerVisualiserWidget->drawFibers = checked;
            this->tracerVisualiserWidget->generateDisplayList();
            this->tracerVisualiserWidget->update();
            });

    connect(checkboxShowFiberSurfaces, &QCheckBox::toggled, [=](bool checked) {
            this->tracerVisualiserWidget->drawFiberSurfaces = checked;
            this->tracerVisualiserWidget->generateDisplayList();
            this->tracerVisualiserWidget->update();
            });

    //connect(checkboxShowFaces, &QCheckBox::toggled, [=](bool checked) {
            //this->tracerVisualiserWidget->drawFaces = checked;

            //this->tracerVisualiserWidget->update();
            //});

    connect(checkboxShowTraces, &QCheckBox::toggled, [=](bool checked) {

            this->plotWidget->fiberPointsTraces.clear();
            this->plotWidget->fiberPointsTraces.shrink_to_fit();
            this->plotWidget->update();

            this->tracerVisualiserWidget->clearFibers = !this->tracerVisualiserWidget->clearFibers;
            this->tracerVisualiserWidget->updateFiber({});
            this->tracerVisualiserWidget->update();
            });

    connect(this->vertexOpacitySlider, &QSlider::valueChanged, plotWidget, [=]() {
        this->tracerVisualiserWidget->vertexOpacity = static_cast<double>(this->vertexOpacitySlider->value()) / 100.0;
        this->tracerVisualiserWidget->update();
    });

    connect(this->faceOpacitySlider, &QSlider::valueChanged, plotWidget, [=]() {
        this->tracerVisualiserWidget->faceOpacity = static_cast<double>(this->faceOpacitySlider->value()) / 100.0;
        this->tracerVisualiserWidget->update();
    });

    connect(this->edgeOpacitySlider, &QSlider::valueChanged, plotWidget, [=]() {
        this->tracerVisualiserWidget->edgeOpacity = static_cast<double>(this->edgeOpacitySlider->value()) / 100.0;
        this->tracerVisualiserWidget->update();
    });
}

TracerVisualiserWindow::~TracerVisualiserWindow()
{
    delete plotWidget;
    delete windowLayout;
    delete tracerVisualiserWidget;
}
