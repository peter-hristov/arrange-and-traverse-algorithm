#include <iostream>
#include <qnamespace.h>
#include <string>
#include <utility>
#include <QCheckBox>
#include <QHBoxLayout>
#include <QGroupBox>
#include <QPlainTextEdit>



#include "./TracerVisualiserWindow.h"
#include "./Data.h"
#include "./io.h"
#include "./ReebSpace2.h"

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
        this->plotWidget->update();
        this->update();
    }

    if (event->key() == Qt::Key_Return || event->key() == Qt::Key_Enter) {
        this->computeFiberSurfaceButton->click();
    }

    if (event->key() == Qt::Key_I) {
        this->plotWidget->mousePoint.setY(this->plotWidget->mousePoint.y() - moveSpeed);

        if (this->buttonShowTraces->isChecked()) 
        {
            if (this->plotWidget->fiberPointsTraces.size() == 0)
            {
                this->plotWidget->fiberPointsTraces.push_back({});
            }
            this->plotWidget->fiberPointsTraces.back().push_back(this->plotWidget->mousePoint);
        }

        this->plotWidget->recomputeFiber = true;
        this->update();
    }
    if (event->key() == Qt::Key_J) {
        this->plotWidget->mousePoint.setX(this->plotWidget->mousePoint.x() - moveSpeed);
        if (this->buttonShowTraces->isChecked()) 
        {
            if (this->plotWidget->fiberPointsTraces.size() == 0)
            {
                this->plotWidget->fiberPointsTraces.push_back({});
            }
            this->plotWidget->fiberPointsTraces.back().push_back(this->plotWidget->mousePoint);
        }
        this->plotWidget->recomputeFiber = true;
        this->update();
    }
    if (event->key() == Qt::Key_K) {
        this->plotWidget->mousePoint.setY(this->plotWidget->mousePoint.y() + moveSpeed);

        if (this->buttonShowTraces->isChecked()) 
        {
            if (this->plotWidget->fiberPointsTraces.size() == 0)
            {
                this->plotWidget->fiberPointsTraces.push_back({});
            }
            this->plotWidget->fiberPointsTraces.back().push_back(this->plotWidget->mousePoint);
        }
        this->plotWidget->recomputeFiber = true;
        this->update();
    }
    if (event->key() == Qt::Key_L) {
        this->plotWidget->mousePoint.setX(this->plotWidget->mousePoint.x() + moveSpeed);

        if (this->buttonShowTraces->isChecked()) 
        {
            if (this->plotWidget->fiberPointsTraces.size() == 0)
            {
                this->plotWidget->fiberPointsTraces.push_back({});
            }
            this->plotWidget->fiberPointsTraces.back().push_back(this->plotWidget->mousePoint);
        }

        this->plotWidget->recomputeFiber = true;
        this->update();
    }

    if (event->key() == Qt::Key_C) {
        this->buttonShowTraces->click();
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
        io::saveFibers(this->data.fibers, this->data.reebSpace2, filename);

        filename = "./output/fiberTraces.vtp";
        std::cout << "Saving fiber to " << filename << std::endl;
        io::saveFiberTraces(this->plotWidget, filename);

        filename = "./output/fiber-surface.vtp";
        std::cout << "Saving surfaces to to " << filename << std::endl;
        io::saveFiberSurface(this->data.fiberSurfaces, filename);
        //io::saveFiberPointsAsTriangleSoup(this->tracerVisualiserWidget->faceFiberSurface, filename);

        filename = "./output/fiber-surface-features.vtp";
        std::cout << "Saving surfaces to to " << filename << std::endl;
        io::saveFiberSurface(this->data.featureSurfaces, filename);
        //io::saveFiberPointsAsTriangleSoup(this->tracerVisualiserWidget->faceFiberSurfaceFeatures, filename);

        filename = "./output/reeb-space.png";
        std::cout << "Saving surfaces to to " << filename << std::endl;

        std::cerr << "This many fiber " << this->data.fibers.size() << " and this many trace points " << plotWidget->fiberPointsTraces.back().size();

        this->plotWidget->saveToFile(filename);
        this->plotWidget->shouldSaveSheets = true;

        this->plotWidget->update();
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

    checkboxShowFeatures = new QCheckBox("Show Features");
    checkboxShowFeatures->setChecked(true);

    //checkboxShowFaces = new QCheckBox("Show Faces");
    //checkboxShowFaces->setChecked(true);


    fiberOpacitySlider = new QSlider(Qt::Horizontal);
    fiberOpacitySlider->setMinimum(0);
    fiberOpacitySlider->setMaximum(100);
    fiberOpacitySlider->setValue(100);

    fiberSurfaceOpacitySlider = new QSlider(Qt::Horizontal);
    fiberSurfaceOpacitySlider->setMinimum(0);
    fiberSurfaceOpacitySlider->setMaximum(100);
    fiberSurfaceOpacitySlider->setValue(100);

    featureSurfaceOpacitySlider = new QSlider(Qt::Horizontal);
    featureSurfaceOpacitySlider->setMinimum(0);
    featureSurfaceOpacitySlider->setMaximum(100);
    featureSurfaceOpacitySlider->setValue(100);

    //fakeSlider = new QSlider(Qt::Horizontal);
    //fakeSlider->setTracking(false);

    //this->checkboxShowTraces = new QCheckBox("Trace Fiber.");


    this->buttonShowTraces = new QPushButton("Trace Fiber", this);
    buttonShowTraces->setCheckable(true);
    buttonShowTraces->setChecked(false);

    this->computeTracedFiberSurfaceButton = new QPushButton("Trace FS", this);
    this->computeFiberSurfaceFeatureButton = new QPushButton("Features", this);
    this->computeFiberSurfaceButton = new QPushButton("FS", this);

    this->clearAllButton = new QPushButton("All", this);

    this->clearFibersButton = new QPushButton("Fibers", this);
    this->clearFiberSurfaceButton = new QPushButton("FS", this);
    this->clearSelectedSheetsButton = new QPushButton("Sheets", this);

    //this->buttonAddNewControlPolygon = new QPushButton("Add FSCP", this);


    // Create widgets
    this->spinBoxAddSheet = new QSpinBox(this);
    spinBoxAddSheet->setRange(1, FiberGraph::componentCount);

    this->buttonAddSheet = new QPushButton("Add", this);

    this->spinBoxAddTopSheets = new QSpinBox(this);
    spinBoxAddTopSheets->setRange(1, 100);
    this->buttonAddTopSheets = new QPushButton("Add", this);

    //
    // Layouts
    //
    optionsLayout = new QGridLayout();
    optionsLayout2 = new QGridLayout();

    QGroupBox* visibilityGroup = new QGroupBox("Visibility");
    auto rowOneLayout = new QGridLayout(visibilityGroup);
    rowOneLayout->setContentsMargins(4, 4, 4, 4);
    rowOneLayout->setSpacing(4);
    rowOneLayout->addWidget(checkboxShowFibers,        0, 0);
    rowOneLayout->addWidget(fiberOpacitySlider,       0, 1);
    rowOneLayout->addWidget(checkboxShowFiberSurfaces, 1, 0);
    rowOneLayout->addWidget(fiberSurfaceOpacitySlider,         1, 1);
    rowOneLayout->addWidget(checkboxShowFeatures,      2, 0);
    rowOneLayout->addWidget(featureSurfaceOpacitySlider,                2, 1);
    optionsLayout->addWidget(visibilityGroup, 0, 1);



    // --- Clear group box (2x2) ---
    QGroupBox* clearGroup = new QGroupBox("Clear");
    QGridLayout* clearGrid = new QGridLayout(clearGroup);
    clearGrid->addWidget(clearFibersButton,        0, 0);
    clearGrid->addWidget(clearFiberSurfaceButton,  0, 1);
    clearGrid->addWidget(clearSelectedSheetsButton,1, 0);
    clearGrid->addWidget(clearAllButton,           1, 1);

    optionsLayout2->addWidget(clearGroup, 0, 0);








    // --- Compute group box (2x2) ---
    QGroupBox* computeGroup = new QGroupBox("Compute");
    QGridLayout* computeGrid = new QGridLayout(computeGroup);
    computeGrid->addWidget(computeFiberSurfaceButton,        0, 0);
    computeGrid->addWidget(computeFiberSurfaceFeatureButton, 0, 1);
    computeGrid->addWidget(computeTracedFiberSurfaceButton,  1, 0);
    //computeGrid->addWidget(checkboxShowTraces,               1, 1);
    computeGrid->addWidget(buttonShowTraces,               1, 1);


    optionsLayout2->addWidget(computeGroup, 0, 1);

    //optionsLayout2->addWidget(clearFibersButton, 0, 0);
    //optionsLayout2->addWidget(clearFiberSurfaceButton, 0, 1);
    //optionsLayout2->addWidget(clearSelectedSheetsButton, 0, 2);
    //optionsLayout2->addWidget(clearAllButton, 0, 3);



    //optionsLayout2->addWidget(checkboxShowTraces, 1, 0);
    //optionsLayout2->addWidget(computeFiberSurfaceButton, 1, 1);
    //optionsLayout2->addWidget(computeFiberSurfaceFeatureButton, 1, 2);
    //optionsLayout2->addWidget(computeTracedFiberSurfaceButton, 1, 3);
    ////optionsLayout2->addWidget(buttonAddNewControlPolygon, 1, 4);

    ////optionsLayout2->addWidget(fakeSlider, 1, 1);




    //optionsLayout2->addWidget(spinBoxAddSheet, 1, 0);
    //optionsLayout2->addWidget(buttonAddSheet, 1, 1);
    //optionsLayout2->addWidget(spinBoxAddTopSheets, 1, 2);
    //optionsLayout2->addWidget(buttonAddTopSheets, 1, 3);

    QGroupBox* addSheetsGroup = new QGroupBox("Select Sheets");
    QGridLayout* addSheetsGrid = new QGridLayout(addSheetsGroup);
    addSheetsGrid->setContentsMargins(4, 4, 4, 4);
    addSheetsGrid->setSpacing(4);

    addSheetsGrid->addWidget(new QLabel("By ID:"), 0, 0);
    addSheetsGrid->addWidget(spinBoxAddSheet,       0, 1);
    addSheetsGrid->addWidget(buttonAddSheet,        0, 2);

    addSheetsGrid->addWidget(new QLabel("Top N:"), 1, 0);
    addSheetsGrid->addWidget(spinBoxAddTopSheets,   1, 1);
    addSheetsGrid->addWidget(buttonAddTopSheets,    1, 2);

    optionsLayout2->addWidget(addSheetsGroup, 0, 2);


    QGroupBox* infoGroup = new QGroupBox("Information");
    QVBoxLayout* infoLayout = new QVBoxLayout(infoGroup);
    infoLayout->setContentsMargins(4, 4, 4, 4);
    infoLayout->setSpacing(2);

    infoBox = new QPlainTextEdit();
    infoBox->setTextInteractionFlags(Qt::TextSelectableByMouse | Qt::TextSelectableByKeyboard);
    infoLayout->addWidget(infoBox);

    infoGroup->setFixedHeight(visibilityGroup->sizeHint().height());
    optionsLayout->addWidget(infoGroup, 0, 0);



    // Set up layout
    windowLayout = new QGridLayout(this);
    QGroupBox* domainGroup = new QGroupBox("Domain View");
    QGroupBox* rangeGroup = new QGroupBox("Range View");

    QVBoxLayout* domainLayout = new QVBoxLayout(domainGroup);
    QVBoxLayout* rangeLayout = new QVBoxLayout(rangeGroup);

    domainLayout->addWidget(tracerVisualiserWidget);
    rangeLayout->addWidget(plotWidget);

    windowLayout->addWidget(domainGroup, 0, 0);
    windowLayout->addWidget(rangeGroup, 0, 1);

    windowLayout->addLayout(optionsLayout, 1, 0);
    windowLayout->addLayout(optionsLayout2, 1, 1);

    // After setting up windowLayout, fix column widths
    windowLayout->setColumnStretch(0, 1);  // left view gets 1 part
    windowLayout->setColumnStretch(1, 1);  // right view gets 1 part

    windowLayout->setRowStretch(0, 1);  // views row takes all extra space
    windowLayout->setRowStretch(1, 0);  // options row stays minimum height










    //connect(spinBoxAddSheet, &QSpinBox::editingFinished, buttonAddSheet, &QPushButton::click);
    //connect(spinBoxAddSheet, &QSpinBox::editingFinished, buttonAddSheet, &QPushButton::click);


    connect(buttonAddSheet, &QPushButton::clicked, this, [this]() {
            const int sheetId = spinBoxAddSheet->value();

            if (this->data.reebSpace2.sheetArea.contains(sheetId))
            {
                this->data.selectedSheetIds.insert(sheetId);

                this->plotWidget->staticReebSpaceCache = nullptr;
                this->plotWidget->update();
                this->tracerVisualiserWidget->update();
            }
            this->updateSelectedSheets(this->data.selectedSheetIds);
            });




    connect(buttonAddTopSheets, &QPushButton::clicked, this, [this]() {
            const size_t numberOfSheets = spinBoxAddTopSheets->value();
            this->data.selectedSheetIds = {};

            for (int i = 0 ; i < std::min(numberOfSheets, data.reebSpace2.orderSheet.size()) ; i++)
            {
                this->data.selectedSheetIds.insert(data.reebSpace2.orderSheet[i]);
            }


            this->plotWidget->staticReebSpaceCache = nullptr;
            this->plotWidget->update();
            this->tracerVisualiserWidget->update();

            this->updateSelectedSheets(this->data.selectedSheetIds);
            });



    connect(this->clearFibersButton, &QPushButton::clicked, this, [this]() {
            this->plotWidget->fiberPointsTraces.clear();
            this->plotWidget->fiberPointsTraces.shrink_to_fit();

            this->tracerVisualiserWidget->clearFibers();

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
            this->plotWidget->featureControlPolygons = {};
            this->data.selectedSheetIds = {};
            this->tracerVisualiserWidget->clearFiberSurfaceSheets();

            this->plotWidget->update();
            this->tracerVisualiserWidget->update();
            this->updateSelectedSheets({});
            });

    connect(this->clearAllButton, &QPushButton::clicked, this, [this]() {

            this->clearFibersButton->click();
            this->clearFiberSurfaceButton->click();
            this->clearSelectedSheetsButton->click();
            this->updateSelectedSheets({});

            //this->plotWidget->fiberPointsTraces.clear();
            //this->plotWidget->fiberPointsTraces.shrink_to_fit();
            //this->plotWidget->controlPoints.clear();
            //this->plotWidget->controlPoints.shrink_to_fit();
            //this->plotWidget->featureControlPolygons = {};
            //this->plotWidget->staticReebSpaceCache = nullptr;

            //this->tracerVisualiserWidget->selectedSheetIds = {};
            //this->tracerVisualiserWidget->clearFiber();
            //this->tracerVisualiserWidget->clearFiberSurface();

            //this->plotWidget->update();
            //this->tracerVisualiserWidget->update();
            });


    connect(this->computeFiberSurfaceButton, &QPushButton::clicked, this, [this]() {
            this->plotWidget->recomputeFiberSurface = true;
            this->plotWidget->update();
            this->update();
            });

    connect(this->computeFiberSurfaceFeatureButton, &QPushButton::clicked, this, [this]() {
            this->plotWidget->recomputeFiberSurfaceFeature = true;
            this->plotWidget->update();
            this->update();
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

    connect(checkboxShowFeatures, &QCheckBox::toggled, [=](bool checked) {
            this->tracerVisualiserWidget->drawFeatureSurfaces = checked;
            this->tracerVisualiserWidget->generateDisplayList();
            this->tracerVisualiserWidget->update();
            });

    //connect(checkboxShowFaces, &QCheckBox::toggled, [=](bool checked) {
            //this->tracerVisualiserWidget->drawFaces = checked;

            //this->tracerVisualiserWidget->update();
            //});

    //connect(checkboxShowTraces, &QCheckBox::toggled, [=](bool checked) {

            //this->plotWidget->fiberPointsTraces.clear();
            //this->plotWidget->fiberPointsTraces.shrink_to_fit();
            //this->plotWidget->update();

            //this->tracerVisualiserWidget->clearFibers = !this->tracerVisualiserWidget->clearFibers;
            //this->tracerVisualiserWidget->updateFiber({});
            //this->tracerVisualiserWidget->update();
            //});


    connect(this->buttonShowTraces, &QPushButton::clicked, this, [this]() {
            this->plotWidget->fiberPointsTraces.clear();
            this->plotWidget->fiberPointsTraces.shrink_to_fit();
            this->plotWidget->update();

            this->tracerVisualiserWidget->clearFibers();
            this->tracerVisualiserWidget->update();
            });

    connect(this->fiberOpacitySlider, &QSlider::valueChanged, plotWidget, [=]() {
        this->tracerVisualiserWidget->generateDisplayList();
        this->tracerVisualiserWidget->update();
    });

    connect(this->fiberSurfaceOpacitySlider, &QSlider::valueChanged, plotWidget, [=]() {
        this->tracerVisualiserWidget->generateDisplayList();
        this->tracerVisualiserWidget->update();
    });

    connect(this->featureSurfaceOpacitySlider, &QSlider::valueChanged, plotWidget, [=]() {
        this->tracerVisualiserWidget->generateDisplayList();
        this->tracerVisualiserWidget->update();
    });
}

TracerVisualiserWindow::~TracerVisualiserWindow()
{
    delete plotWidget;
    delete windowLayout;
    delete tracerVisualiserWidget;
}

void TracerVisualiserWindow::updateSelectedSheets(const std::set<int>& selectedSheetIds)
{
    if (selectedSheetIds.empty())
    {
        infoBox->setPlainText("");
        return;
    }

    std::ostringstream info;
    info << "Selected sheets: ";
    for (const int id : selectedSheetIds)
    {
        info << id << " (area " << std::fixed << std::setprecision(2) << 100.0 * data.reebSpace2.sheetAreaProportion[id] << ")  ";
    }
    infoBox->setPlainText(QString::fromStdString(info.str()));

    this->update();
}
