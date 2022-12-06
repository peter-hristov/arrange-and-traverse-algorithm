#include <iostream>
#include <utility>
#include <QCheckBox>
#include <QHBoxLayout>


#include "./TracerVisualiserWindow.h"
#include "src/Data.h"

using namespace std;
void
TracerVisualiserWindow::keyPressEvent(QKeyEvent* event)
{
    if (this->data->mousePoints.empty())
        return;

    int speed = 2;

    if (event->key() == Qt::Key_U) {
        QApplication::quit();
    }

    if (event->key() == Qt::Key_I) {
        this->data->mousePoints[0].setY(this->data->mousePoints[0].y() - speed);
        this->update();
    }
    if (event->key() == Qt::Key_J) {
        this->data->mousePoints[0].setX(this->data->mousePoints[0].x() - speed);
        this->update();
    }
    if (event->key() == Qt::Key_K) {
        this->data->mousePoints[0].setY(this->data->mousePoints[0].y() + speed);
        this->update();
    }
    if (event->key() == Qt::Key_L) {
        this->data->mousePoints[0].setX(this->data->mousePoints[0].x() + speed);
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

    if (event->key() == Qt::Key_O) {
        this->tracerVisualiserWidget->clearFibers = !this->tracerVisualiserWidget->clearFibers;
        this->tracerVisualiserWidget->update();
    }
}

TracerVisualiserWindow::TracerVisualiserWindow(QWidget* parent,
                                               Data* _data
                                               )
  : QWidget(parent)
{
    this->data = _data;

    // Initialise Widgets
    plotWidget = new PlotWidget(this, data, "none");
    tracerVisualiserWidget = new TracerVisualiserWidget(this, plotWidget, data);

    checkboxShowVertices = new QCheckBox("Show Vertices");
    checkboxShowVertices->setChecked(true);

    checkboxShowEdges = new QCheckBox("Show Edges");
    checkboxShowEdges->setChecked(true);

    checkboxShowFaces = new QCheckBox("Show Faces");
    checkboxShowFaces->setChecked(true);


    vertexOpacitySlider = new QSlider(Qt::Horizontal);
    vertexOpacitySlider->setValue(this->tracerVisualiserWidget->vertexOpacity * 100);

    edgeOpacitySlider = new QSlider(Qt::Horizontal);
    edgeOpacitySlider->setValue(this->tracerVisualiserWidget->edgeOpacity * 100);

    faceOpacitySlider = new QSlider(Qt::Horizontal);
    faceOpacitySlider->setValue(this->tracerVisualiserWidget->faceOpacity * 100);

    fakeSlider = new QSlider(Qt::Horizontal);
    fakeSlider->setMaximum(1000);
    fakeSlider->setTracking(false);

    fakeSlider2 = new QSlider(Qt::Horizontal);
    fakeSlider2->setMaximum(1000);
    fakeSlider2->setTracking(false);

    checkbox2 = new QCheckBox("Function Value");

    //
    // Layouts
    //
    optionsLayout = new QGridLayout();
    optionsLayout2 = new QGridLayout();

    auto rowOneLayout = new QGridLayout();
    rowOneLayout->addWidget(checkboxShowVertices,0, 0);
    rowOneLayout->addWidget(vertexOpacitySlider, 0, 1);

    rowOneLayout->addWidget(checkboxShowEdges,1, 0);
    rowOneLayout->addWidget(edgeOpacitySlider, 1, 1);

    rowOneLayout->addWidget(checkboxShowFaces,2, 0);
    rowOneLayout->addWidget(faceOpacitySlider, 2, 1);


    optionsLayout->addLayout(rowOneLayout, 0, 0);

    optionsLayout2->addWidget(checkbox2, 0, 0);
    optionsLayout2->addWidget(fakeSlider, 0, 1);
    optionsLayout2->addWidget(fakeSlider2, 1, 1);


    // Set up layout
    windowLayout = new QGridLayout(this);
    windowLayout->addWidget(tracerVisualiserWidget, 0, 0);
    windowLayout->addWidget(plotWidget, 0, 1);

    windowLayout->addLayout(optionsLayout, 1, 0);
    windowLayout->addLayout(optionsLayout2, 1, 1);


    connect(checkboxShowVertices, &QCheckBox::toggled, [=](bool checked) {
            this->tracerVisualiserWidget->drawVertices = checked;
            this->tracerVisualiserWidget->update();
            });

    connect(checkboxShowEdges, &QCheckBox::toggled, [=](bool checked) {
            this->tracerVisualiserWidget->drawEdges = checked;
            this->tracerVisualiserWidget->update();
            });

    connect(checkboxShowFaces, &QCheckBox::toggled, [=](bool checked) {
            this->tracerVisualiserWidget->drawFaces = checked;

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

    connect(this->fakeSlider, &QSlider::valueChanged, plotWidget, [=]() {
            double sliderRatio =  static_cast<double>(this->fakeSlider->value()) / 100.0;
            double isovalue = this->data->minF + sliderRatio * (this->data->maxF - this->data->minF);

            this->tracerVisualiserWidget->generateDisplayListTriangles(isovalue, this->data->vertexCoordinatesF, 1);
            this->tracerVisualiserWidget->update();

            //cout << "Isovalue is " isovalue << endl;
            this->checkbox2->setText(QString::number(isovalue));

        //this->tracerVisualiserWidget->edgeOpacity = static_cast<double>(this->fakeSlider->value()) / 100.0;
        //this->tracerVisualiserWidget->update();
    });

    connect(this->fakeSlider2, &QSlider::valueChanged, plotWidget, [=]() {

            double sliderRatio =  static_cast<double>(this->fakeSlider2->value()) / 100.0;
            double isovalue = this->data->minG + sliderRatio * (this->data->maxG - this->data->minG);

            this->tracerVisualiserWidget->generateDisplayListTriangles(isovalue, this->data->vertexCoordinatesG, 2);
            this->tracerVisualiserWidget->update();

        //this->tracerVisualiserWidget->edgeOpacity = static_cast<double>(this->fakeSlider->value()) / 100.0;
        //this->tracerVisualiserWidget->update();
    });
}

TracerVisualiserWindow::~TracerVisualiserWindow()
{
    delete plotWidget;
    delete windowLayout;
    delete tracerVisualiserWidget;
}
