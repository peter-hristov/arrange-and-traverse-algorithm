#pragma once

#include <QKeyEvent>
#include <QWidget>
#include <QGridLayout>
#include <QSlider>
#include <QComboBox>
#include <QSpinBox>
#include <QPushButton>


#include "./Data.h"
#include "./PlotWidget.h"
#include "./TracerVisualiserWidget.h"

class TracerVisualiserWindow : public QWidget
{
    Q_OBJECT
  public:
    Data &data;
    TracerVisualiserWindow(QWidget*, Data &);

    ~TracerVisualiserWindow();

    void keyPressEvent(QKeyEvent* event);


    PlotWidget* plotWidget;
    TracerVisualiserWidget* tracerVisualiserWidget;

    QGridLayout* windowLayout;

    QGridLayout *optionsLayout;
    QGridLayout * optionsLayout2;

    QCheckBox *checkboxShowFibers;
    QCheckBox *checkboxShowFiberSurfaces;
    QCheckBox *checkboxShowFaces;

    QSlider *faceOpacitySlider;
    QSlider *edgeOpacitySlider;
    QSlider *vertexOpacitySlider;

    QSlider *fakeSlider;
    QPushButton *computeTracedFiberSurfaceButton; 
    QPushButton *computeFiberSurfaceButton; 
    QPushButton *computeFiberSurfaceFeatureButton; 


    QPushButton *clearFibersButton; 
    QPushButton *clearFiberSurfaceButton; 
    QPushButton *clearSelectedSheetsButton; 
    QPushButton *clearAllButton; 


    QCheckBox *checkboxShowTraces;

    QSpinBox* spinBoxAddSheet;
    QPushButton* buttonAddSheet;

    QSpinBox* spinBoxAddTopSheets;
    QPushButton* buttonAddTopSheets;
};
