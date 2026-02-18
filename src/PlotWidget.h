#pragma once

#include <QApplication>
#include <QDesktopWidget>
#include <QLabel>
#include <QMainWindow>
#include <QMessageBox>
#include <QObject>
#include <QPainter>
#include <QVector>
#include <QtGui>
#include <QtMath>
#include <algorithm>
#include <assert.h>
#include <iostream>
#include <limits>
#include <map>
#include <qpoint.h>
#include <qvector.h>
#include <utility>
#include <cmath>

#include "./Data.h"
#include "./TracerVisualiserWidget.h"

class PlotWidget : public QWidget
{
    Q_OBJECT

    public:
        PlotWidget(QWidget*, Data&);

        // Add some padding to the bounding box of the image of the domain so it does't fill the widget.
        float paddingScalingFactor = 0.1;
        float paddedMinF;
        float paddedMaxF;
        float paddedMinG;
        float paddedMaxG;

        // Cache the Reeb space so we don't have to draw it every referesh
        std::unique_ptr<QPixmap> staticReebSpaceCache;
        void generateStaticReebSpaceCache();
        void drawReebSpaceBackground(QPainter &p);

        // Indicates that we need to recompute the fiber, clicked/dragged the mouse/fiber point
        bool recomputeFiber = false;

        Data &data;
        TracerVisualiserWidget *sibling;

        const float resolution = 2000;

        QPointF rescalePoint(const float&, const float&);


        std::vector<float> verticalLineNumbers;
        std::vector<float> horizontalLineNumbers;

        void paintEvent(QPaintEvent* event);
        void resizeEvent(QResizeEvent* event);

        // Fiber point stuff

        // Initial position of where we click in the widget
        QPointF mousePointInitialPos;

        // Position of where we click and drag on the widget
        QPointF mousePoint;

        std::vector<QPointF> controlPoints;

        // Whether the mouse is being dragged over the widget
        bool dragging = false;

        const qreal dragThreshold = 4.0;

        // This is the radius of the sphere around the fiber point
        const float sphereRadius = this->resolution / 100.0;
        const float controlPointRadious = this->resolution / 200.0;

    protected:
        void mouseMoveEvent(QMouseEvent* event);
        void mousePressEvent(QMouseEvent* event);
        void drawInteriorPointsImages(QPainter& p);
        void drawAxisLabels(QPainter& p);
};
