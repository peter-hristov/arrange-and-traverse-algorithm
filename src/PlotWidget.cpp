#include "./CGALTypedefs.h"

#include <QApplication>
#include <QDesktopWidget>
#include <QGraphicsScene>
#include <QLabel>
#include <QMainWindow>
#include <QMessageBox>
#include <QObject>
#include <QPainter>
#include <QVector>
#include <QtGui>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <qcolor.h>
#include <qnamespace.h>
#include <qpoint.h>
#include <qtransform.h>
#include <utility>


#include "./PlotWidget.h"
#include "./Timer.h"
#include "./Fiber.h"
#include "./utility/Geometry.h"
#include "./TracerVisualiserWindow.h"

using namespace std;

const bool DRAW_GRIDLINES = true;


PlotWidget::PlotWidget(QWidget *parent, Data &_data)
  : QWidget(parent)
  , data(_data)
{
    setMouseTracking(true);
    setEnabled(true);

    paddedMinF = data.tetMesh.minF - paddingScalingFactor * (data.tetMesh.maxF - data.tetMesh.minF);
    paddedMaxF = data.tetMesh.maxF + paddingScalingFactor * (data.tetMesh.maxF - data.tetMesh.minF);
    paddedMinG = data.tetMesh.minG - paddingScalingFactor * (data.tetMesh.maxG - data.tetMesh.minG);
    paddedMaxG = data.tetMesh.maxG + paddingScalingFactor * (data.tetMesh.maxG - data.tetMesh.minG);


    paddedMinF = -0.006111;
    paddedMaxF = 0.047295;

    paddedMinG = 0.09477;
    paddedMaxG = 0.135989;
}

void PlotWidget::mousePressEvent(QMouseEvent* event)
{
    if (event->button() == Qt::LeftButton) 
    {
        mousePointInitialPos = event->localPos();
        mousePoint = mousePointInitialPos;
        dragging = false;
        recomputeFiber = true;
        update();
    }
}

void PlotWidget::mouseMoveEvent(QMouseEvent* event)
{
    if (event->buttons() & Qt::LeftButton)
    {
        QPointF currentPos = event->localPos();

        if (!dragging)
        {
            if ((currentPos - mousePointInitialPos).manhattanLength() > dragThreshold)
            {
                dragging = true;
            }
        }

        if (dragging)
        {
            mousePoint = currentPos;
            recomputeFiber = true;
            update();
        }
    }
}

void PlotWidget::drawReebSpaceBackground(QPainter &p)
{
    //for (const auto &[faceHandle, componentIds] : data.reebSpace2.correspondenceGraph)
    for (auto faceHandle = data.singularArrangement.arr.faces_begin(); faceHandle != data.singularArrangement.arr.faces_end(); ++faceHandle) 
    {
        if (faceHandle->is_unbounded()) { continue; }



        //
        // Assemble the polygon
        //
        QVector<QPointF> points;

        //printf("\nFace with ID = %d has these points\n", data.singularArrangement.arrangementFacesIdices[faceHandle]);
        typename Arrangement_2::Ccb_halfedge_const_circulator circ = faceHandle->outer_ccb();
        typename Arrangement_2::Ccb_halfedge_const_circulator curr = circ;
        do {
            typename Arrangement_2::Halfedge_const_handle e = curr;

            // Get point from CGAL (and convert to double )
            const float u = CGAL::to_double(e->source()->point().x());
            const float v = CGAL::to_double(e->source()->point().y());

            // Add to the polygon
            points << rescalePoint(u, v);

            //std::cout << "   (" << e->source()->point() << ")  -> " << "(" << e->target()->point() << ")" << std::endl;
        } while (++curr != circ);

        QPolygonF qPolygon(points);


        //
        // Draw a polygon per face
        //
        for (const int componentId : data.reebSpace2.correspondenceGraph[faceHandle->data()])
        {
            const int sheetId = data.reebSpace2.correspondenceGraphDS.parent[componentId];

            const array<float, 3> colorF = fiber::fiberColours[sheetId % fiber::fiberColours.size()];

            p.setBrush(QColor::fromRgbF(colorF[0], colorF[1], colorF[2], 0.292f));
            p.setPen(Qt::NoPen);
            p.drawPolygon(qPolygon);
        }
    }









    //for (const auto &[faceHandle, componentIds] : data.reebSpace2.correspondenceGraph)

    // Draw polygonds from the regular reeb space
    //for (const auto &[sheetId, polygon] : data.reebSpace.sheetPolygon)
    //{
        //QVector<QPointF> points;

        //// If the sheet is incomplete, the polygon will not be corret, just draww all the faces manually
        //if (data.reebSpace.incompleteSheets.contains(sheetId))
        //{

            //// Loop through all faces to see which ones are in the sheet
            //for (auto f = data.arrangement.arr.faces_begin(); f != data.arrangement.arr.faces_end(); ++f) 
            //{
                //const int currentFaceID = data.arrangement.arrangementFacesIdices[f];

                //// For each fiber component in the face, see if one of those is in our sheet
                //for (const auto &[triangleId, fiberComponentId] : this->data.reebSpace.fiberSeeds[currentFaceID])
                //{
                    //const int componentSheetId = data.reebSpace.correspondenceGraph.findElement({currentFaceID, fiberComponentId});

                    //// Now we can add the polygon
                    //if (componentSheetId == sheetId)
                    //{
                        //typename Arrangement_2::Ccb_halfedge_const_circulator circ = f->outer_ccb();
                        //typename Arrangement_2::Ccb_halfedge_const_circulator curr = circ;
                        //do {
                            //typename Arrangement_2::Halfedge_const_handle e = curr;

                            //// Get point from CGAL (and convert to double )
                            //const float u = CGAL::to_double(e->source()->point().x());
                            //const float v = CGAL::to_double(e->source()->point().y());

                            //// Add to the polygon
                            //points << rescalePoint(u, v);

                            ////std::cout << "   (" << e->source()->point() << ")  -> " << "(" << e->target()->point() << ")" << std::endl;
                        //} while (++curr != circ);
                    //}
                //}
            //}
        //}

        //// If the sheet is no incomplete, the polygon is valid, draw the directly
        //else
        //{
            //for (const CartesianPoint &point : polygon) 
            //{
                //// Get point from CGAL (and convert to double )
                //const float u = point.x();
                //const float v = point.y();

                //// Add to the polygon
                //points << rescalePoint(u, v);
            //}
        //}

        //QPolygonF qPolygon(points);

        //const int colourID = data.reebSpace.sheetConsequitiveIndices[sheetId];
        //const array<float, 3> colorF = fiber::fiberColours[colourID % fiber::fiberColours.size()];

        //p.setBrush(QColor::fromRgbF(colorF[0], colorF[1], colorF[2], 0.392f));
        //p.setPen(Qt::NoPen);
        //p.drawPolygon(qPolygon);
    //}









    //for (int i = 0 ; i < this->data.fiberSeeds.size() ; i++) 
    //{
        //const int currentFaceID = i;

        //// For each fiber component
        //for (int j = 0 ; j < this->data.fiberSeeds[i].size() ; j++) 
        //{
            ////const auto &[]
            //const auto &[triangleId, fiberComponentId] = this->data.fiberSeeds[i][j];
            //const int sheetId = this->data.reebSpace.findTriangle({currentFaceID, fiberComponentId});
            //const int colourID = data.sheetConsequitiveIndices[sheetId];
            //const vector<float> colorF = data.fiberColours[colourID % data.fiberColours.size()];

            //p.setBrush(QColor::fromRgbF(colorF[0], colorF[1], colorF[2], 0.392f));
            //p.setPen(Qt::NoPen);
            //p.drawPolygon(this->arrangementPolygons[i]);
        //}
    //}

    // We assume that the fiber seeds per face are sorted by their sheetId
    //for (int i = 0 ; i < this->data.fiberSeeds.size() ; i++) 
    //{
        //const int currentFaceID = i;

        //// For each fiber component
        //for (int j = 0 ; j < this->data.fiberSeeds[i].size() ; j++) 
        //{
            ////const auto &[]
            //const auto &[triangleId, fiberComponentId] = this->data.fiberSeeds[i][j];
            //const int sheetId = this->data.reebSpace.findTriangle({currentFaceID, fiberComponentId});
            //const int colourID = data.sheetConsequitiveIndices[sheetId];
            //const vector<float> colorF = data.fiberColours[colourID % data.fiberColours.size()];

            //p.setBrush(QColor::fromRgbF(colorF[0], colorF[1], colorF[2], 0.392f));
            //p.setPen(Qt::NoPen);
            //p.drawPolygon(this->arrangementPolygons[i]);
        //}
    //}


    //
    // Draw the polygons of the Reeb space
    //
    //for (int i = 0 ; i < this->arrangementPolygons.size() ; i++) 
    //{
        //// Set the random color for filling the polygon
        //p.setBrush(this->arrangementPolygonColours[i]);
        //p.setPen(Qt::NoPen);

        //// Draw the filled polygon with the random color
        //p.drawPolygon(this->arrangementPolygons[i]);

        ////qDebug() << "New polygon --- ";
        ////for (const QPoint& point : this->arrangementPolygons[i]) {
            ////qDebug() << "(" << point.x() << ", " << point.y() << ")";
        ////}
    //}

    //
    // Draw the Jacobi set
    //
    //for (const auto &[edge, type] : data.reebSpace.jacobiType)
    //{
        //if (type != 1)
        //{
            //if (type == 0)
            //{
                ////p.setPen(QPen(Qt::black, 1, Qt::DashLine));
                //p.setPen(QPen(Qt::black, 0.2));
            //}
            //else
            //{
                //p.setPen(QPen(Qt::black, 0.2));
            //}

            //float x1 = (resolution / (data.tetMesh.maxF - data.tetMesh.minF)) * (this->data.tetMesh.vertexCoordinatesF[edge.first] - data.tetMesh.minF);
            //float y1 = (resolution / (data.tetMesh.maxG - data.tetMesh.minG)) * (this->data.tetMesh.vertexCoordinatesG[edge.first] - data.tetMesh.minG);

            //float x2 = (resolution / (data.tetMesh.maxF - data.tetMesh.minF)) * (this->data.tetMesh.vertexCoordinatesF[edge.second] - data.tetMesh.minF);
            //float y2 = (resolution / (data.tetMesh.maxG - data.tetMesh.minG)) * (this->data.tetMesh.vertexCoordinatesG[edge.second] - data.tetMesh.minG);

            //p.setRenderHint(QPainter::Antialiasing, true);
            //p.drawLine(x1, y1, x2, y2);
        //}
    //}






    // Draw all edges
    //
    for (const auto &[edge, type] : data.tetMesh.edgeSingularTypes)
    {
        const float u1 = this->data.tetMesh.vertexCoordinatesF[edge[0]];
        const float v1 = this->data.tetMesh.vertexCoordinatesG[edge[0]];

        const float u2 = this->data.tetMesh.vertexCoordinatesF[edge[1]];
        const float v2 = this->data.tetMesh.vertexCoordinatesG[edge[1]];

        if (type == 0)
        {
            p.setPen(QPen(Qt::black, 4.2, Qt::DashLine));
        }
        else if (type == 1)
        {
            //continue;
            //p.setPen(QPen(Qt::black, 3.2, Qt::SolidLine));
            p.setPen(QPen(Qt::black, 3.0, Qt::SolidLine));
        }
        else
        {
            p.setPen(QPen(Qt::black, 10.2, Qt::SolidLine));
        }

        p.setRenderHint(QPainter::Antialiasing, true);
        p.drawLine(rescalePoint(u1, v1), rescalePoint(u2, v2));
    }

    //// Draw all the vertex coordinates
    //for(size_t i = 0 ; i <  this->data.tetMesh.vertexCoordinatesF.size() ; i++)
    //{
        //float u = this->data.tetMesh.vertexCoordinatesF[i];
        //float v = this->data.tetMesh.vertexCoordinatesG[i];

        //p.setPen(QPen(Qt::black, 6, Qt::SolidLine));
        //p.setBrush(Qt::white);           // Fill color
        //p.drawEllipse(rescalePoint(u, v), 20, 20);

        //QFont font = p.font();
        //font.setPointSize(70);
        //p.setFont(font);
        //p.drawText(rescalePoint(u, v), QString::number(i));
    //}


    // Draw all singular vertices
    //for (auto vit = data.singularArrangement.arr.vertices_begin(); vit != data.singularArrangement.arr.vertices_end(); ++vit) 
    //{
        //const float u = CGAL::to_double(vit->point().x());
        //const float v = CGAL::to_double(vit->point().y());

        //p.setPen(QPen(Qt::black, 6, Qt::SolidLine));
        //p.setBrush(Qt::white);           // Fill color
        //p.drawEllipse(rescalePoint(u, v), 20, 20);

        //QFont font = p.font();
        //font.setPointSize(70);
        //p.setFont(font);
        ////p.drawText(rescalePoint(u, v), QString::number(1));
    //}

    //// Draw all singular arrangement half-edges
    //for (const auto &[halfEdge, vertices] : data.singularArrangement.halfEdgePoints)
    //{
        //for (const auto &vertex : vertices)
        //{
            //const float u = CGAL::to_double(vertex.x());
            //const float v = CGAL::to_double(vertex.y());

            //p.setPen(QPen(Qt::black, 5, Qt::SolidLine));
            //p.setBrush(Qt::white);           // Fill color
            //p.drawEllipse(rescalePoint(u, v), 15, 15);
        //}
    //}
}

void PlotWidget::generateStaticReebSpaceCache()
{
    if (!staticReebSpaceCache) 
    {
        staticReebSpaceCache = std::make_unique<QPixmap>(resolution, resolution);
        staticReebSpaceCache->fill(Qt::white);
        qDebug() << "Redrawing REEB SPACE ...";
        QPainter p(staticReebSpaceCache.get());
        this->drawReebSpaceBackground(p);
    }
}

void PlotWidget::resizeEvent(QResizeEvent* event)
{
    generateStaticReebSpaceCache();
}

void PlotWidget::paintEvent(QPaintEvent*)
{
    QPainter p(this);
    p.setWindow(QRect(0, 0, resolution, resolution));
    // p.setViewport(QRect(0, 0, resolution, resolution));

    p.save();
    p.setTransform(QTransform(1., 0., 0., -1., 0., resolution));
    p.setPen(Qt::gray);

    generateStaticReebSpaceCache();
    p.drawPixmap(0, 0, *(this->staticReebSpaceCache));

    QPointF fiberPoint = p.combinedTransform().inverted().map(mousePoint);

    // Draw fiber point
    auto penGrey = QPen(QColor(0, 0, 0, 250));
    penGrey.setWidthF(8.0);
    p.setPen(penGrey);
    p.drawEllipse(QPointF(fiberPoint.x(), fiberPoint.y()), sphereRadius, sphereRadius);

    // Crosshair around fiber point
    penGrey.setWidthF(1.0);
    p.setPen(penGrey);
    p.drawLine(fiberPoint.x(), fiberPoint.y() - resolution, fiberPoint.x(), fiberPoint.y() + resolution);
    p.drawLine(fiberPoint.x() - resolution, fiberPoint.y(), fiberPoint.x() + resolution, fiberPoint.y());

    if (this->recomputeFiber == true)
    {
        this->recomputeFiber = false;

        const float u = this->paddedMinF + (fiberPoint.x() / resolution) * (this->paddedMaxF - this->paddedMinF);
        const float v = this->paddedMinG + (fiberPoint.y() / resolution) * (this->paddedMaxG - this->paddedMinG);

        qDebug() << "Computing fiber (" << u << ", " << v << ")";

        //const std::vector<FiberPoint> fiber = fiber::computeFiber(data.tetMesh, data.arrangement, data.reebSpace, {u, v}, -1);
        //const std::vector<FiberPoint> fiber = fiber::computeFiberFromFiberGraph(data.tetMesh, data.singularArrangement, data.reebSpace2, {u, v});
        const std::vector<FiberPoint> fiber = fiber::computeFiberSAT(data.tetMesh, data.singularArrangement, data.reebSpace2, {u, v});

        sibling->updateFiber(fiber);
    }

    p.restore();
    drawAxisLabels(p);
}



QPointF PlotWidget::rescalePoint(const float &u, const GLfloat &v)
{
    const float rescaledU = (resolution / (paddedMaxF - paddedMinF)) * (u - paddedMinF);
    const float rescaledV = (resolution / (paddedMaxG - paddedMinG)) * (v - paddedMinG);
    return QPointF(rescaledU, rescaledV);
}

void PlotWidget::drawAxisLabels(QPainter& p)
{
    auto penGrey = QPen(QColor(0, 0, 0, 50));
    penGrey.setWidth(0.5);

    auto penBlack = QPen(Qt::black);
    penBlack.setWidth(0.5);
    p.setPen(penBlack);

    float boxOffset = 15;

    // X Axis
    p.drawLine(boxOffset, resolution - boxOffset, resolution - boxOffset + 100, resolution - boxOffset);
    // Y Axis
    p.drawLine(boxOffset, resolution - boxOffset, boxOffset, boxOffset - 100);

    // Write out numbers
    QFont font = p.font();
    font.setPixelSize(6);
    // font.setWeight(20);
    p.setFont(font);

    //
    // Write out numbers on axis
    //
    float step = resolution / 15;

    float xMin = paddedMinF;
    float yMin = paddedMinG;
    float xMax = paddedMaxF;
    float yMax = paddedMaxG;
    float xRange = xMax - xMin;
    float yRange = yMax - yMin;
    float a = 10.0;
    int minSteps = 6;
    float xStepSize = pow(a, round(log(xRange) / log(a)) - 1);
    float yStepSize = pow(a, round(log(yRange) / log(a)) - 1);

    if (xRange / xStepSize < minSteps) {
        xStepSize /= 4.0;
    }
    if (yRange / yStepSize < minSteps) {
        yStepSize /= 4.0;
    }

    // ranges for labels and gridlines
    float xMinPlot = ceil(xMin / xStepSize) * xStepSize;
    float yMinPlot = ceil(yMin / yStepSize) * yStepSize;
    float xMaxPlot = floor(xMax / xStepSize) * xStepSize;
    float yMaxPlot = floor(yMax / yStepSize) * yStepSize;

    //
    // Draw Carthesian Grid and labels
    //

    // X labels and gridlines
    float xCurrent = xMinPlot;
    do {
        int i = int(resolution * (xCurrent - xMin) / xRange);
        if (i < boxOffset || i > (resolution - boxOffset)) {
            // skip this point
        } else {
            // grid line
            if (DRAW_GRIDLINES) {
                p.setPen(penGrey);
                p.drawLine(i, resolution - boxOffset - 20000, i, resolution - boxOffset);
            }
            // label
            p.setPen(penBlack);
            p.drawText(i - 5, resolution - boxOffset - 5, QString::number(xCurrent).mid(0, 6));
        }
        xCurrent += xStepSize;
        xCurrent = round(xCurrent / xStepSize) * xStepSize;
    } while (xCurrent < xMaxPlot);

    // Y labels and gridlines
    float yCurrent = yMinPlot;
    do {
        int i = int(resolution * (1.0 - (yCurrent - yMin) / yRange));
        if (i < boxOffset || i > (resolution - boxOffset)) {
            // skip this point
        } else {
            // grid line
            if (DRAW_GRIDLINES) {
                p.setPen(penGrey);
                p.drawLine(boxOffset, i, boxOffset + 2000, i);
            }
            // label
            p.setPen(penBlack);
            p.drawText(boxOffset + 5, i + 2, QString::number(yCurrent).mid(0, 6));
        }
        yCurrent += yStepSize;
        yCurrent = round(yCurrent / yStepSize) * yStepSize;
    } while (yCurrent < yMaxPlot);

    //
    // Draw custom lines
    //
    auto pen3 = QPen(QColor(0, 100, 0, 100));
    pen3.setWidth(1);
    p.setPen(pen3);

    // Vertical (constant X)
    for (const auto number : this->verticalLineNumbers) {
        if (paddedMinF < number && number < paddedMaxF) {
            float ratio = (number - paddedMinF) / ((paddedMaxF - paddedMinF));
            float plotPosition = ratio * (resolution);
            p.drawLine(plotPosition, resolution - boxOffset - 20000, plotPosition, resolution - boxOffset);
        }
    }

    // Horizontal (constant Y)
    for (const auto number : this->horizontalLineNumbers) {
        if (paddedMinG < number && number < paddedMaxG) {
            // Invert the y axis
            float ratio = 1.0 - (number - paddedMinG) / ((paddedMaxG - paddedMinG));
            float plotPosition = ratio * (resolution);
            p.drawLine(boxOffset, plotPosition, boxOffset + 2000, plotPosition);
        }
    }

    // Draw Labels
    p.setPen(penBlack);

    font.setPixelSize(7);
    p.setFont(font);

    // x label
    p.drawText(resolution / 2 - 30,
            resolution - boxOffset + 10,
            QString::fromStdString(data.tetMesh.longnameF) + " (" + QString::fromStdString(data.tetMesh.units) +
            ")");

    p.translate(boxOffset - 5, resolution / 2 + 20);
    p.rotate(-90);

    // y label
    p.drawText(
            0, 0, QString::fromStdString(data.tetMesh.longnameG) + " (" + QString::fromStdString(data.tetMesh.units) + ")");
}
