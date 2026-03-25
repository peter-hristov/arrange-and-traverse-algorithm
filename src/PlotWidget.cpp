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

#include <filesystem>
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
#include "./io.h"
#include "./Fiber.h"
#include "./utility/Geometry.h"
#include "./TracerVisualiserWindow.h"
#include "src/SurfaceMesh.h"

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

    // min -0.0140185 ,  0.106133 
    // max 0.0405787 ,  0.125916

    //paddedMinF = -0.0180185;
    //paddedMinG = 0.106133;

    //paddedMaxF = 0.0405787;
    //paddedMaxG = 0.125916;
}




void PlotWidget::mousePressEvent(QMouseEvent* event)
{
    if (event->button() == Qt::RightButton) 
    {
        this->controlPoints.push_back(event->localPos());
        update();
    }

    else if (event->button() == Qt::LeftButton && event->modifiers() == Qt::ShiftModifier)
    {
        controlPointSheetSelection = event->localPos();
        update();
    }
    else if (event->button() == Qt::LeftButton) 
    {
        mousePointInitialPos = event->localPos();
        mousePoint = mousePointInitialPos;

        if (false == this->sibling->clearFibers)
        {
            fiberPointsTraces.push_back({mousePoint});
        }

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
            if (false == this->sibling->clearFibers)
            {
                if (fiberPointsTraces.size() == 0)
                {
                    fiberPointsTraces.push_back({});
                }
                fiberPointsTraces.back().push_back(mousePoint);
            }
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

            //if (!this->sibling->selectedSheetIds.empty() && !this->sibling->selectedSheetIds.contains(sheetId))
            //{
                //continue;
            //}

            //if (desiredSheetId != -1 && sheetId != desiredSheetId)
            //{
                //continue;
            //}

            float alpha;

            if (this->sibling->selectedSheetIds.empty())
            {
                alpha = 0.3;
            }
            else if (this->sibling->selectedSheetIds.contains(sheetId))
            {
                alpha = 0.7;
            }
            else
            {
                alpha = 0.05;
            }

            const array<float, 3> colorF = fiber::fiberColours[sheetId % fiber::fiberColours.size()];

            p.setBrush(QColor::fromRgbF(colorF[0], colorF[1], colorF[2], alpha));
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
    //for (const auto &[edge, type] : data.tetMesh.edgeSingularTypes)
    //{
        //const float u1 = this->data.tetMesh.vertexCoordinatesF[edge[0]];
        //const float v1 = this->data.tetMesh.vertexCoordinatesG[edge[0]];

        //const float u2 = this->data.tetMesh.vertexCoordinatesF[edge[1]];
        //const float v2 = this->data.tetMesh.vertexCoordinatesG[edge[1]];

        //if (type == 0)
        //{
            //p.setPen(QPen(Qt::black, 4.2, Qt::DashLine));
        //}
        //else if (type == 1)
        //{
            ////continue;
            ////p.setPen(QPen(Qt::black, 3.2, Qt::SolidLine));
            //p.setPen(QPen(Qt::black, 3.0, Qt::SolidLine));
        //}
        //else
        //{
            //p.setPen(QPen(Qt::black, 10.2, Qt::SolidLine));
        //}

        //p.setRenderHint(QPainter::Antialiasing, true);
        //p.drawLine(rescalePoint(u1, v1), rescalePoint(u2, v2));
    //}

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
        //qDebug() << "Redrawing REEB SPACE ...";
        QPainter p(staticReebSpaceCache.get());

        Timer::start();
        this->drawReebSpaceBackground(p);
        Timer::stop("Rendered Reeb space                    :");
    }
}

void PlotWidget::resizeEvent(QResizeEvent* event)
{
    generateStaticReebSpaceCache();
}

void PlotWidget::saveToFile(const std::string& filename)
{
    std::filesystem::path filePath(filename);
    if (filePath.has_parent_path())
        std::filesystem::create_directories(filePath.parent_path());

    // Force a repaint so everything is current, then grab
    this->repaint();
    QPixmap pixmap = this->grab();
    pixmap = pixmap.scaled(1000, 1000, Qt::IgnoreAspectRatio, Qt::SmoothTransformation);

    pixmap.save(QString::fromStdString(filename));
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

    auto penBlack = QPen(QColor(0, 0, 0, 250));
    penBlack.setWidthF(8.0);
    p.setPen(penBlack);

    QPointF fiberPoint = p.combinedTransform().inverted().map(mousePoint);
    p.drawEllipse(fiberPoint, sphereRadius, sphereRadius);

    // Crosshair around fiber point
    //p.drawLine(fiberPoint.x(), fiberPoint.y() - resolution, fiberPoint.x(), fiberPoint.y() + resolution);
    //p.drawLine(fiberPoint.x() - resolution, fiberPoint.y(), fiberPoint.x() + resolution, fiberPoint.y());

    if (false == this->sibling->clearFibers)
    {
        penBlack.setWidthF(8.0);
        p.setPen(penBlack);

        for (const QVector<QPointF> &fiberPointsTrace : this->fiberPointsTraces)
        {
            QVector<QPointF> fiberPointsTraceTransformed(fiberPointsTrace.size());
            for (int i = 0 ; i < fiberPointsTrace.size() ; i++)
            {
                const QPointF &controlPoint = fiberPointsTrace[i];
                const QPointF controlPointTransformed = p.combinedTransform().inverted().map(controlPoint);
                fiberPointsTraceTransformed[i] = controlPointTransformed;
            }
            p.drawPolyline(QPolygonF(fiberPointsTraceTransformed));
        }
    }

    if (controlPointSheetSelection.has_value())
    {
        QPointF pos = p.combinedTransform().inverted().map(controlPointSheetSelection.value());
        const float u = this->paddedMinF + (pos.x() / resolution) * (this->paddedMaxF - this->paddedMinF);
        const float v = this->paddedMinG + (pos.y() / resolution) * (this->paddedMaxG - this->paddedMinG);
        
        Face_const_handle activeFace = data.singularArrangement.getActiveFace(std::array<double, 2>{u, v});
        for (const int componentId : data.reebSpace2.correspondenceGraph[activeFace->data()])
        {
            const int sheetId = data.reebSpace2.correspondenceGraphDS.find(componentId);
            this->sibling->selectedSheetIds.insert(sheetId);

            printf("Adding sheet %d with area %.2f (which is %.2f%%).\n", sheetId, data.reebSpace2.sheetArea[sheetId], 100.0 * data.reebSpace2.sheetArea[sheetId]);
        }

        controlPointSheetSelection.reset();

        this->staticReebSpaceCache = nullptr;
        this->update();
    }




    // Draw fiber point


    // ----------------------------------------------------------------
    // Fiber Drawing
    // ----------------------------------------------------------------

    if (this->recomputeFiber == true)
    {



        //p.drawLine(fiberPoint.x(), fiberPoint.y(), fiberPoint.x(), fiberPoint.y() + 1000);

        //const Point_2 endPoint(controlPoint[0], controlPoint[1] + tetMesh.maxG + 10.0);

        this->recomputeFiber = false;

        const float u = this->paddedMinF + (fiberPoint.x() / resolution) * (this->paddedMaxF - this->paddedMinF);
        const float v = this->paddedMinG + (fiberPoint.y() / resolution) * (this->paddedMaxG - this->paddedMinG);

        //const double u = -0.0734849;
        //const double v = -0.0625043;

        //qDebug() << "Computing fiber (" << u << ", " << v << ")";

        //const std::vector<FiberPoint> fiber = fiber::computeFiber(data.tetMesh, data.arrangement, data.reebSpace, {u, v}, -1);
        //const std::vector<FiberPoint> fiber = fiber::computeFiberFromFiberGraph(data.tetMesh, data.singularArrangement, data.reebSpace2, {u, v});
        
        const std::vector<FiberPoint> fiber = fiber::computeFiberSAT(data.tetMesh, data.singularArrangement, data.reebSpace2, {u, v}, this->sibling->selectedSheetIds);

        sibling->updateFiber(fiber);
    }
    


    //if (this->recomputeFiber == true && controlPointsTransformed.size() >= 2)










    // Draw the control points and control polygon
    //
    QVector<QPointF> controlPointsTransformed(this->controlPoints.size());

    penBlack.setWidthF(8.0);
    p.setPen(penBlack);
    for (int i = 0 ; i < this->controlPoints.size() ; i++)
    {
        const QPointF &controlPoint = this->controlPoints[i];
        const QPointF controlPointTransformed = p.combinedTransform().inverted().map(controlPoint);
        controlPointsTransformed[i] = controlPointTransformed;

        p.drawEllipse(controlPointTransformed, controlPointRadious, controlPointRadious);
    }
    p.drawPolygon(QPolygonF(controlPointsTransformed));




    std::vector<std::array<double, 2>> controlPointsInternal;
    controlPointsInternal.reserve(controlPointsTransformed.size());


    for (const QPointF &controlPoint : controlPointsTransformed)
    {
        const double u = this->paddedMinF + (controlPoint.x() / resolution) * (this->paddedMaxF - this->paddedMinF);
        const double v = this->paddedMinG + (controlPoint.y() / resolution) * (this->paddedMaxG - this->paddedMinG);

        controlPointsInternal.emplace_back(std::array<double, 2>{u, v});
    }






    // ----------------------------------------------------------------
    // Fiber Surface drawing
    // ----------------------------------------------------------------


    if (this->recomputeFiberSurface == true && controlPoints.size() >= 2)
    {

        std::cout << "The control polygon is :\n";
        for (int i = 0 ; i < controlPointsInternal.size() ; i++)
        {
            qDebug() << controlPointsInternal[i][0] << ", " << controlPointsInternal[i][1];
        }

        // TTK FS
        //
        //Start point : 0.033263165761928781 0.14687746720889849 end point 0.21512107934692115 0.10335287336486487
        //Start point : 0.0013081379514506108 0.15420964485311991 end point 0.14955198551200763 -0.0054015285199353466

        // Long vertical
        // Start point : 0.0018143038060367939, 0.20159421194824584 end point -0.00063704487815769751, -0.19743372148439267
        //controlPointsInternal = {{0.0018143038060367939, 0.20159421194824584}, {-0.00063704487815769751, -0.19743372148439267}};

        //Start point : -0.021533425147216016, 0.17978634037260727 end point -0.026869776102463102, -0.17316282244046066
        //controlPointsInternal = {{-0.021533425147216016, 0.17978634037260727}, {-0.026869776102463102, -0.17316282244046066}};


        //controlPointsInternal = {{0.0, 0.2}, {0.0, -0.2}};


        // ET diagonal
        // Start point : -2.1, 0.9 end point 2.1, -0.6
        //controlPointsInternal = {{-2.1, 0.9}, {2.1, -0.6}};

        //std::vector<FiberPoint> fibersAll = fiber::computeFiberSurfaceOld(data.tetMesh, data.singularArrangement, data.reebSpace2, {controlPointsInternal[0], controlPointsInternal[1]}, desiredSheetId);

        this->data.surfaceMeshes.clear();
        this->data.surfaceMeshes.shrink_to_fit();

        std::vector<FiberPoint> fibersAll;
        if (controlPointsTransformed.size() == 2)
        {
            this->data.surfaceMeshes.push_back(fiber::computeFiberSurfaceSingularSegment(data.tetMesh, data.singularArrangement, data.reebSpace2, {controlPointsInternal[0], controlPointsInternal[1]}, desiredSheetId));

            const std::vector<FiberPoint> fibers = fiber::computeFiberPointsFromSurfaceMesh(this->data.surfaceMeshes.back(), {});
            //const std::vector<FiberPoint> fibers = fiber::computeFiberPointsFromSurfaceMesh(this->data.surfaceMeshes.back(), this->sibling->selectedSheetIds);

            fibersAll.insert(
                    fibersAll.end(), 
                    std::make_move_iterator(fibers.begin()), 
                    std::make_move_iterator(fibers.end())
                    );

            //this->data.surfaceMeshes.push_back(std::move(mesh));
        }
        else
        {
            this->data.surfaceMeshes.reserve(controlPointsInternal.size() + 1);

            for (int i = 0 ; i < controlPointsInternal.size() ; i++)
            {
                this->data.surfaceMeshes.emplace_back(fiber::computeFiberSurfaceSingularSegment(data.tetMesh, data.singularArrangement, data.reebSpace2, {controlPointsInternal[i], controlPointsInternal[(i+1) % controlPointsInternal.size()]}, desiredSheetId));

                const std::vector<FiberPoint> fibers = fiber::computeFiberPointsFromSurfaceMesh(this->data.surfaceMeshes.back(), {});
                //const std::vector<FiberPoint> fibers = fiber::computeFiberPointsFromSurfaceMesh(this->data.surfaceMeshes.back(), this->sibling->selectedSheetIds);
                fibersAll.insert(
                        fibersAll.end(), 
                        std::make_move_iterator(fibers.begin()), 
                        std::make_move_iterator(fibers.end())
                        );

            }
        }





        // Time to beat - 7s
        //std::vector<FiberPoint> fibersAll;
        //if (controlPointsTransformed.size() == 2)
        //{
        //const std::vector<FiberPoint> fibers = fiber::computeFiberSurface(data.tetMesh, data.singularArrangement, data.reebSpace2, {controlPointsInternal[0], controlPointsInternal[1]}, desiredSheetId);

        //fibersAll.insert(
        //fibersAll.end(), 
        //std::make_move_iterator(fibers.begin()), 
        //std::make_move_iterator(fibers.end())
        //);
        //}
        //else
        //{
        //for (int i = 0 ; i < controlPointsInternal.size() ; i++)
        //{
        //const std::vector<FiberPoint> fibers = fiber::computeFiberSurface(data.tetMesh, data.singularArrangement, data.reebSpace2, {controlPointsInternal[i], controlPointsInternal[(i+1) % controlPointsInternal.size()]}, desiredSheetId);

        //fibersAll.insert(
        //fibersAll.end(), 
        //std::make_move_iterator(fibers.begin()), 
        //std::make_move_iterator(fibers.end())
        //);

        //}
        //}



        //std::vector<FiberPoint> fibersAll = io::readDataVtp("/home/peter/Projects/data/reeb-space-test-data/nana/trajectories/State_2/fiberSurfaceExample.vtp").getFiberPoints();


        this->recomputeFiberSurface = false;
        sibling->updateFiberSurface(fibersAll);
    }


























    // ----------------------------------------------------------------
    // Feature Drawing
    // ----------------------------------------------------------------

    if (recomputeFiberSurfaceFeature && this->sibling->selectedSheetIds.size() == 1)
    {
        const int desiredSheetId = *this->sibling->selectedSheetIds.begin();

        const std::vector<std::vector<std::array<double, 2>>> sheetPolygons = data.reebSpace2.computeSheetControlPolygons(desiredSheetId);

        // Draw the control polygons 
        //
        QVector<QVector<QPointF>> controlPointsTransformed(sheetPolygons.size());
        for (int i = 0 ; i < sheetPolygons.size() ; i++)
        {
            const auto &sheetPolygon = sheetPolygons[i];
            controlPointsTransformed[i].resize(sheetPolygon.size());

            for (int j = 0 ; j < sheetPolygon.size() ; j++)
            {
                controlPointsTransformed[i][j] = rescalePoint(sheetPolygon[j][0], sheetPolygon[j][1]);
                p.drawEllipse(controlPointsTransformed[i][j], 20, 20);
            }

            p.drawPolygon(QPolygonF(controlPointsTransformed[i]));
        }




        // Compute the fiber surface
        //
        std::vector<FiberPoint> fibersAll;

        for (const auto &sheetPolygon : sheetPolygons)
        {
            std::vector<std::array<double, 2>> controlPointsInternal = sheetPolygon;

            for (int i = 0 ; i < controlPointsInternal.size(); i++)
            {
                this->data.surfaceMeshes.emplace_back(fiber::computeFiberSurfaceSingularSegment(data.tetMesh, data.singularArrangement, data.reebSpace2, {controlPointsInternal[i], controlPointsInternal[(i+1) % controlPointsInternal.size()]}, desiredSheetId));

                const std::vector<FiberPoint> fibers = fiber::computeFiberPointsFromSurfaceMesh(this->data.surfaceMeshes.back(), this->sibling->selectedSheetIds);

                fibersAll.insert(
                        fibersAll.end(), 
                        std::make_move_iterator(fibers.begin()), 
                        std::make_move_iterator(fibers.end())
                        );

                //const std::vector<FiberPoint> fibers = fiber::computeFiberSurface(
                        //data.tetMesh, 
                        //data.singularArrangement, 
                        //data.reebSpace2, 
                        //{controlPointsInternal[i], controlPointsInternal[(i+1) % controlPointsInternal.size()]}, 
                        //desiredSheetId
                        //);

                //fibersAll.insert(
                        //fibersAll.end(), 
                        //std::make_move_iterator(fibers.begin()), 
                        //std::make_move_iterator(fibers.end())
                        //);

            }

        }

        // Update the fiber
        this->recomputeFiberSurfaceFeature = false;
        sibling->updateFiberSurface(fibersAll);
    }















    p.restore();
    drawAxisLabels2(p);
}



QPointF PlotWidget::rescalePoint(const float &u, const GLfloat &v)
{
    const float rescaledU = (resolution / (paddedMaxF - paddedMinF)) * (u - paddedMinF);
    const float rescaledV = (resolution / (paddedMaxG - paddedMinG)) * (v - paddedMinG);
    return QPointF(rescaledU, rescaledV);
}

void PlotWidget::drawAxisLabels2(QPainter& p)
{
    auto font = p.font();
    auto penBlack = QPen(Qt::black);
    penBlack.setWidthF(5.0);
    p.setPen(penBlack);

    font.setPixelSize(70);
    p.setFont(font);

    float boxOffset = 15;


    // X Axis
    //p.drawLine(boxOffset, resolution - boxOffset, resolution - boxOffset + 100, resolution - boxOffset);
    // Y Axis
    //p.drawLine(boxOffset, resolution - boxOffset, boxOffset, boxOffset - 100);


    float fZero = (resolution / (data.tetMesh.maxF - data.tetMesh.minF)) * (0.0 - data.tetMesh.minF);
    float gZero = (resolution / (data.tetMesh.maxG - data.tetMesh.minG)) * (0.0 - data.tetMesh.minG);

    p.drawLine(fZero, -resolution, fZero, resolution);
    p.drawLine(-resolution, gZero, resolution, gZero);


    // x label
    p.drawText(resolution / 2 - 30, resolution - boxOffset + 5, QString::fromStdString(data.tetMesh.longnameF));

    p.translate(boxOffset + 60, resolution / 2 + 20);
    p.rotate(-90);

    // y label
    p.drawText(0, 0, QString::fromStdString(data.tetMesh.longnameG));
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
