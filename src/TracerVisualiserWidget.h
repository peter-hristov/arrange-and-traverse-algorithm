#pragma once

#ifdef __APPLE__
#include <GLUT/glut.h>
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#endif

#include <QOpenGLWidget>

#include "./Data.h"
#include "./FiberPoint.h"
#include "./ArcBall/Ball.h"

class TracerVisualiserWidget : public QOpenGLWidget
{
  Q_OBJECT 

  public:
    Data &data;
    TracerVisualiserWidget(QWidget*, Data&);
    GLfloat scale = 0;

    //bool drawEdges = false;
    //bool drawFaces = false;
    //bool drawVertices = false;

    // Higher is slower zoom out.
    bool drawEdges = 0;
    bool drawFaces = 0;
    bool drawVertices = 0;

    bool enableLighting = false;

    // Depricated
    bool showUIsosurface = false;
    // Depricated
    bool showVIsosurface = false;

    float edgeOpacity = 0.8;
    float faceOpacity = 0.2;
    float vertexOpacity = 0.8;

    int fiberColour = 0;

    //std::vector<std::vector<float>> fiberColours = {
        //{1, 0, 0},
        //{0, 1, 0},
        //{0, 0, 1}
    //};

    //std::vector<std::vector<float>> fiberColours2 = {
        //{0.0, 0.45, 0.7},   // Blue
        //{0.8, 0.4, 0.0},    // Orange
        //{0.0, 0.6, 0.5},    // Teal
        //{0.9, 0.6, 0.0},    // Gold
        //{0.8, 0.6, 0.7},    // Pink/Mauve
        //{0.35, 0.7, 0.9},   // Sky Blue
        //{0.0, 0.6, 0.3},    // Green
        //{0.9, 0.2, 0.2},    // Red
        //{0.5, 0.5, 0.0}    // Olive
    //};

    bool clearFibers = 1;

    GLfloat isovalueMult = -1.0;

    std::vector<FiberPoint> faceFibers;

    void updateFiber(const std::vector<FiberPoint>&);

    int displayListIndex = 0;
    void generateDisplayList();

    int displayListIndexTriangles = 0;
    int displayListIndexTrianglesG = 0;

    QWidget* sibling;

  protected:
    void initializeGL();
    void resizeGL(int w, int h);
    void paintGL();

    // Events
    void mousePressEvent(QMouseEvent* event);
    void mouseReleaseEvent(QMouseEvent* event);
    void mouseMoveEvent(QMouseEvent* event);
    void wheelEvent(QWheelEvent* event);
    void mouseDoubleClickEvent(QMouseEvent*);
    void keyPressEvent(QKeyEvent* event);


  private:

    // Arcball stuff
    QPointF position;
    BallData theBall;

    // Render Triangles
    void drawSolidTriangle(GLfloat vertices[3][3]);

    // Render Various Functions
    void cube();
    void drawAxis(GLfloat, GLfloat);
    void drawWiredCube(const GLfloat vertices[8][3]);

    void drawScene();

    float translateX = 0.;
    float translateY = 0.;

    float initialX = -1.;
    float initialY = -1.;

    // Utility
    void setMaterial(GLfloat, GLfloat, GLfloat, GLfloat, GLfloat);

}; // class GLPolygonWidget
