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
#include "./ArcBall/Ball.h"

class TracerVisualiserWidget : public QOpenGLWidget
{
  Q_OBJECT 

  public:
    TracerVisualiserWidget(QWidget*, QWidget*, Data*);
    GLfloat scale = 1;

    //bool drawEdges = false;
    //bool drawFaces = false;
    //bool drawVertices = false;

    bool drawEdges = 0;
    bool drawFaces = 0;
    bool drawVertices = 0;

    bool showUIsosurface = true;
    bool showVIsosurface = true;

    float edgeOpacity = 0.8;
    float faceOpacity = 0.2;
    float vertexOpacity = 0.8;

    int fiberColour = 0;
    std::vector<std::vector<float>> fiberColours = {
        {1, 0, 0},
        {0, 1, 0},
        {0, 0, 1}
    };

    bool clearFibers = 1;

    GLfloat isovalueMult = -1.0;

    int displayListIndex = 0;
    void generateDisplayList();

    int displayListIndexTriangles = 0;
    int displayListIndexTrianglesG = 0;
    void generateDisplayListTriangles(const float isovalue, const std::vector<float> &values, const int type);

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
    Data* data;
    QWidget* sibling;

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
