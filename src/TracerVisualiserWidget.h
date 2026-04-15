#pragma once

#include "./CGALTypedefs.h"

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

    bool drawFibers = 1;
    bool drawFiberSurfaces = 1;
    bool drawFeatureSurfaces = 1;

    bool enableLighting = 1;

    // Depricated
    bool showUIsosurface = false;
    // Depricated
    bool showVIsosurface = false;

    float fiberOpcity = 1.8;
    float fsOpacity = 1.0;
    float featureOpacity = 1.8;

    int fiberColour = 0;

    std::set<int> selectedSheetIds;

    bool traceFibers = true;

    GLfloat isovalueMult = -1.0;

    void updateFiber();
    void updateFiberSurface();
    void updateFiberSurfaceFeatures();

    void clearFibers();
    void clearFiberSurface();
    void clearFiberSurfaceSheets();

    void renderSurface(std::vector<FiberSurface> &surfaceMesh);


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
    std::array<GLfloat, 3> computeTriangleNormal(const GLfloat* v0, const GLfloat* v1, const GLfloat* v2);

    void renderMolecule();

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


    // Store these as members
    std::vector<CartesianTriangle_3> pickingTriangles;
    std::vector<int>                 pickingSheetIds;

    std::vector<TriangleTree> aabbTriangleTrees;

    void buildAABBTree();

    int pickSegment(int mouseX, int mouseY);
}; // class GLPolygonWidget
