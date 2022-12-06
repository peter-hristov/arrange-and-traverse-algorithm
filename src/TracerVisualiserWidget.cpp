//#include <qpoint.h>
#ifdef __APPLE__
#include <GLUT/glut.h>
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#endif

#include <vector>

#include "./utility/Geometry.h"
#include "./TracerVisualiserWidget.h"
#include "./TracerVisualiserWindow.h"

using namespace std;

void
TracerVisualiserWidget::setMaterial(GLfloat red, GLfloat green, GLfloat blue, GLfloat alpha, GLfloat shininess)
{
    glMaterialfv(GL_FRONT, GL_AMBIENT, &(vector<GLfloat>({ red, green, blue, alpha })[0]));
    glMaterialfv(GL_FRONT, GL_DIFFUSE, &(vector<GLfloat>({ red, green, blue, alpha })[0]));
    glMaterialfv(GL_FRONT, GL_SPECULAR, &(vector<GLfloat>({ red, green, blue, alpha })[0]));
    glMaterialf(GL_FRONT, GL_SHININESS, shininess);
}

TracerVisualiserWidget::TracerVisualiserWidget(QWidget* parent,
                                               QWidget* _sibling,
                                               Data* _data)
  : QOpenGLWidget(parent)
{
    // Set points to other singletons
    this->data = _data;
    this->sibling = _sibling;

    // Default values for paraters
    this->scale = (data->maxZ - data->minZ) * 2;

    // Initialise Arcball
    Ball_Init(&theBall);
    Ball_Place(&theBall, qOne, 1);
}


void
TracerVisualiserWidget::initializeGL()
{
    glClearColor(0.3, 0.3, 0.3, 0.0);

    //glEnable(GL_LIGHTING);
    //glEnable(GL_LIGHT0);

    // Enable alpha blending
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    // Set the projection matrix
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(90, 1, 0.01, 10000);
}

void
TracerVisualiserWidget::resizeGL(int w, int h)
{
    glViewport(0, 0, w, h);
}


// Draw fibers
void
TracerVisualiserWidget::generateDisplayList()
{
    glDeleteLists(displayListIndex, 1);
    displayListIndex = glGenLists(1);
    glNewList(displayListIndex, GL_COMPILE);

    //setMaterial(1, 0, 0, 1.0, 0.0);

    // Draw Fiber
    glBegin(GL_LINES);
    {
        for(const auto &faceFiber : this->data->faceFibers)
        {
            glColor3fv(faceFiber.colour.data());
            glVertex3fv(faceFiber.point.data());
        }
    }
    glEnd();

    // Draw fiber endpoints (in every tet)
    for(const auto &faceFiber : this->data->faceFibers)
    {
        glColor3fv(faceFiber.colour.data());

        glPushMatrix();
        {
            glTranslatef(faceFiber.point[0], faceFiber.point[1], faceFiber.point[2]);
            GLUquadric* sphere = gluNewQuadric();
            gluSphere(sphere, 0.01, 10, 10);
            delete sphere;
        }
        glPopMatrix();

    }

    glEndList();
}

void
TracerVisualiserWidget::drawSolidTriangle(GLfloat vertices[3][3])
{
    GLfloat a[] = { vertices[1][0] - vertices[0][0], vertices[1][1] - vertices[0][1], vertices[1][2] - vertices[0][2] };

    GLfloat b[] = { vertices[2][0] - vertices[0][0], vertices[2][1] - vertices[0][1], vertices[2][2] - vertices[0][2] };

    // Invert Normal for sub/super level sets
    GLfloat normal[] = { isovalueMult * (a[2] * b[1] - b[1] * b[2]),
                         isovalueMult * (a[0] * b[2] - a[2] * b[0]),
                         isovalueMult * (a[1] * b[0] - a[0] * b[1]) };

    GLfloat normalFlipped[] = { -1 * normal[0], -1 * normal[1], -1 * normal[2] };

    glNormal3fv(normalFlipped);
    glVertex3fv(vertices[0]);
    glVertex3fv(vertices[1]);
    glVertex3fv(vertices[2]);
}

void
TracerVisualiserWidget::drawAxis(GLfloat length, GLfloat width)
{
    GLUquadric* line = gluNewQuadric();

    // Z Axis
    setMaterial(0., 0., 1., 1.0, 30.);
    glColor3f(0, 0, 1);
    gluCylinder(line, width, width, length, 100, 100);

    // Y Axis
    setMaterial(0., 1., 0., 1.0, 30.);
    glColor3f(0, 1, 0);
    glPushMatrix();
    glRotatef(-90, 1, 0, 0);
    gluCylinder(line, width, width, length, 100, 100);
    glPopMatrix();

    // X Axis
    setMaterial(1., 0., 0., 1.0, 30.);
    glColor3f(1, 0, 0);
    glPushMatrix();
    glRotatef(90, 0, 1, 0);
    gluCylinder(line, width, width, length, 100, 100);
    glPopMatrix();

    gluDeleteQuadric(line);
}

void
TracerVisualiserWidget::drawWiredCube(const GLfloat vertices[8][3])
{
    glBegin(GL_LINES);
    {
        // Front Side
        glVertex3fv(vertices[0]);
        glVertex3fv(vertices[2]);

        glVertex3fv(vertices[2]);
        glVertex3fv(vertices[3]);

        glVertex3fv(vertices[1]);
        glVertex3fv(vertices[3]);

        glVertex3fv(vertices[0]);
        glVertex3fv(vertices[1]);

        // Back Side
        glVertex3fv(vertices[4]);
        glVertex3fv(vertices[6]);

        glVertex3fv(vertices[6]);
        glVertex3fv(vertices[7]);

        glVertex3fv(vertices[5]);
        glVertex3fv(vertices[7]);

        glVertex3fv(vertices[4]);
        glVertex3fv(vertices[5]);

        // Sides
        glVertex3fv(vertices[0]);
        glVertex3fv(vertices[4]);

        glVertex3fv(vertices[2]);
        glVertex3fv(vertices[6]);

        glVertex3fv(vertices[3]);
        glVertex3fv(vertices[7]);

        glVertex3fv(vertices[1]);
        glVertex3fv(vertices[5]);
    }
    glEnd();
}

void
TracerVisualiserWidget::drawScene()
{
    glEnable(GL_DEPTH_TEST);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // Rescale normals when scaling
    glEnable(GL_RESCALE_NORMAL);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    GLfloat light_pos[] = { 0, 0, 10000, 1. };
    glLightfv(GL_LIGHT0, GL_POSITION, light_pos);
    glLightf(GL_LIGHT0, GL_SPOT_CUTOFF, 180.);

    // Scrolling
    glTranslatef(0.0, 0.0, -1 * scale / 2.6);

    // Offset along x, y (with the right mouse button)
    glTranslatef(translateX, translateY, 0.0);

    // Arcball
    GLfloat mNow[16];
    Ball_Value(&theBall, mNow);
    glMultMatrixf(mNow);

    // GLUquadric* sphere = gluNewQuadric();
    // gluSphere(sphere, 2, 3, 3);

    // Data Center
    glTranslatef((-1.0 * (this->data->xDim - 1)) / 2.0, 0, 0);
    glTranslatef(0, 0, ((this->data->yDim) - 1) / 2.0);
    glTranslatef(0, -1.0 * (this->data->zDim - 1) / 2.0, 0);
    //glTranslatef(-1.0 * (this->data->maxX - this->data->minX) / 2, 0, 0);
    //glTranslatef(0, 0, (this->data->maxY - this->data->minY) / 2);
    //glTranslatef(0, -1.0 * (this->data->maxZ - this->data->minZ) / 2, 0);

    // Data center
    //glTranslatef(-1.0 * this->data->vertexDomainCoordinates[this->data->vertexDomainCoordinates.size() - 1][0] , 0, 0);
    //glTranslatef(0, 0, this->data->vertexDomainCoordinates[this->data->vertexDomainCoordinates.size() - 1][1]);
    //glTranslatef(0, -1.0 * this->data->vertexDomainCoordinates[this->data->vertexDomainCoordinates.size() - 1][2], 0);

    glRotatef(-90., 1., 0., 0.);

    this->drawAxis(1000., 1.0 * (this->data->maxX - this->data->minX) / 1800.0);

    glColor3f(1, 1, 1);

    // Draw bounding cube
    glPushMatrix();
    {
        glScalef(this->data->xDim - 1, this->data->yDim - 1, this->data->zDim - 1);
        setMaterial(255, 255, 255, 100, 30.0);
        this->drawWiredCube(tv9k::geometry::cubeVertices);
    }
    glPopMatrix();

    if (true == this->drawEdges)
    {
        glColor4f(1, 1, 1, this->edgeOpacity);

        // Tet Edges
        glBegin(GL_LINES);
        {
            for(int t = 0 ; t < this->data->tetrahedra.size() ; t++)
            {
                //if (this->data->tetsWithFibers[t] == false) {continue;}
                const auto tet = this->data->tetrahedra[t];

                for(int i = 0 ; i < 4 ; i++)
                {
                    for(int j = i + 1 ; j < 4 ; j++)
                    {
                        GLfloat pointA[3], pointB[3];

                        pointA[0] = this->data->vertexDomainCoordinates[tet[i]][0];
                        pointA[1] = this->data->vertexDomainCoordinates[tet[i]][1];
                        pointA[2] = this->data->vertexDomainCoordinates[tet[i]][2];

                        pointB[0] = this->data->vertexDomainCoordinates[tet[j]][0];
                        pointB[1] = this->data->vertexDomainCoordinates[tet[j]][1];
                        pointB[2] = this->data->vertexDomainCoordinates[tet[j]][2];

                        glVertex3fv(pointA);
                        glVertex3fv(pointB);
                    }
                }
                //cout << endl;
            }
        }
        glEnd();
    }


    glPushMatrix();
    {

        glDeleteLists(displayListIndex, 1);
        glCallList(displayListIndex);
    }
    glPopMatrix();


    if (this->showVIsosurface == true)
    {
        glPushMatrix();
        {
            glCallList(displayListIndexTrianglesG);
        }
        glPopMatrix();
    }

    if (this->showUIsosurface == true)
    {
        glPushMatrix();
        {
            glCallList(displayListIndexTriangles);
        }
        glPopMatrix();
    }

    if (true == this->drawVertices)
    {
        // Draw Vertices
        {
            //for (const auto &vertex : this->data->vertexDomainCoordinates) 
            for (int i = 0 ; i < this->data->vertexDomainCoordinates.size() ; i++) 
            {
                const auto &vertex = this->data->vertexDomainCoordinates[i];

                glColor4f(1, 1, 1, this->vertexOpacity);
                if (dynamic_cast<PlotWidget*>(sibling)->dragMode ==  PlotWidget::MouseDragMode::DataPoint && i == dynamic_cast<PlotWidget*>(sibling)->movePoint)
                {
                    glColor4f(1, 0, 0, 0.2);
                }
                glPushMatrix();
                {
                    glTranslatef(vertex[0], vertex[1], vertex[2]);
                    GLUquadric* sphere = gluNewQuadric();
                    gluSphere(sphere, 0.01, 10, 10);
                    delete sphere;
                }
                glPopMatrix();
            }

        }
    }

    // Tet Faces
    if (true == drawFaces)
    {
        glColor4f(1, 1, 1, this->faceOpacity);

        int centerVertexId = this->data->vertexDomainCoordinates.size() - 1;

        int triangles = 0;

        glBegin(GL_TRIANGLES);
        {
            //for(const auto &tet : this->data->tetrahedra)
            for(int t = 0 ; t < this->data->tetrahedra.size() ; t++)
            {
                //if (this->data->tetsWithFibers[t] == false) {continue;}
                const auto tet = this->data->tetrahedra[t];

                bool isCorrectTet = true;
                //bool isCorrectTet = false;

                for(int i = 0 ; i < 4 ; i++)
                {
                    if (tet[i] == centerVertexId)
                    {
                        isCorrectTet = true;
                    }
                }

                if (false == isCorrectTet)
                {
                    continue;
                }

                for(int i = 0 ; i < 4 ; i++)
                {
                    for(int j = i + 1 ; j < 4 ; j++)
                    {
                        for(int k = j + 1 ; k < 4 ; k++)
                        {
                            if (tet[i] == centerVertexId || tet[j] == centerVertexId || tet[k] == centerVertexId)
                            //if (!(tet[i] == centerVertexId || tet[j] == centerVertexId || tet[k] == centerVertexId))
                            {
                                continue;
                            }

                            GLfloat pointA[3], pointB[3], pointC[3];

                            pointA[0] = this->data->vertexDomainCoordinates[tet[i]][0];
                            pointA[1] = this->data->vertexDomainCoordinates[tet[i]][1];
                            pointA[2] = this->data->vertexDomainCoordinates[tet[i]][2];

                            pointB[0] = this->data->vertexDomainCoordinates[tet[j]][0];
                            pointB[1] = this->data->vertexDomainCoordinates[tet[j]][1];
                            pointB[2] = this->data->vertexDomainCoordinates[tet[j]][2];

                            pointC[0] = this->data->vertexDomainCoordinates[tet[k]][0];
                            pointC[1] = this->data->vertexDomainCoordinates[tet[k]][1];
                            pointC[2] = this->data->vertexDomainCoordinates[tet[k]][2];

                            glVertex3fv(pointA);
                            glVertex3fv(pointB);
                            glVertex3fv(pointC);

                            triangles++;
                        }
                    }
                }
                //cout << endl;
            }
        }
    }
    glEnd();

    glFlush();
}

void
TracerVisualiserWidget::paintGL()
{
    this->drawScene();
}

void
TracerVisualiserWidget::wheelEvent(QWheelEvent* event)
{
    this->scale -= event->angleDelta().y() / 100.0;

    if (scale <= 0) {
        scale = 0;
    }

    this->update();
}

void
TracerVisualiserWidget::mousePressEvent(QMouseEvent* event)
{
    if (event->button() == Qt::LeftButton) {
        HVect vNow;
        vNow.x = (2.0 * event->x() - width()) / width();
        vNow.y = (height() - 2.0 * event->y()) / height();

        Ball_Mouse(&theBall, vNow);
        Ball_BeginDrag(&theBall);

        this->update();
    }
    if (event->button() == Qt::RightButton) {
        initialX = event->localPos().x();
        initialY = event->localPos().y();
    }
}

void
TracerVisualiserWidget::mouseMoveEvent(QMouseEvent* event)
{
    if (event->buttons() == Qt::LeftButton) {
        HVect vNow;
        vNow.x = (2.0 * event->localPos().x() - width()) / width();
        vNow.y = (height() - 2.0 * event->y()) / height();

        Ball_Mouse(&theBall, vNow);
        Ball_Update(&theBall);

        this->update();
    } else if (event->buttons() == Qt::RightButton) {
        float x = event->localPos().x();
        float y = event->localPos().y();

        translateX -= (initialX - x) / 10;
        translateY += (initialY - y) / 10;

        initialX = event->localPos().x();
        initialY = event->localPos().y();

        this->update();
    }
}

void
TracerVisualiserWidget::mouseReleaseEvent(QMouseEvent* event)
{
    if (event->button() == Qt::LeftButton) {
        Ball_EndDrag(&theBall);
        this->update();
    }
}

// Depricated
void
TracerVisualiserWidget::keyPressEvent(QKeyEvent* event)
{
    if (event->key() == Qt::Key_W) {
        //translateX += 0.1 * this->data->xdim;
        this->update();
    }
    if (event->key() == Qt::Key_A) {
        //translateY += 0.1 * this->data->xdim;
        this->update();
    }
    if (event->key() == Qt::Key_S) {
        //translateX -= 0.1 * this->data->xdim;
        this->update();
    }
    if (event->key() == Qt::Key_D) {
        //translateY -= 0.1 * this->data->xdim;
        this->update();
    }
}

void
TracerVisualiserWidget::mouseDoubleClickEvent(QMouseEvent* event)
{
    //this->data->faceFibers.clear();
    this->generateDisplayList();
    this->update();
}

void
TracerVisualiserWidget::generateDisplayListTriangles(const float isovalue, const vector<float> &values, const int type)
{
    this->data->meshTriangles.clear();

    // For every tet do marching tets
    for(size_t tetId = 0 ; tetId < this->data->tetrahedra.size(); tetId++)
    {
        const auto tet = this->data->tetrahedra[tetId];

        vector<size_t> smallerVertices;
        vector<size_t> biggerOrEqualVertices;

        // For every vertex in the tet
        for(int i = 0 ; i < 4 ; i++)
        {
            const float vertexValue = values[tet[i]];

            if (vertexValue < isovalue)
            {
                smallerVertices.push_back(tet[i]);
            }
            else
            {
                biggerOrEqualVertices.push_back(tet[i]);
            }
        }

        if (smallerVertices.size() == 0)
        {
        }
        else if (smallerVertices.size() == 1)
        {

            Data::MeshTriangle triangle;
            // s[0] -> b[0]
            float t0 = 
                (isovalue - values[smallerVertices[0]]) / 
                (values[biggerOrEqualVertices[0]] - values[smallerVertices[0]]);

            float t1 = 
                (isovalue - values[smallerVertices[0]]) / 
                (values[biggerOrEqualVertices[1]] - values[smallerVertices[0]]);

            float t2 = 
                (isovalue - values[smallerVertices[0]]) / 
                (values[biggerOrEqualVertices[2]] - values[smallerVertices[0]]);

            triangle.vertixA = {
                (1 - t0) * this->data->vertexDomainCoordinates[smallerVertices[0]][0] + t0 * this->data->vertexDomainCoordinates[biggerOrEqualVertices[0]][0],
                (1 - t0) * this->data->vertexDomainCoordinates[smallerVertices[0]][1] + t0 * this->data->vertexDomainCoordinates[biggerOrEqualVertices[0]][1],
                (1 - t0) * this->data->vertexDomainCoordinates[smallerVertices[0]][2] + t0 * this->data->vertexDomainCoordinates[biggerOrEqualVertices[0]][2],
            };

            triangle.vertixB = {
                (1 - t1) * this->data->vertexDomainCoordinates[smallerVertices[0]][0] + t1 * this->data->vertexDomainCoordinates[biggerOrEqualVertices[1]][0],
                (1 - t1) * this->data->vertexDomainCoordinates[smallerVertices[0]][1] + t1 * this->data->vertexDomainCoordinates[biggerOrEqualVertices[1]][1],
                (1 - t1) * this->data->vertexDomainCoordinates[smallerVertices[0]][2] + t1 * this->data->vertexDomainCoordinates[biggerOrEqualVertices[1]][2],
            };

            triangle.vertixC = {
                (1 - t2) * this->data->vertexDomainCoordinates[smallerVertices[0]][0] + t2 * this->data->vertexDomainCoordinates[biggerOrEqualVertices[2]][0],
                (1 - t2) * this->data->vertexDomainCoordinates[smallerVertices[0]][1] + t2 * this->data->vertexDomainCoordinates[biggerOrEqualVertices[2]][1],
                (1 - t2) * this->data->vertexDomainCoordinates[smallerVertices[0]][2] + t2 * this->data->vertexDomainCoordinates[biggerOrEqualVertices[2]][2],
            };

            this->data->meshTriangles.push_back(triangle);
        }
        else if (smallerVertices.size() == 2)
        {
            float t0 = 
                (isovalue - values[smallerVertices[0]]) / 
                (values[biggerOrEqualVertices[0]] - values[smallerVertices[0]]);

            float t1 = 
                (isovalue - values[smallerVertices[0]]) / 
                (values[biggerOrEqualVertices[1]] - values[smallerVertices[0]]);

            float t2 = 
                (isovalue - values[smallerVertices[1]]) / 
                (values[biggerOrEqualVertices[0]] - values[smallerVertices[1]]);

            float t3 = 
                (isovalue - values[smallerVertices[1]]) / 
                (values[biggerOrEqualVertices[1]] - values[smallerVertices[1]]);



            Data::MeshTriangle triangle1;

            triangle1.vertixA = {
                (1 - t0) * this->data->vertexDomainCoordinates[smallerVertices[0]][0] + t0 * this->data->vertexDomainCoordinates[biggerOrEqualVertices[0]][0],
                (1 - t0) * this->data->vertexDomainCoordinates[smallerVertices[0]][1] + t0 * this->data->vertexDomainCoordinates[biggerOrEqualVertices[0]][1],
                (1 - t0) * this->data->vertexDomainCoordinates[smallerVertices[0]][2] + t0 * this->data->vertexDomainCoordinates[biggerOrEqualVertices[0]][2],
            };

            triangle1.vertixB = {
                (1 - t1) * this->data->vertexDomainCoordinates[smallerVertices[0]][0] + t1 * this->data->vertexDomainCoordinates[biggerOrEqualVertices[1]][0],
                (1 - t1) * this->data->vertexDomainCoordinates[smallerVertices[0]][1] + t1 * this->data->vertexDomainCoordinates[biggerOrEqualVertices[1]][1],
                (1 - t1) * this->data->vertexDomainCoordinates[smallerVertices[0]][2] + t1 * this->data->vertexDomainCoordinates[biggerOrEqualVertices[1]][2],
            };

            triangle1.vertixC = {
                (1 - t2) * this->data->vertexDomainCoordinates[smallerVertices[1]][0] + t2 * this->data->vertexDomainCoordinates[biggerOrEqualVertices[0]][0],
                (1 - t2) * this->data->vertexDomainCoordinates[smallerVertices[1]][1] + t2 * this->data->vertexDomainCoordinates[biggerOrEqualVertices[0]][1],
                (1 - t2) * this->data->vertexDomainCoordinates[smallerVertices[1]][2] + t2 * this->data->vertexDomainCoordinates[biggerOrEqualVertices[0]][2],
            };

            Data::MeshTriangle triangle2;

            triangle2.vertixA = {
                (1 - t2) * this->data->vertexDomainCoordinates[smallerVertices[1]][0] + t2 * this->data->vertexDomainCoordinates[biggerOrEqualVertices[0]][0],
                (1 - t2) * this->data->vertexDomainCoordinates[smallerVertices[1]][1] + t2 * this->data->vertexDomainCoordinates[biggerOrEqualVertices[0]][1],
                (1 - t2) * this->data->vertexDomainCoordinates[smallerVertices[1]][2] + t2 * this->data->vertexDomainCoordinates[biggerOrEqualVertices[0]][2],
            };

            triangle2.vertixB = {
                (1 - t3) * this->data->vertexDomainCoordinates[smallerVertices[1]][0] + t3 * this->data->vertexDomainCoordinates[biggerOrEqualVertices[1]][0],
                (1 - t3) * this->data->vertexDomainCoordinates[smallerVertices[1]][1] + t3 * this->data->vertexDomainCoordinates[biggerOrEqualVertices[1]][1],
                (1 - t3) * this->data->vertexDomainCoordinates[smallerVertices[1]][2] + t3 * this->data->vertexDomainCoordinates[biggerOrEqualVertices[1]][2],
            };

            triangle2.vertixC = {
                (1 - t1) * this->data->vertexDomainCoordinates[smallerVertices[0]][0] + t1 * this->data->vertexDomainCoordinates[biggerOrEqualVertices[1]][0],
                (1 - t1) * this->data->vertexDomainCoordinates[smallerVertices[0]][1] + t1 * this->data->vertexDomainCoordinates[biggerOrEqualVertices[1]][1],
                (1 - t1) * this->data->vertexDomainCoordinates[smallerVertices[0]][2] + t1 * this->data->vertexDomainCoordinates[biggerOrEqualVertices[1]][2],
            };

            this->data->meshTriangles.push_back(triangle1);
            this->data->meshTriangles.push_back(triangle2);
        }
        else if (smallerVertices.size() == 3)
        {
            Data::MeshTriangle triangle;
            // Interpolate along all edges
            // s[0] -> b[0]
            float t0 = 
                (isovalue - values[biggerOrEqualVertices[0]]) / 
                (values[smallerVertices[0]] - values[biggerOrEqualVertices[0]]);

            float t1 = 
                (isovalue - values[biggerOrEqualVertices[0]]) / 
                (values[smallerVertices[1]] - values[biggerOrEqualVertices[0]]);

            float t2 = 
                (isovalue - values[biggerOrEqualVertices[0]]) / 
                (values[smallerVertices[2]] - values[biggerOrEqualVertices[0]]);

            triangle.vertixA = {
                (1 - t0) * this->data->vertexDomainCoordinates[biggerOrEqualVertices[0]][0] + t0 * this->data->vertexDomainCoordinates[smallerVertices[0]][0],
                (1 - t0) * this->data->vertexDomainCoordinates[biggerOrEqualVertices[0]][1] + t0 * this->data->vertexDomainCoordinates[smallerVertices[0]][1],
                (1 - t0) * this->data->vertexDomainCoordinates[biggerOrEqualVertices[0]][2] + t0 * this->data->vertexDomainCoordinates[smallerVertices[0]][2],
            };

            triangle.vertixB = {
                (1 - t1) * this->data->vertexDomainCoordinates[biggerOrEqualVertices[0]][0] + t1 * this->data->vertexDomainCoordinates[smallerVertices[1]][0],
                (1 - t1) * this->data->vertexDomainCoordinates[biggerOrEqualVertices[0]][1] + t1 * this->data->vertexDomainCoordinates[smallerVertices[1]][1],
                (1 - t1) * this->data->vertexDomainCoordinates[biggerOrEqualVertices[0]][2] + t1 * this->data->vertexDomainCoordinates[smallerVertices[1]][2],
            };

            triangle.vertixC = {
                (1 - t2) * this->data->vertexDomainCoordinates[biggerOrEqualVertices[0]][0] + t2 * this->data->vertexDomainCoordinates[smallerVertices[2]][0],
                (1 - t2) * this->data->vertexDomainCoordinates[biggerOrEqualVertices[0]][1] + t2 * this->data->vertexDomainCoordinates[smallerVertices[2]][1],
                (1 - t2) * this->data->vertexDomainCoordinates[biggerOrEqualVertices[0]][2] + t2 * this->data->vertexDomainCoordinates[smallerVertices[2]][2],
            };

            this->data->meshTriangles.push_back(triangle);
        }
    }

    if (type == 1)
    {
        glDeleteLists(displayListIndexTriangles, 1);
        displayListIndexTriangles = glGenLists(1);
        glNewList(displayListIndexTriangles, GL_COMPILE);
    }
    else{
        glDeleteLists(displayListIndexTrianglesG, 1);
        displayListIndexTrianglesG = glGenLists(1);
        glNewList(displayListIndexTrianglesG, GL_COMPILE);
    }

    //setMaterial(1, 0, 0, 1.0, 0.0);
    if (type == 1)
    {
        glColor4f(0, 1, 0, 0.5);
    }
    else
    {
        glColor4f(0, 0, 1, 0.5);
    }

    // Draw Fiber
    glBegin(GL_TRIANGLES);
    {
        for (const auto &point : this->data->meshTriangles)
        {
            GLfloat pointA[3] = { point.vertixA[0], point.vertixA[1], point.vertixA[2] };
            GLfloat pointB[3] = { point.vertixB[0], point.vertixB[1], point.vertixB[2] };
            GLfloat pointC[3] = { point.vertixC[0], point.vertixC[1], point.vertixC[2] };

            //GLfloat pointA[3] = { 0, 0, 0 };
            //GLfloat pointB[3] = { 0, 10, 0 };
            //GLfloat pointC[3] = { 0, 0, 10 };

            //printf("Triangle : (%f, %f, %f) | (%f, %f, %f) | (%f, %f, %f).\n", pointA[0], pointA[1], pointA[2], pointB[0], pointB[1], pointB[2], pointC[0], pointC[1], pointC[2]);

            glVertex3fv(pointA);
            glVertex3fv(pointB);
            glVertex3fv(pointC);
        }
    }
    glEnd();

    glEndList();
}
