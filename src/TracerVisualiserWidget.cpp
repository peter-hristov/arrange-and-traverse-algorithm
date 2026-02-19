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

TracerVisualiserWidget::TracerVisualiserWidget(QWidget* parent, Data &_data)
  : QOpenGLWidget(parent)
  , data(_data)
{
    // Default values for paraters
    this->scale = (data.tetMesh.maxZ - data.tetMesh.minZ) * 2;

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
    //glEnable(GL_COLOR_MATERIAL);
    //glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);

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
    //glBegin(GL_LINES);
    //{
        //for(const auto &faceFiber : this->faceFibers)
        //{
            //glColor3fv(faceFiber.colour.data());
            //glVertex3fv(faceFiber.point.data());
        //}
    //}
    //glEnd();


    glBegin(GL_TRIANGLES);
    {
        for(const auto &faceFiber : this->faceFibers)
        {
            glColor3fv(faceFiber.colour.data());
            glVertex3fv(faceFiber.point.data());
        }
    }
    glEnd();





    // Draw fiber endpoints (in every tet)
    //for(const auto &faceFiber : this->data.faceFibers)
    //{
        //glColor3fv(faceFiber.colour.data());

        //glPushMatrix();
        //{
            //glTranslatef(faceFiber.point[0], faceFiber.point[1], faceFiber.point[2]);
            //GLUquadric* sphere = gluNewQuadric();
            //gluSphere(sphere, 0.01, 10, 10);
            //delete sphere;
        //}
        //glPopMatrix();

    //}

    glEndList();
}

void TracerVisualiserWidget::drawSolidTriangle(GLfloat vertices[3][3])
{
    // Compute edges
    GLfloat a[3] = {
        vertices[1][0] - vertices[0][0],
        vertices[1][1] - vertices[0][1],
        vertices[1][2] - vertices[0][2]
    };

    GLfloat b[3] = {
        vertices[2][0] - vertices[0][0],
        vertices[2][1] - vertices[0][1],
        vertices[2][2] - vertices[0][2]
    };

    // Correct cross product for the normal
    GLfloat normal[3] = {
        a[1]*b[2] - a[2]*b[1],
        a[2]*b[0] - a[0]*b[2],
        a[0]*b[1] - a[1]*b[0]
    };

    // Optionally scale by isovalueMult
    //normal[0] *= isovalueMult;
    //normal[1] *= isovalueMult;
    //normal[2] *= isovalueMult;

    // Normalize the normal
    GLfloat len = std::sqrt(normal[0]*normal[0] +
                            normal[1]*normal[1] +
                            normal[2]*normal[2]);
    if (len != 0.0f) {
        normal[0] /= len;
        normal[1] /= len;
        normal[2] /= len;
    }

    // Flip normal if needed
    GLfloat normalFlipped[3] = { -normal[0], -normal[1], -normal[2] };

    // Set normal for OpenGL
    glNormal3fv(normalFlipped);

    // Draw the triangle
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
    glTranslatef(-1.0 * (this->data.tetMesh.maxX + this->data.tetMesh.minX) / 2.0, 0, 0);
    glTranslatef(0, 0, (this->data.tetMesh.maxY + this->data.tetMesh.minY) / 2.0);
    glTranslatef(0, -1.0 * (this->data.tetMesh.maxZ + this->data.tetMesh.minZ) / 2.0, 0);

    // Don't remember why I need this
    glRotatef(-90., 1., 0., 0.);

    //this->drawAxis(1000., 1.0 * (this->data.tetMesh.maxX - this->data.tetMesh.minX) / 1800.0);


    //
    // Data bounding box
    //
    GLfloat vertices[8][3] = {
        {this->data.tetMesh.minX, this->data.tetMesh.minY, this->data.tetMesh.minZ},
        {this->data.tetMesh.minX, this->data.tetMesh.minY, this->data.tetMesh.maxZ},
        {this->data.tetMesh.minX, this->data.tetMesh.maxY, this->data.tetMesh.minZ},
        {this->data.tetMesh.minX, this->data.tetMesh.maxY, this->data.tetMesh.maxZ},

        {this->data.tetMesh.maxX, this->data.tetMesh.minY, this->data.tetMesh.minZ},
        {this->data.tetMesh.maxX, this->data.tetMesh.minY, this->data.tetMesh.maxZ},
        {this->data.tetMesh.maxX, this->data.tetMesh.maxY, this->data.tetMesh.minZ},
        {this->data.tetMesh.maxX, this->data.tetMesh.maxY, this->data.tetMesh.maxZ},
    };

    setMaterial(0., 0., 1., 1.0, 30.);
    glColor3f(1.0, 1.0, 1.0);
    drawWiredCube(vertices);


    glColor3f(1, 1, 1);

    if (true == this->drawEdges)
    {
        glColor4f(1, 1, 1, this->edgeOpacity);

        // Tet Edges
        glBegin(GL_LINES);
        {
            for(int t = 0 ; t < this->data.tetMesh.tetrahedra.size() ; t++)
            {
                //if (this->data.tetsWithFibers[t] == false) {continue;}
                const auto tet = this->data.tetMesh.tetrahedra[t];

                for(int i = 0 ; i < 4 ; i++)
                {
                    for(int j = i + 1 ; j < 4 ; j++)
                    {
                        // Get the indices of the vertices for the edge
                        int aIndex = tet[i];
                        int bIndex = tet[j];

                        // Make sure the vertices of the edge are in sorted order to have consistent orientation
                        if (aIndex > bIndex)
                        {
                            std::swap(aIndex, bIndex);
                        }

                        const int edgeType = data.tetMesh.edgeSingularTypes.at({aIndex, bIndex});

                        if (edgeType == 0)
                        {
                            glColor3f(0.6f, 0.85f, 0.85f);
                        }
                        else if(edgeType == 2)
                        {
                            glColor3f(1.0f, 0.647f, 0.0f);
                        }

                        if (edgeType !=1)
                        {
                            GLfloat pointA[3], pointB[3];

                            pointA[0] = this->data.tetMesh.vertexDomainCoordinates[tet[i]][0];
                            pointA[1] = this->data.tetMesh.vertexDomainCoordinates[tet[i]][1];
                            pointA[2] = this->data.tetMesh.vertexDomainCoordinates[tet[i]][2];

                            pointB[0] = this->data.tetMesh.vertexDomainCoordinates[tet[j]][0];
                            pointB[1] = this->data.tetMesh.vertexDomainCoordinates[tet[j]][1];
                            pointB[2] = this->data.tetMesh.vertexDomainCoordinates[tet[j]][2];

                            glVertex3fv(pointA);
                            glVertex3fv(pointB);
                        }
                    }
                }
            }
        }
        glEnd();
    }

    glPushMatrix();
    {
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
            for (int i = 0 ; i < this->data.tetMesh.vertexDomainCoordinates.size() ; i++) 
            {
                const auto &vertex = this->data.tetMesh.vertexDomainCoordinates[i];

                glColor4f(1, 1, 1, this->vertexOpacity);
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

        int centerVertexId = this->data.tetMesh.vertexDomainCoordinates.size() - 1;

        int triangles = 0;

        glBegin(GL_TRIANGLES);
        {
            //for(const auto &tet : this->data.tetrahedra)
            for(int t = 0 ; t < this->data.tetMesh.tetrahedra.size() ; t++)
            {
                //if (this->data.tetsWithFibers[t] == false) {continue;}
                const auto tet = this->data.tetMesh.tetrahedra[t];

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

                            pointA[0] = this->data.tetMesh.vertexDomainCoordinates[tet[i]][0];
                            pointA[1] = this->data.tetMesh.vertexDomainCoordinates[tet[i]][1];
                            pointA[2] = this->data.tetMesh.vertexDomainCoordinates[tet[i]][2];

                            pointB[0] = this->data.tetMesh.vertexDomainCoordinates[tet[j]][0];
                            pointB[1] = this->data.tetMesh.vertexDomainCoordinates[tet[j]][1];
                            pointB[2] = this->data.tetMesh.vertexDomainCoordinates[tet[j]][2];

                            pointC[0] = this->data.tetMesh.vertexDomainCoordinates[tet[k]][0];
                            pointC[1] = this->data.tetMesh.vertexDomainCoordinates[tet[k]][1];
                            pointC[2] = this->data.tetMesh.vertexDomainCoordinates[tet[k]][2];

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





    // Draw features

    //glBegin(GL_TRIANGLES);


    //for (const int triangleId : this->data.reebSpace2.trianglesPerSheet[4])
    //{
        //const std::set<int>& triangle = this->data.tetMesh.triangles[triangleId];

        //// Collect vertices into a GLfloat[3][3] array
        //GLfloat verts[3][3];
        //int idx = 0;
        //for (const int vertexId : triangle)
        //{
            //const std::array<GLfloat, 3>& vertexCoordinates = this->data.tetMesh.vertexDomainCoordinates[vertexId];

            //verts[idx][0] = vertexCoordinates[0];
            //verts[idx][1] = vertexCoordinates[1];
            //verts[idx][2] = vertexCoordinates[2];

            //idx++;
            //if (idx >= 3) break; // just in case
        //}

        //// Draw the triangle with proper normal
        //drawSolidTriangle(verts);
    //}


    //glEnd();


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
    float delta = -event->angleDelta().y();  // e.g. +120 or -120 per notch
    float zoomFactor = 1.0f + delta * 0.001f;  // tweak 0.001f to taste

    this->scale *= zoomFactor;

    if (this->scale < 1e-4f)  // Clamp to reasonable minimum
        this->scale = 1e-4f;

    update();
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
        //translateX += 0.1 * this->data.xdim;
        this->update();
    }
    if (event->key() == Qt::Key_A) {
        //translateY += 0.1 * this->data.xdim;
        this->update();
    }
    if (event->key() == Qt::Key_S) {
        //translateX -= 0.1 * this->data.xdim;
        this->update();
    }
    if (event->key() == Qt::Key_D) {
        //translateY -= 0.1 * this->data.xdim;
        this->update();
    }
}

void
TracerVisualiserWidget::mouseDoubleClickEvent(QMouseEvent* event)
{
    //this->data.faceFibers.clear();
    this->generateDisplayList();
    this->update();
}

void TracerVisualiserWidget::updateFiber(const std::vector<FiberPoint> &newFiberPoints)
{
    if (true == clearFibers)
    {
        this->faceFibers.clear();
    }

    this->faceFibers.insert(this->faceFibers.end(), newFiberPoints.begin(), newFiberPoints.end());
    this->generateDisplayList();
    this->update();
}

