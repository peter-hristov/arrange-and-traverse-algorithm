//#include <qpoint.h>
#include "src/CGALTypedefs.h"
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
#include "./Fiber.h"

#include <vtkPointData.h>
#include <vtkFloatArray.h>

using namespace std;


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

    if (this->enableLighting)
    {
        glEnable(GL_LIGHTING);
        glEnable(GL_LIGHT0);
        glDisable(GL_CULL_FACE);
        glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

        //glEnable(GL_COLOR_MATERIAL);
        //glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
    }

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

void
TracerVisualiserWidget::setMaterial(GLfloat red, GLfloat green, GLfloat blue, GLfloat alpha, GLfloat shininess)
{
    //glMaterialfv(GL_FRONT, GL_AMBIENT, &(vector<GLfloat>({ red, green, blue, alpha })[0]));
    //glMaterialfv(GL_FRONT, GL_DIFFUSE, &(vector<GLfloat>({ red, green, blue, alpha })[0]));
    //glMaterialfv(GL_FRONT, GL_SPECULAR, &(vector<GLfloat>({ red, green, blue, alpha })[0]));
    //glMaterialf(GL_FRONT, GL_SHININESS, shininess);

    GLfloat mat[4] = { red, green, blue, alpha };

    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat);

    // Use a low white specular instead of the full material colour
    GLfloat specular[4] = { 0.1f, 0.1f, 0.1f, 1.0f };
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular);
    //glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat); // optional

    glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, shininess);
}

// Draw fibers
void
TracerVisualiserWidget::generateDisplayList()
{
    glDeleteLists(displayListIndex, 1);
    displayListIndex = glGenLists(1);
    glNewList(displayListIndex, GL_COMPILE);

    //setMaterial(1, 0, 0, 1.0, 0.0);

    if (this->drawFiberSurfaceFeatures)
    {
        glBegin(GL_TRIANGLES);
        {
            //for(const auto &faceFiber : this->faceFibers)
            for(int i = 0 ; i < this->faceFiberSurfaceFeatures.size() ; i+=3)
            {
                const auto &faceFiber = this->faceFiberSurfaceFeatures[i];
                const auto &faceFiber2 = this->faceFiberSurfaceFeatures[i+1];
                const auto &faceFiber3 = this->faceFiberSurfaceFeatures[i+2];

                if (this->enableLighting)
                {
                    setMaterial(faceFiber.colour[0], faceFiber.colour[1], faceFiber.colour[2], this->featureOpacity, 1.0);

                    GLfloat vertices[3][3] = {
                        {faceFiber.point[0], faceFiber.point[1], faceFiber.point[2]}, 
                        {faceFiber2.point[0], faceFiber2.point[1], faceFiber2.point[2]}, 
                        {faceFiber3.point[0], faceFiber3.point[1], faceFiber3.point[2]}, 

                    };

                    std::array<GLfloat, 3> normal = this->computeTriangleNormal(faceFiber.point.data(), faceFiber2.point.data(), faceFiber3.point.data());


                    // Set normal for OpenGL
                    glNormal3fv(normal.data());
                }
                else
                {
                    glColor3fv(faceFiber.colour.data());

                }

                glVertex3fv(faceFiber.point.data());
                glVertex3fv(faceFiber2.point.data());
                glVertex3fv(faceFiber3.point.data());

            }
        }
        glEnd();

    }


    if (this->drawFiberSurfaces)
    {
        glBegin(GL_TRIANGLES);
        {
            //for(const auto &faceFiber : this->faceFibers)
            for (const auto fiberSurface : this->data.surfaceMeshes)
            {
                const auto sheetIdMap = fiberSurface.sheetId();
                for (const auto triangle : fiberSurface.mesh.faces())
                {

                    // Get the vertex coordinates
                    std::vector<std::array<GLfloat, 3>> vertices;
                    for (auto vertex : fiberSurface.mesh.vertices_around_face(fiberSurface.mesh.halfedge(triangle))) 
                    {
                        auto& p = fiberSurface.mesh.point(vertex);
                        vertices.push_back({static_cast<GLfloat>(p.x()), static_cast<GLfloat>(p.y()), static_cast<GLfloat>(p.z())});
                    }

                    // Get the colour
                    const int sheetId = sheetIdMap[triangle];
                    std::array<float, 3> triangleColour;
                    if (sheetId == -1)
                    {
                        triangleColour = {1.0, 1.0, 0.0};
                    }
                    else
                    {
                        const int sheetSortId = data.reebSpace2.sheetOrder.at(sheetId);
                        triangleColour = fiber::fiberColours[sheetSortId % fiber::fiberColours.size()];
                    }


                    // Set colour and compute normal
                    if (this->enableLighting)
                    {
                        setMaterial(triangleColour[0], triangleColour[1], triangleColour[2], fsOpacity, 1.0);
                        std::array<GLfloat, 3> normal = this->computeTriangleNormal(vertices[0].data(), vertices[1].data(), vertices[2].data());
                        glNormal3fv(normal.data());
                    }
                    else
                    {
                        glColor3f(triangleColour[0], triangleColour[1], triangleColour[2]);
                    }

                    // Push out the vertices
                    glVertex3fv(vertices[0].data());
                    glVertex3fv(vertices[1].data());
                    glVertex3fv(vertices[2].data());
                }

            }
        }
        glEnd();

    }




    if (this->drawFibers)
    {
        glDisable(GL_LIGHTING);
        // Draw Fiber
        glBegin(GL_LINES);
        {
            for(const auto &faceFiber : this->faceFibers)
            {
                if (this->enableLighting)
                {
                    //glColor3fv(faceFiber.colour.data());
                    glColor4f(faceFiber.colour[0], faceFiber.colour[1], faceFiber.colour[2], this->fiberOpcity);

                }
                else
                {
                    setMaterial(faceFiber.colour[0], faceFiber.colour[1], faceFiber.colour[2], 1.0, 1.0);
                }

                glVertex3fv(faceFiber.point.data());
            }
        }
        glEnd();
        glEnable(GL_LIGHTING);

    }




    // Draw fiber endpoints (in every tet)
    //for(const auto &faceFiber : this->faceFibers)
    //{
        //glColor3f(1.0f, 1.0f, 1.0f);

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

std::array<GLfloat, 3> TracerVisualiserWidget::computeTriangleNormal(const GLfloat* v0, const GLfloat* v1, const GLfloat* v2)
{
    // Edge vectors
    GLfloat ax = v1[0] - v0[0];
    GLfloat ay = v1[1] - v0[1];
    GLfloat az = v1[2] - v0[2];

    GLfloat bx = v2[0] - v0[0];
    GLfloat by = v2[1] - v0[1];
    GLfloat bz = v2[2] - v0[2];

    // Cross product a × b
    GLfloat nx = ay * bz - az * by;
    GLfloat ny = az * bx - ax * bz;
    GLfloat nz = ax * by - ay * bx;

    // Normalize
    GLfloat len = std::sqrt(nx*nx + ny*ny + nz*nz);
    if (len > 0.0f)
    {
        nx /= len;
        ny /= len;
        nz /= len;
    }

    return { nx, ny, nz };
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
    glTranslatef(0.0, 0.0, -1 * scale / 1.2);

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

    glDisable(GL_LIGHTING);
    glColor3f(0.5, 0.5, 0.5);
    drawWiredCube(vertices);
    glEnable(GL_LIGHTING);


    glColor3f(1, 1, 1);

        glDisable(GL_LIGHTING);
    if (true == this->drawEdges)
    {
        //glColor4f(1, 1, 1, this->edgeOpacity);
        glColor3f(0.5, 0.5, 0.5);

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

                        //const int edgeType = data.tetMesh.edgeSingularTypes.at({aIndex, bIndex});

                        //if (edgeType == 0)
                        //{
                            //glColor3f(0.6f, 0.85f, 0.85f);
                        //}
                        //else if(edgeType == 2)
                        //{
                            //glColor3f(1.0f, 0.647f, 0.0f);
                        //}

                        //if (edgeType !=1)
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
        glEnable(GL_LIGHTING);

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

                //glColor4f(1, 1, 1, this->vertexOpacity);
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
        //glColor4f(1, 1, 1, this->faceOpacity);

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
        glEnd();
    }


    this->renderMolecule();



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
        //initialX = event->localPos().x();
        //initialY = event->localPos().y();

        const int sheetId = pickSegment(event->x(), event->y());

        if (sheetId >= 0)
        {
            this->selectedSheetIds.insert(sheetId);

            for (const int id : this->selectedSheetIds)
            {
                //std::cout << "Sheet " << id << " has area " << data.reebSpace2.sheetArea[id] << " which is a ratio of : " << 100.0 * data.reebSpace2.sheetArea[id] <<  std::endl;
                printf("Sheet %d had area %.2f (which is %.2f%%).\n", id, data.reebSpace2.sheetArea[id], 100.0 * data.reebSpace2.sheetAreaProportion[id]);
            }

            static_cast<PlotWidget*>(this->sibling)->staticReebSpaceCache = nullptr;
            this->sibling->update();

            if (auto *window = qobject_cast<TracerVisualiserWindow*>(this->parent()->parent())) {
                window->updateSelectedSheets(this->selectedSheetIds);
            }

            update();
        }
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
        //float x = event->localPos().x();
        //float y = event->localPos().y();

        //translateX -= (initialX - x) / 10;
        //translateY += (initialY - y) / 10;

        //initialX = event->localPos().x();
        //initialY = event->localPos().y();

        //this->update();
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

void TracerVisualiserWidget::updateFiberSurface()
{
    this->buildAABBTree();
    this->generateDisplayList();
    this->update();
}

void TracerVisualiserWidget::updateFiberSurfaceFeatures(const std::vector<FiberPoint> &newFiberPoints)
{
    this->faceFiberSurfaceFeatures = newFiberPoints;
    this->buildAABBTree();

    this->generateDisplayList();
    this->update();
}

void TracerVisualiserWidget::clearFiber()
{
    this->faceFibers = {};
    this->generateDisplayList();
    this->update();
}

void TracerVisualiserWidget::clearFiberSurface()
{
    this->data.surfaceMeshes.clear();
    this->data.surfaceMeshes.shrink_to_fit();
    this->generateDisplayList();
    this->update();
}

void TracerVisualiserWidget::clearFiberSurfaceSheets()
{
    this->faceFiberSurfaceFeatures = {};
    this->generateDisplayList();
    this->update();
}

void TracerVisualiserWidget::renderMolecule()
{
    if (!data.molecule)
    {
        return;
    }

    vtkPoints*    points = data.molecule->GetPoints();
    vtkCellArray* lines  = data.molecule->GetLines();

    glDisable(GL_LIGHTING);
    glColor3f(1.0f, 1.0f, 1.0f);



    glBegin(GL_LINES);

    lines->InitTraversal();
    vtkNew<vtkIdList> idList;

    while (lines->GetNextCell(idList))
    {
        for (vtkIdType i = 0; i < idList->GetNumberOfIds() - 1; i++)
        {
            double p0[3], p1[3];
            points->GetPoint(idList->GetId(i),     p0);
            points->GetPoint(idList->GetId(i + 1), p1);

            glVertex3d(p0[0], p0[1], p0[2]);
            glVertex3d(p1[0], p1[1], p1[2]);
        }
    }

    glEnd();

    // Get the point data arrays
    vtkFloatArray* colourArray = vtkFloatArray::SafeDownCast(data.molecule->GetPointData()->GetArray("atom_color"));
    vtkFloatArray* radiusArray = vtkFloatArray::SafeDownCast(data.molecule->GetPointData()->GetArray("atom_radius")); // keeping your spelling
    if (!colourArray || !radiusArray) { return; }

    for (vtkIdType i = 0; i < points->GetNumberOfPoints(); ++i)
    {
        // Read radius; skip if zero
        const float radius = static_cast<float>(radiusArray->GetValue(i)) * 0.2;
        if (radius < 0.000001f) { continue; }

        // Read position
        double pos[3];
        points->GetPoint(i, pos);

        // Read colour (3 floats)
        const float r = static_cast<float>(colourArray->GetComponent(i, 0));
        const float g = static_cast<float>(colourArray->GetComponent(i, 1));
        const float b = static_cast<float>(colourArray->GetComponent(i, 2));

        // Render sphere
        glColor3f(r, g, b);
        //setMaterial(r, g, b, 1.0f, 20.0f);

        glPushMatrix();
        {
            glTranslatef(
                static_cast<float>(pos[0]),
                static_cast<float>(pos[1]),
                static_cast<float>(pos[2]));

            GLUquadric* sphere = gluNewQuadric();
            gluSphere(sphere, radius, 16, 16);
            gluDeleteQuadric(sphere);   // use gluDeleteQuadric, not delete
        }
        glPopMatrix();
    }
    glEnable(GL_LIGHTING);

}

void TracerVisualiserWidget::buildAABBTree()
{
    aabbTriangleTrees.clear();
    aabbTriangleTrees.reserve(this->data.surfaceMeshes.size());
    for (const auto& fiberSurface : this->data.surfaceMeshes)
    {
        aabbTriangleTrees.emplace_back(
                faces(fiberSurface.mesh).first,
                faces(fiberSurface.mesh).second,
                fiberSurface.mesh
                );
        aabbTriangleTrees.back().accelerate_distance_queries();
    }
    //qDebug() << "AABB tree built with" << pickingTriangles.size() << "triangles";
}

// Claud generated code
int TracerVisualiserWidget::pickSegment(int mouseX, int mouseY)
{
    if (aabbTriangleTrees.empty()) return -1;

    makeCurrent();

    GLint    viewport[4];
    GLdouble mv[16], proj[16];
    glGetIntegerv(GL_VIEWPORT,        viewport);
    glGetDoublev(GL_MODELVIEW_MATRIX,  mv);
    glGetDoublev(GL_PROJECTION_MATRIX, proj);

    double winY = viewport[3] - mouseY;

    GLdouble nx, ny, nz, fx, fy, fz;
    gluUnProject(mouseX, winY, 0.0, mv, proj, viewport, &nx, &ny, &nz);
    gluUnProject(mouseX, winY, 1.0, mv, proj, viewport, &fx, &fy, &fz);

    CartesianPoint_3 origin(nx, ny, nz);
    CartesianPoint_3 target(fx, fy, fz);
    CartesianKernel::Ray_3 ray(origin, target);

    double bestDist = std::numeric_limits<double>::max();
    int bestSheetId = -1;

    // Find the closest one from all the meshes
    for (int m = 0; m < aabbTriangleTrees.size(); ++m)
    {
        auto hit = aabbTriangleTrees[m].first_intersection(ray);
        if (hit)
        {
            // intersection point is in hit->first
            const CartesianPoint_3* p = std::get_if<CartesianPoint_3>(&(hit->first));
            if (p)
            {
                double dist = CGAL::squared_distance(origin, *p);
                if (dist < bestDist)
                {
                    std::cerr << "Chaing closest best\n";
                    bestDist = dist;
                    bestSheetId = this->data.surfaceMeshes[m].sheetId()[hit->second];
                }
            }
        }
    }
    return bestSheetId;
}

