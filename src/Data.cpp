#include "Data.h"

#include <fstream>
#include <sstream>
#include <utility>

using namespace std;

void Data::computeTetExitPoints(const GLfloat u, const GLfloat v, const std::vector<float> color)
{
    this->tetsWithFibers = vector<bool>(this->tetrahedra.size(), false);

    // For every tet, compute the two exit points
    for(size_t tetId = 0 ; tetId < this->tetrahedra.size(); tetId++)
    {
        const auto tet = this->tetrahedra[tetId];

        // For every triangle in every tet, get the fiber in it
        for(int i = 0 ; i < 4 ; i++)
        {
            for(int j = i + 1 ; j < 4 ; j++)
            {
                for(int k = j + 1 ; k < 4 ; k++)
                {
                    float x1 = this->vertexCoordinatesF[tet[i]];
                    float y1 = this->vertexCoordinatesG[tet[i]];

                    float x2 = this->vertexCoordinatesF[tet[j]];
                    float y2 = this->vertexCoordinatesG[tet[j]];

                    float x3 = this->vertexCoordinatesF[tet[k]];
                    float y3 = this->vertexCoordinatesG[tet[k]];

                    float det = (x1 - x3) * (y2 - y3) - (x2 - x3) * (y1 - y3);

                    float alpha = ((y2 - y3) * (u - x3) + (x3 - x2) * (v - y3)) / det;
                    float beta = ((y3 - y1) * (u - x3) + (x1 - x3) * (v - y3)) / det;
                    float gamma = 1 - alpha - beta;

                    // Are we inside the triangle. We exclude the 0 and 1 because weird things happen there
                    if (alpha > 0 && beta > 0 && gamma > 0 && alpha < 1 && beta < 1 && gamma < 1)
                    {
                        //printf("In triangle %ld, %ld, %ld in tet %ld.\n", tet[i], tet[j], tet[k], tetId);
                        //printf("In triangle (%f, %f) | (%f, %f) | (%f, %f) comparing with point (%f, %f) and alpha = %f, betta = %f, gamma = %f.\n", x1, y1, x2, y2, x3, y3, x, y, alpha, betta, gamma);

                        FaceFiberPoint fb(alpha, beta, {
                                this->vertexDomainCoordinates[tet[i]],
                                this->vertexDomainCoordinates[tet[j]],
                                this->vertexDomainCoordinates[tet[k]],
                                },
                                color);

                        this->faceFibers.push_back(fb);
                        this->tetsWithFibers[tetId] = true;
                    }
                }
            }
        }
    }
}

void
Data::computeMinMaxRangeDomainCoordinates()
{

    // Compute the min/max domain coordinates
    this->minX = this->vertexDomainCoordinates[0][0];
    this->maxX = this->vertexDomainCoordinates[0][0];

    this->minY = this->vertexDomainCoordinates[0][1];
    this->maxY = this->vertexDomainCoordinates[0][1];

    this->minZ = this->vertexDomainCoordinates[0][2];
    this->maxZ = this->vertexDomainCoordinates[0][2];

    for (int i = 0 ; i < this->vertexDomainCoordinates.size() ; i++)
    {
        this->minX = std::min(this->minX, this->vertexDomainCoordinates[i][0]);
        this->maxX = std::max(this->maxX, this->vertexDomainCoordinates[i][0]);

        this->minY = std::min(this->minY, this->vertexDomainCoordinates[i][1]);
        this->maxY = std::max(this->maxY, this->vertexDomainCoordinates[i][1]);

        this->minZ = std::min(this->minZ, this->vertexDomainCoordinates[i][2]);
        this->maxZ = std::max(this->maxZ, this->vertexDomainCoordinates[i][2]);
    }


    // Compute the min/max range coordinates
    this->minF = this->vertexCoordinatesF[0];
    this->maxF = this->vertexCoordinatesF[0];

    this->minG = this->vertexCoordinatesG[0];
    this->maxG = this->vertexCoordinatesG[0];

    for (int i = 0 ; i < this->vertexCoordinatesF.size() ; i++)
    {
        this->minF = std::min(this->minF, this->vertexCoordinatesF[i]);
        this->maxF = std::max(this->maxF, this->vertexCoordinatesF[i]);

        this->minG = std::min(this->minG, this->vertexCoordinatesG[i]);
        this->maxG = std::max(this->maxG, this->vertexCoordinatesG[i]);
    }

    // Add some padding to the range coordinates for better visibility
    //this->minF -= .2;
    //this->maxF += .2;
    //this->minG -= .2;
    //this->maxG += .2;
}


void
Data::readData(string filename)
{
    // Set deault names for the range axis
    this->longnameF = "f";
    this->longnameG = "g";

    // Open data file
    std::ifstream dataFile (filename);
    if (false == dataFile.is_open()) { throw "Could not open data file."; }

    // Read in data in a string and skip the comments
    string rawStringData;
    string myline;
    while (dataFile) {
        std::getline (dataFile, myline);
        if (myline[0] == '#')
        {
            //std::cout << myline << '\n';
        }
        else
        {
            rawStringData += " " + myline;
        }
    }

    // Set up the inputstream from the string
    std::istringstream dataStream(rawStringData);

    // Read in the number of vertices and tets
    int numVertices, numTets;
    dataStream >> numVertices >> numTets;

    // Initialize all the data arrays
    this->vertexCoordinatesF = std::vector<GLfloat>(numVertices, 0);
    this->vertexCoordinatesG = std::vector<GLfloat>(numVertices, 0);
    this->tetrahedra = std::vector<std::vector<size_t>>(numTets, {0, 0, 0, 0});
    this->vertexDomainCoordinates = std::vector<std::vector<GLfloat>>(numVertices, {0, 0, 0});

    // Read in the domain coordinates
    for  (int i = 0 ; i < numVertices ; i++)
    {
        dataStream >> this->vertexDomainCoordinates[i][0];
        dataStream >> this->vertexDomainCoordinates[i][1];
        dataStream >> this->vertexDomainCoordinates[i][2];
    }

    // Read in the range coordinates
    for  (int i = 0 ; i < numVertices ; i++)
    {
        dataStream >> this->vertexCoordinatesF[i];
        dataStream >> this->vertexCoordinatesG[i];
    }
    
    // Read in the tetrahedron configuration
    for  (int i = 0 ; i < numTets ; i++)
    {
        dataStream >> this->tetrahedra[i][0];
        dataStream >> this->tetrahedra[i][1];
        dataStream >> this->tetrahedra[i][2];
        dataStream >> this->tetrahedra[i][3];
    }
}

void
Data::readDataGrid(const string filename)
{
    // Set deault names for the range axis
    this->longnameF = "f";
    this->longnameG = "g";

    // Open data file
    std::ifstream dataFile (filename);
    if (false == dataFile.is_open()) { throw "Could not open data file."; }

    dataFile >> xDim >> yDim >> zDim;

    // Add vertex domain coordinates
    //for (int k = 0 ; k < this->zdim ; k++)
    //{
        //for (int j = 0 ; j < this->ydim ; j++)
        //{
            //for (int i = 0 ; i < this->xdim ; i++)
            //{
                //this->vertexDomainCoordinates.push_back({static_cast<float>(i), static_cast<float>(j), static_cast<float>(k)});
            //}
        //}
    //}

    // Add vertex domain coordinates
    for (int k = 0 ; k < this->zDim ; k++)
    {
        for (int j = 0 ; j < this->yDim ; j++)
        {
            for (int i = 0 ; i < this->xDim ; i++)
            {
                // Interpolated the coordinates
                this->vertexDomainCoordinates.push_back({static_cast<float>(i), static_cast<float>(j), static_cast<float>(k)});
            }
        }
    }

    // Add vertex range coordinates
    for (int k = 0 ; k < this->zDim ; k++)
    {
        for (int j = 0 ; j < this->yDim ; j++)
        {
            for (int i = 0 ; i < this->xDim ; i++)
            {
                float value;
                dataFile >> value;

                this->vertexCoordinatesF.push_back(value);
                this->vertexCoordinatesG.push_back(static_cast<float>(j));
            }
        }
    }



    // Add tets
    for (int i = 0 ; i < this->xDim - 1 ; i++)
    {
        for (int j = 0 ; j < this->yDim - 1 ; j++)
        {
            for (int k = 0 ; k < this->zDim - 1 ; k++)
            {
                this->addTetsForCube(i, j, k);
            }
        }
    }
}

size_t Data::trippleToIndex(const size_t i, const size_t j, const size_t k)
{
    return this->xDim * this->yDim * k + this->xDim * j + i;
}

void Data::addTetsForCube(const size_t i, const size_t j, const size_t k)
{
    this->tetrahedra.push_back({this->trippleToIndex(i, j, k), this->trippleToIndex(i+1, j+1, k+1), this->trippleToIndex(i, j+1, k), this->trippleToIndex(i+1, j+1, k)});
    this->tetrahedra.push_back({this->trippleToIndex(i, j, k), this->trippleToIndex(i+1, j+1, k+1), this->trippleToIndex(i, j+1, k), this->trippleToIndex(i, j+1, k+1)});
    this->tetrahedra.push_back({this->trippleToIndex(i, j, k), this->trippleToIndex(i+1, j+1, k+1), this->trippleToIndex(i, j, k+1), this->trippleToIndex(i, j+1, k+1)});

    this->tetrahedra.push_back({this->trippleToIndex(i, j, k), this->trippleToIndex(i+1, j+1, k+1), this->trippleToIndex(i+1, j, k), this->trippleToIndex(i+1, j+1, k)});
    this->tetrahedra.push_back({this->trippleToIndex(i, j, k), this->trippleToIndex(i+1, j+1, k+1), this->trippleToIndex(i+1, j, k), this->trippleToIndex(i+1, j, k+1)});
    this->tetrahedra.push_back({this->trippleToIndex(i, j, k), this->trippleToIndex(i+1, j+1, k+1), this->trippleToIndex(i, j, k+1), this->trippleToIndex(i+1, j, k+1)});
}
