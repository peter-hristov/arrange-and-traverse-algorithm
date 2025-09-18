#include "./CGALTypedefs.h"

#include "./TetMesh.h"
#include <random>
#include <queue>

void TetMesh::perturbRangeValues(const double &epsilon)
{
    static std::mt19937 gen(std::random_device{}());
    std::uniform_real_distribution<double> dist(-epsilon, epsilon);

    for (int i = 0; i < this->vertexCoordinatesF.size(); i++) 
    {
        this->vertexCoordinatesF[i] += dist(gen);
        this->vertexCoordinatesG[i] += dist(gen);
    }
}

void TetMesh::computeDomainBoundingBox()
{
    assert(!vertexDomainCoordinates.empty() && vertexDomainCoordinates[0].size() == 3);

    const auto &first = vertexDomainCoordinates[0];
    minX = maxX = first[0];
    minY = maxY = first[1];
    minZ = maxZ = first[2];

    for (const auto &v : vertexDomainCoordinates)
    {
        minX = std::min(minX, v[0]);
        maxX = std::max(maxX, v[0]);

        minY = std::min(minY, v[1]);
        maxY = std::max(maxY, v[1]);

        minZ = std::min(minZ, v[2]);
        maxZ = std::max(maxZ, v[2]);
    }
}

void TetMesh::computeRangeBoundingBox()
{
    assert(!vertexCoordinatesF.empty() && !vertexCoordinatesG.empty());

    const auto [minFIt, maxFIt] = std::minmax_element(vertexCoordinatesF.begin(), vertexCoordinatesF.end());
    const auto [minGIt, maxGIt] = std::minmax_element(vertexCoordinatesG.begin(), vertexCoordinatesG.end());

    minF = *minFIt;
    maxF = *maxFIt;
    minG = *minGIt;
    maxG = *maxGIt;
}

void TetMesh::computeBoundingBoxes()
{
    computeDomainBoundingBox();
    computeRangeBoundingBox();
}

void TetMesh::sortVertices()
{
    // Save the unsorted index of every points so that we can rename later
    std::vector<std::pair<CartesianPoint, int>> pointWithIndices(this->vertexCoordinatesG.size());

    for (int i = 0 ; i < pointWithIndices.size() ; i++)
    {
        pointWithIndices[i] = {CartesianPoint(this->vertexCoordinatesF[i], this->vertexCoordinatesG[i]), i};
    }

    // Sort all points
    std::sort(pointWithIndices.begin(), pointWithIndices.end(),
            [](const std::pair<CartesianPoint, int> &a, const std::pair<CartesianPoint, int>& b) {
            return CGAL::compare_xy(a.first, b.first) == CGAL::SMALLER;
            });


    // Set up the inverse index search
    std::vector<int> meshIDtoSortIndex(pointWithIndices.size());
    for (int i = 0 ; i < pointWithIndices.size() ; i++)
    {
        meshIDtoSortIndex[pointWithIndices[i].second] = i;
    }


    // Now we can swap things around

    // Set up copies of the originals for the swap, otherwise editin in place causes errors
    std::vector<std::array<float, 3>> vertexDomainCoordinatesOriginal = this->vertexDomainCoordinates;
    std::vector<double> vertexCoordinatesFOriginal = this->vertexCoordinatesF;
    std::vector<double> vertexCoordinatesGOriginal = this->vertexCoordinatesG;

    // Swap tet indices
    for (int i = 0 ; i < this->tetrahedra.size() ; i++)
    {
        for (int j = 0 ; j < this->tetrahedra[i].size() ; j++)
        {
            this->tetrahedra[i][j] = meshIDtoSortIndex[this->tetrahedra[i][j]];
        }
    }

    // Swap domain positions
    for (int i = 0 ; i < this->vertexDomainCoordinates.size() ; i++)
    {
        this->vertexDomainCoordinates[i] = vertexDomainCoordinatesOriginal[pointWithIndices[i].second];
    }
    
    // Swap range positions
    for (int i = 0 ; i < this->vertexCoordinatesF.size() ; i++)
    {
        this->vertexCoordinatesF[i] = vertexCoordinatesFOriginal[pointWithIndices[i].second];
        this->vertexCoordinatesG[i] = vertexCoordinatesGOriginal[pointWithIndices[i].second];
    }
}

void TetMesh::printMesh()
{
    // Print vertex domain coordinates
    std::cout << "Vertex domain coordinates: " << std::endl;
    for (int i = 0 ; i < this->vertexDomainCoordinates.size() ; i++)
    {
        printf("%d - (%f, %f, %f)\n", i, this->vertexDomainCoordinates[i][0], this->vertexDomainCoordinates[i][1], this->vertexDomainCoordinates[i][2]);
    }
    // Print vertex range coordinates
    std::cout << "Vertex range coordinates: " << std::endl;
    for (int i = 0 ; i < this->vertexDomainCoordinates.size() ; i++)
    {
        printf("%d - (%f, %f)\n", i, this->vertexCoordinatesF[i], this->vertexCoordinatesG[i]);
    }

    // Print tetrahedra
    std::cout << "Tetrahedra: " << std::endl;
    for (int i = 0 ; i < this->tetrahedra.size() ; i++)
    {
        printf("%d - (%d, %d, %d, %d)\n", i, this->tetrahedra[i][0], this->tetrahedra[i][1], this->tetrahedra[i][2], this->tetrahedra[i][3]);
    }
}


void TetMesh::computeSingularEdgeTypes()
{
    // Initialization
    for (const auto &edge : this->edges)
    {
        this->edgeSingularTypes[edge] = {};
    }

    std::map<std::array<int, 2>, DisjointSet<int>> upperLinkComponentsDS;
    for (const auto &[edge, vertices] : this->upperLink)
    {
        upperLinkComponentsDS[edge].initialize(vertices);
    }

    std::map<std::array<int, 2>, DisjointSet<int>> lowerLinkComponentsDS;
    for (const auto &[edge, vertices] : this->lowerLink)
    {
        lowerLinkComponentsDS[edge].initialize(vertices);
    }

    // otherwise the lower and upper link flip around
    for (int i = 0 ; i <  this->tetrahedra.size() ; i++)
    {
        const auto &tet = this->tetrahedra[i];

        // All pairs give you all six edges
        for (int a = 0 ; a < 4 ; a++)
        {
            for (int b = a + 1 ; b < 4 ; b++)
            {
                // Get the indices of the vertices for the edge
                int aIndex = tet[a];
                int bIndex = tet[b];

                // Make sure the vertices of the edge are in sorted order to have consistent orientation
                if (aIndex > bIndex)
                {
                    std::swap(aIndex, bIndex);
                }

                std::array<int, 2> edge = {aIndex, bIndex};

                std::array<int, 2> linkEdge;
                int k = 0;
                for (int v = 0; v < 4; ++v)
                {
                    if (v != a && v != b)
                    {
                        linkEdge[k++] = tet[v];
                    }
                }

                const auto &upperLink = this->upperLink.at(edge);
                const auto &lowerLink = this->lowerLink.at(edge);

                const bool isEdgeInUpperLink = upperLink.contains(linkEdge[0]) && upperLink.contains(linkEdge[1]);
                const bool isEdgeInLowerLink = lowerLink.contains(linkEdge[0]) && lowerLink.contains(linkEdge[1]);

                if (true == isEdgeInUpperLink)
                {
                    upperLinkComponentsDS[edge].unionElements(linkEdge[0], linkEdge[1]);
                }
                else if (true == isEdgeInLowerLink)
                {
                    lowerLinkComponentsDS[edge].unionElements(linkEdge[0], linkEdge[1]);
                }

            }
        }
    }

    this->regularEdgesNumber = 0;

    for (auto &[edge, type] : this->edgeSingularTypes)
    {
        const int upperLinkComponents = upperLinkComponentsDS[edge].getComponentRepresentatives().size();
        const int lowerLinkComponents = lowerLinkComponentsDS[edge].getComponentRepresentatives().size();

        // Definite edge
        if (0 == upperLinkComponents || 0 == lowerLinkComponents)
        {
            type = 0;
        }
        // Regular edge
        else if (1 == upperLinkComponents && 1 == lowerLinkComponents)
        {
            type = 1;
            this->regularEdgesNumber++;
        }
        // Indefinite edge
        else
        {
            type = 2;
        }

        //printf("Currently at the edge [%d, %d].\n", edge[0], edge[1]);

        //printf("The upper link has %d components and these vertices: ", upperLinkComponents);
        //for (const auto &v : this->upperLink[edge])
        //{
            //printf("%d ", v);
        //}
        //printf("\n");

        //printf("The lower link has %d components and these vertices: ", lowerLinkComponents);
        //for (const auto &v : this->lowerLink[edge])
        //{
            //printf("%d ", v);
        //}
        //printf("\n");
    }

    this->singularEdgesNumber = this->edges.size() - this->regularEdgesNumber;


    printf("There are %d singular edges out of %ld edges. That is %.2f%%.\n", this->singularEdgesNumber, edges.size(), 100.0 * (double)this->singularEdgesNumber / (double)edges.size());
}


void TetMesh::computeUpperLowerLinkAndStar()
{
    this->upperStarTrianglesNew.resize(this->edges.size());
    this->lowerStarTrianglesNew.resize(this->edges.size());

    // Initialize all maps with the empty set, so that all edges are covered
    for (auto &edge : this->edges)
    {
        this->upperLink[edge] = {};
        this->lowerLink[edge] = {};
        this->upperStarTriangles[edge] = {};
        this->lowerStarTriangles[edge] = {};
    }

    for (int i = 0 ; i <  this->tetrahedra.size() ; i++)
    {
        const auto &tet = this->tetrahedra[i];

        // All pairs give you all six edges
        for (int a = 0 ; a < 4 ; a++)
        {
            for (int b = a + 1 ; b < 4 ; b++)
            {
                // Get the indices of the vertices for the edge
                int aIndex = tet[a];
                int bIndex = tet[b];

                // Make sure the vertices of the edge are in sorted order to have consistent orientation
                if (aIndex > bIndex)
                {
                    std::swap(aIndex, bIndex);
                }

                // Search though the other two unused vertices
                for (int v = 0 ; v < 4 ; v++)
                {
                    // If the currect vertex is not one of two used to define the edge
                    if (v != a && v != b) 
                    {
                        const int vIndex = tet[v];
                        const std::array<int, 2> edge = {aIndex, bIndex};
                        const int edgeId = this->edgeIndices.at(edge);
                        const bool isUpperLink = this->isUpperLinkEdgeVertex(aIndex, bIndex, vIndex);

                        if (true == isUpperLink && false == this->upperLink[edge].contains(vIndex)) 
                        {
                            this->upperLink[edge].insert(vIndex);
                            this->upperStarTriangles[edge].push_back(triangleIndices.at({aIndex, bIndex, vIndex}));

                            this->upperStarTrianglesNew[edgeId].push_back(triangleIndices.at({aIndex, bIndex, vIndex}));
                        } 
                        else if (false == isUpperLink && false == this->lowerLink[edge].contains(vIndex)) 
                        {
                            this->lowerLink[edge].insert(vIndex);
                            this->lowerStarTriangles[edge].push_back(triangleIndices.at({aIndex, bIndex, vIndex}));

                            this->lowerStarTrianglesNew[edgeId].push_back(triangleIndices.at({aIndex, bIndex, vIndex}));
                        }
                    }
                }
            }
        }
    }
}


bool TetMesh::isUpperLinkEdgeVertex(int aIndex, int bIndex, int vIndex)
{
    // Make sure the vertices of the edge are in sorted order to have consistent orientation
    if (aIndex > bIndex)
    {
        std::swap(aIndex, bIndex);
    }

    // Define the two points that form the line
    const Point_2 a(this->vertexCoordinatesF[aIndex], this->vertexCoordinatesG[aIndex]);
    const Point_2 b(this->vertexCoordinatesF[bIndex], this->vertexCoordinatesG[bIndex]);

    // Define the test point
    const Point_2 v(this->vertexCoordinatesF[vIndex], this->vertexCoordinatesG[vIndex]);  // Change this to test different locations

    // Determine which half-plane r is in
    const CGAL::Orientation result = CGAL::orientation(a, b, v);

    //printf("Checking line (%d, %d) against vertex %d\n", aIndex, bIndex, vIndex);

    //std::cout << "a coords = " << a << std::endl;
    //std::cout << "b coords = " << b << std::endl;
    //std::cout << "v coords = " << v << std::endl;


    // Upper link = left
    if (result == CGAL::LEFT_TURN) {
        return true;
        //std::cout << "Point r is in the LEFT half-plane.\n";
    // Lower link = right
    } else if (result == CGAL::RIGHT_TURN) {
        return false;
        //std::cout << "Point r is in the RIGHT half-plane.\n";
    // This should not happen for generic maps
    } else {
        //std::cout << "Point r is on the line.\n";
        throw std::runtime_error("Input data is degeneate, a triangle is mapped to a line.");       
        //assert(false);
    }

    // Paranoid assert
    assert(false);
}



void TetMesh::computeCombinatorialStructure()
{
    // Compute all unique edges
    std::set<std::array<int, 2>> allEdges;
    this->skeletonGraph.resize(this->vertexCoordinatesF.size());

    for (int i = 0 ; i <  this->tetrahedra.size() ; i++)
    {
        const auto &tet = this->tetrahedra[i];

        for (int a = 0 ; a < 4 ; a++)
        {
            for (int b = a + 1 ; b < 4 ; b++)
            {
                int aIndex = tet[a];
                int bIndex = tet[b];

                // Make sure the vertices of the edge are in sorted order to have consistent orientation
                if (aIndex > bIndex)
                {
                    std::swap(aIndex, bIndex);
                }

                allEdges.insert({aIndex, bIndex});
                skeletonGraph[aIndex].insert(bIndex);
                skeletonGraph[bIndex].insert(aIndex);
            }
        }
    }

    for (const auto &edge: allEdges)
    {
        this->edges.push_back(edge);
        this->edgeIndices[edge] = this->edges.size() - 1;
    }


    // Compute all unique triangles
    std::set<std::set<int>> allTriangles;

    for (const std::array<int, 4> tet : this->tetrahedra)
    {
        std::set<std::set<int>> triangles;

        for (int a = 0 ; a < 4 ; a++)
        {
            for (int b = a + 1 ; b < 4 ; b++)
            {
                for (int c = b + 1 ; c < 4 ; c++)
                {
                    int aIndex = tet[a];
                    int bIndex = tet[b];
                    int cIndex = tet[c];
                    allTriangles.insert({aIndex, bIndex, cIndex});
                }
            }
        }
    }

    for (const std::set<int> triangle : allTriangles)
    {
        this->triangles.push_back(triangle);
        this->triangleIndices[triangle] = this->triangles.size() - 1;
    }


    // Compute which triangles are tet-adjacent (are both the faces of the same tet)
    this->tetIncidentTriangles.resize(allTriangles.size());

    for (const std::array<int, 4> tet : this->tetrahedra)
    {
        std::set<std::set<int>> triangles;

        for (int a = 0 ; a < 4 ; a++)
        {
            for (int b = a + 1 ; b < 4 ; b++)
            {
                for (int c = b + 1 ; c < 4 ; c++)
                {
                    const int aIndex = tet[a];
                    const int bIndex = tet[b];
                    const int cIndex = tet[c];
                    triangles.insert({aIndex, bIndex, cIndex});
                }
            }
        }

        // Connect all the triangles together
        for(const std::set<int> t1 : triangles)
        {
            for(const std::set<int> t2 : triangles)
            {
                int t1Index = this->triangleIndices[t1];
                int t2Index = this->triangleIndices[t2];

                if (t1Index == t2Index) { continue; }

                //this->tetIncidentTriangles[t1Index].push_back(t2Index);
                //this->tetIncidentTriangles[t2Index].push_back(t1Index);

                this->tetIncidentTriangles[t1Index].insert(t2Index);
                this->tetIncidentTriangles[t2Index].insert(t1Index);
            }
        }
    }
}

int TetMesh::computeSingularSetConnectivity()
{
    DisjointSet<int> singularComponents;

    for (const auto &[edge, type] : this->edgeSingularTypes)
    {
        if (type != 1)
        {
            if (false == singularComponents.data.contains(edge[0]))
            {
                singularComponents.addElement(edge[0]);
            }

            if (false == singularComponents.data.contains(edge[1]))
            {
                singularComponents.addElement(edge[1]);
            }
        }
    }

    for (const auto &[edge, type] : this->edgeSingularTypes)
    {
        if (type != 1)
        {
            singularComponents.unionElements(edge[0], edge[1]);
        }
    }

    std::map<int, int> componentSize;

    for (const auto &[vertexId, internalId] : singularComponents.data)
    {
        componentSize[singularComponents.findElement(vertexId)]++;
    }

    //printf("Ther are %d singular vertices.\n", skeletonGraph.data.size());
    printf("The singular set has %d connected components.\n", singularComponents.countConnectedComponents());

    return singularComponents.countConnectedComponents();

    //int totalSize = 0;
    //for (const auto &[rootId, size] : componentSize)
    //{
        //totalSize += size;
        //printf("Component with ID %d has %d components.\n", rootId, size);
    //}

    //if (totalSize != skeletonGraph.data.size())
    //{
        //std::cerr << "The sizes do not match up!";
    //}
}

void TetMesh::singularTraversalBFS(const std::vector<int> &roots, std::vector<bool> &visited)
{
    std::queue<int> q;

    for (int r : roots)
    {
        q.push(r);
        visited[r] = true;
    }

    //printf("--------------------------------------\n");
    //printf("New component\n");

    while (false == q.empty())
    {
        int currentVertexId = q.front();
        q.pop();

        //printf("At vertex (%f, %f)\n", this->vertexCoordinatesF[currentVertexId], this->vertexCoordinatesG[currentVertexId]);

        for (int neighbourVertexId : this->skeletonGraph[currentVertexId])
        {
            if (false == visited[neighbourVertexId] && this->getEdgeType({currentVertexId, neighbourVertexId}) != 1)
            {
                visited[neighbourVertexId] = true;
                q.push(neighbourVertexId);
            }
        }
    }
}

int TetMesh::computeSingularSetConnectivity2()
{
    std::vector<bool> visited(this->vertexCoordinatesF.size(), true);

    for (const auto &[edge, type] : this->edgeSingularTypes)
    {
        if (type != 1)
        {
            visited[edge[0]] = false;
            visited[edge[1]] = false;
        }
    }

    int componentCount = 0;
    for (int i = 0 ; i < visited.size() ; i++)
    {
        if (visited[i] == false)
        {
            componentCount++;
            singularTraversalBFS({i}, visited);
        }
    }

    printf("The singular set has %d connected components.\n", componentCount);

    return componentCount;
}




const std::vector<int> TetMesh::findShortestPath(const std::vector<int> &source, const std::set<int> &sink) const
{
    std::vector<int> parent(this->vertexCoordinatesF.size(), -1);

    printf("\nThe source is              : ");
    std::queue<int> q;
    for (int s : source)
    {
        printf("%d ", s);
        q.push(s);
        parent[s] = s;
    }

    printf("\nThe sink is                : ");
    for (int s : sink)
    {
        printf("%d ", s);
    }

    int firstFound = -1;

    while (false == q.empty())
    {
        int currentVertexId = q.front();
        q.pop();
        
        for (int neighbourVertexId : this->skeletonGraph[currentVertexId])
        {
            if (sink.contains(neighbourVertexId))
            {
                firstFound = neighbourVertexId;
                parent[firstFound] = currentVertexId;
                break;
            }

            //if (-1 == parent[neighbourVertexId] && this->getEdgeType({currentVertexId, neighbourVertexId}) == 1)
            if (-1 == parent[neighbourVertexId])
            {
                parent[neighbourVertexId] = currentVertexId;
                q.push(neighbourVertexId);
            }

        }

        if (firstFound != -1)
        {
            break;
        }
    }

    if (firstFound == -1)
    {
        throw std::runtime_error("Could not find a path between the source and the sink.");       
    }

    // Assemble path
    std::vector<int> shortestPath;

    int currentVertexId = firstFound;

    printf("\nThe path is                : ");
    while (parent[currentVertexId] != currentVertexId)
    {
        printf("%d ", currentVertexId);
        shortestPath.push_back(currentVertexId);
        currentVertexId = parent[currentVertexId];
    }
    printf("%d", currentVertexId);

    shortestPath.push_back(currentVertexId); // finally push the root
    std::reverse(shortestPath.begin(), shortestPath.end());

    return shortestPath;
}

void TetMesh::markPseudoSingularEdges(const std::vector<int> &shortestPath)
{
    printf("\nMarking the following edges as pseudo-singular : ");
    for (int i = 0 ; i < shortestPath.size() - 1 ; i++)
    {
        int vertexA = shortestPath[i];
        int vertexB = shortestPath[i+1];

        if (vertexA > vertexB)
        {
            std::swap(vertexA, vertexB);
        }

        if (this->edgeSingularTypes[{vertexA, vertexB}] == 1)
        {
            printf("[%d, %d] ", vertexA, vertexB);
            this->edgeSingularTypes[{vertexA, vertexB}] = -1;

            // A regular edge becomes pseudo-singular
            this->pseudoSingularEdgesNumber++;
            this->regularEdgesNumber--;
        }

    }
    printf("\n");
    printf("\n");
    printf("\n");
}



const std::vector<int>& TetMesh::getMinusTriangles(const int &edgeId, const bool &isDirectionLowerToUpper) const
{
    if (isDirectionLowerToUpper)
    {
        return this->lowerStarTrianglesNew[edgeId];
    }
    else
    {
        return this->upperStarTrianglesNew[edgeId];
    }
}

const std::vector<int>& TetMesh::getPlusTriangles(const int &edgeId, const bool &isDirectionLowerToUpper) const
{
    if (isDirectionLowerToUpper)
    {
        return this->upperStarTrianglesNew[edgeId];
    }
    else
    {
        return this->lowerStarTrianglesNew[edgeId];
    }
}
