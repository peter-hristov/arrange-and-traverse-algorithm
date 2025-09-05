#pragma once

#include <cassert>
#include <vector>
//#include <unordered_set>
#include <set>
#include <queue>

#include "./TetMesh.h"
#include "./DisjointSet.h"

class PreimageGraph
{
    public:
        std::map<int, int> componentRoot;

        PreimageGraph() { }

        // Initialize vertices (triangles), the edges are implicit
        PreimageGraph(const std::map<int, int> &_componentRoot)  : componentRoot(_componentRoot) { }

        // @TODO Refactor with a better clear
        void clear()
        {
            this->componentRoot  = std::map<int, int>();
        }


        bool areEqual(PreimageGraph &p)
        {
            // Group elements with the same root
            std::map<int, std::set<int>> groupedOurs = this->groupComponents();
            std::map<int, std::set<int>> groupedTheirs = p.groupComponents();

            // Make sure there is an equal number of roots
            if (groupedOurs.size() != groupedTheirs.size())
            {
                std::cerr << "Size of preimage graphs is not equal!.";
                return false;
            }

            for (const auto &[root, elements] : groupedOurs)
            {
                const int firstElement = *elements.begin();

                // Find which component it's in
                const int componentId = p.componentRoot.at(firstElement);

                if (elements != groupedTheirs[componentId])
                {
                    std::cerr << "Corresponding components are not equal!\n\n";
                    std::cerr << "Our Component root             " << root << ": ";
                    for (int k : elements)
                    {
                        std::cerr << k << " ";
                    }
                    std::cerr << "\n\n";

                    std::cerr << "THEIR (Regular) Component root " << componentId << ": ";
                    for (int k : groupedTheirs[componentId])
                    {
                        std::cerr << k << " ";
                    }
                    std::cerr << "\n\n\n";

                    return false;
                }


            }

            return true;
        }

        bool areEqual(DisjointSet<int> &ds)
        {
            // Group elements with the same root
            std::map<int, std::set<int>> groupedOurs = this->groupComponents();
            std::map<int, std::set<int>> groupedTheirs = ds.groupComponents();

            // Make sure there is an equal number of roots
            if (groupedOurs.size() != groupedTheirs.size())
            {
                std::cerr << "Size of preimage graphs is not equal!.";
                return false;
            }

            for (const auto &[root, elements] : groupedOurs)
            {
                const int firstElement = *elements.begin();

                // Find which component it's in
                const int componentId = ds.findElement(firstElement);

                if (elements != groupedTheirs[componentId])
                {
                    std::cerr << "Corresponding components are not equal!\n\n";
                    std::cerr << "Our Component root             " << root << ": ";
                    for (int k : elements)
                    {
                        std::cerr << k << " ";
                    }
                    std::cerr << "\n\n";

                    std::cerr << "THEIR (Regular) Component root " << componentId << ": ";
                    for (int k : groupedTheirs[componentId])
                    {
                        std::cerr << k << " ";
                    }
                    std::cerr << "\n\n\n";

                    return false;
                }


            }

            return true;
        }

        std::map<int, std::set<int>> groupComponents()
        {
            std::map<int, std::set<int>> grouped;

            for (const auto &[key, value] : componentRoot)
            {
                grouped[value].insert(key);
            }

            return grouped;
        }

        void printByRoot()
        {
            const std::map<int, std::set<int>> grouped = groupComponents();

            // Now grouped[v] contains all keys with that value v
            for (const auto &[value, keys] : grouped)
            {
                std::cout << "Component root " << value << ": ";
                for (int k : keys)
                    std::cout << k << " ";
                std::cout << "\n";
            }

        }

        std::vector<int> getNeighbours(const int &triangleId, const TetMesh &tetMesh)
        {
            // Get the two neighbours of the current triangle in the fiber component
            std::vector<int> neighbours;

            // For all neighbours in the mesh
            for (const int &neighbourTriangleId : tetMesh.tetIncidentTriangles[triangleId])
            {
                //std::cout << "Neighbour ID is " << neighbourTriangleId << std::endl;

                // Is this neighbour in the fiber component?
                if (componentRoot.contains(neighbourTriangleId))
                {
                    neighbours.push_back(neighbourTriangleId);
                }
            }

            if (neighbours.size() == 0 || neighbours.size() > 2)
            {
                throw std::runtime_error( "Triangle has " + std::to_string(neighbours.size()) + " neighbours in its preimage graph.");
            }

            return neighbours;
        }


        void bfsSearch(const int &root, std::set<int> &visited, const TetMesh &tetMesh)
        {
            std::queue<int> q;
            q.push(root);

            visited.insert(root);
            this->componentRoot[root] = root;

            //std::cout << "\nConnected component with root " << root << "...\n";

            while (false == q.empty())
            {
                const int currentTriangleId = q.front();
                q.pop();

                for (const int &neighbourTriangleId : this->getNeighbours(currentTriangleId, tetMesh))
                {
                    if (false == visited.contains(neighbourTriangleId))
                    {
                        //printf("%d -> from %d.\n", currentTriangleId, neighbourTriangleId);
                        visited.insert(neighbourTriangleId);
                        q.push(neighbourTriangleId);
                        this->componentRoot[neighbourTriangleId] = root;
                    }
                }
            }
        }


        // Compute Connected Components from scratch
        void computeConnectedComponents(const TetMesh &tetMesh)
        {
            // Reset all the roots
            for (auto &[triangleId, rootId] : this->componentRoot)
            {
                rootId = triangleId;
            }

            // @TODO You can move this to the bfsSearch function
            std::set<int> visited;
            for (const auto &[triangleId, rootId] : this->componentRoot)
            {
                if (false == visited.contains(triangleId)) 
                {  
                    this->bfsSearch(triangleId, visited, tetMesh);
                }
            }
        }

        void updateConnectedComponents(const TetMesh &tetMesh, const std::vector<int> &minusTriangles, const std::vector<int> &plusTriangles, const PreimageGraph &preimageGraphPrevious)
        {
            this->componentRoot = preimageGraphPrevious.componentRoot;

            for (auto &triangleId : minusTriangles)
            {
                this->componentRoot.erase(triangleId);
            }

            // Add the plus triangles with a sentinel value
            for (auto &triangleId : plusTriangles)
            {
                this->componentRoot[triangleId] = triangleId;
            }

            // @TODO You can move this to the bfsSearch function
            std::set<int> visited;
            for (auto &triangleId : plusTriangles)
            {
                if (false == visited.contains(triangleId)) 
                {  
                    this->bfsSearch(triangleId, visited, tetMesh);
                }
            }

        }

};
