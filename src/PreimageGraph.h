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
        // For each vertex in the preimage graph, store its associated triangleId
        std::set<int> triangleIds;
        std::map<int, int> componentRoot;

        PreimageGraph() { }

        // Initialize vertices (triangles), the edges are implicit
        PreimageGraph(const std::set<int> &_vertexTriangleId) : triangleIds(_vertexTriangleId) { }

        // @TODO Refactor with a better clear
        void clear()
        {
            this->triangleIds  = std::set<int>();
            this->componentRoot  = std::map<int, int>();
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
                if (triangleIds.contains(neighbourTriangleId))
                {
                    neighbours.push_back(neighbourTriangleId);
                }
            }

            assert(neighbours.size() == 1 || neighbours.size() == 2);

            return neighbours;
        }

        bool areEqual(DisjointSet<int> &ds)
        {
            std::map<int, std::set<int>> groupedOurs = this->groupComponents();
            std::map<int, std::set<int>> groupedTheirs = ds.groupComponents();

            assert(groupedOurs.size() == groupedTheirs.size());

            for (const auto &[root, elements] : groupedOurs)
            {
                const int firstElement = *elements.begin();

                // Find which component it's in
                const int componentId = ds.findElement(firstElement);

                if (elements != groupedTheirs[componentId])
                {
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
            std::map<int, std::set<int>> grouped = groupComponents();

            // Now grouped[v] contains all keys with that value v
            for (const auto &[value, keys] : grouped)
            {
                std::cout << "Component root " << value << ": ";
                for (int k : keys)
                    std::cout << k << " ";
                std::cout << "\n";
            }

        }

        void bfsSearch(const int &root, std::set<int> &visited, const TetMesh &tetMesh)
        {
            std::queue<int> q;
            q.push(root);

            visited.insert(root);
            this->componentRoot[root] = root;

            while (false == q.empty())
            {
                const int currentTriangleId = q.front();
                q.pop();

                for (const int &neighbourTriangleId : this->getNeighbours(currentTriangleId, tetMesh))
                {
                    if (false == visited.contains(neighbourTriangleId))
                    {
                        visited.insert(neighbourTriangleId);
                        q.push(neighbourTriangleId);
                        this->componentRoot[neighbourTriangleId] = root;
                    }
                }
            }
        }


        // Compute Connected Components
        void computeConnectedComponents(const TetMesh &tetMesh)
        {


            // Print all the triangles
            //std::cout << "Here are the preimage graph triangles" << std::endl;
            //for (const int &triangleId : this->triangleIds)
            //{
                //std::cout << triangleId << std::endl;

            //}



            std::set<int> visited;
            for (const int &triangleId : this->triangleIds)
            {

                //std::cout << "Currently at triangle ID : " << triangleId << std::endl;


                if (visited.contains(triangleId))
                {
                    continue;
                }

                this->bfsSearch(triangleId, visited, tetMesh);
            }


            //std::cout << std::endl << std::endl;
        }

        // Get connected components

        // Find (which cc am I in)
};
