#pragma once

#include <cassert>
#include <unordered_map>
#include <unordered_set>
#include <vector>
//#include <set>
#include <queue>
#include <limits> 

#include "./TetMesh.h"
#include "./DisjointSet.h"


class PreimageGraph
{
    public:
        std::unordered_map<int, int> componentRoot;
        std::unordered_set<int> uniqueComponentIds;

        inline static int componentCount = 0; // declaration inside class

        PreimageGraph() { }

        // Initialize vertices (triangles), the edges are implicit
        PreimageGraph(const std::unordered_map<int, int> &_componentRoot)  : componentRoot(_componentRoot) { }


        // @TODO Refactor with a better clear
        void clear()
        {
            this->componentRoot  = std::unordered_map<int, int>();
            this->uniqueComponentIds  = std::unordered_set<int>();
        }

        bool areEqual(PreimageGraph &p)
        {
            // Group elements with the same root
            std::unordered_map<int, std::set<int>> groupedOurs = this->groupComponents();
            std::unordered_map<int, std::set<int>> groupedTheirs = p.groupComponents();

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
            std::unordered_map<int, std::set<int>> groupedOurs = this->groupComponents();
            std::unordered_map<int, std::set<int>> groupedTheirs = ds.groupComponents();

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

        std::vector<int> getUniqueComponents()
        {
            return std::vector<int>(this->uniqueComponentIds.begin(), this->uniqueComponentIds.end());
        }

        void unitTestGetUniqueComponents()
        {
            std::set<int> componentIdsA;

            for (const auto &[triangleId, componentId] : this->componentRoot)
            {
                componentIdsA.insert(componentId);
            }

            std::set<int> componentIdsB;
            for (const auto &componentId : this->uniqueComponentIds)
            {
                componentIdsB.insert(componentId);

            }

            if (componentIdsA != componentIdsB)
            {
                throw std::runtime_error( "Unique componentIds not computed correctly.");
            }
        }

        std::unordered_map<int, std::set<int>> groupComponents()
        {
            std::unordered_map<int, std::set<int>> grouped;

            for (const auto &[key, value] : componentRoot)
            {
                grouped[value].insert(key);
            }

            return grouped;
        }

        void printByRoot()
        {
            const std::unordered_map<int, std::set<int>> grouped = groupComponents();

            // Now grouped[v] contains all keys with that value v
            for (const auto &[value, keys] : grouped)
            {
                std::cout << "Component root " << value << ": ";
                for (int k : keys)
                    std::cout << k << " ";
                std::cout << "\n";
            }

        }

        // Get the two neighbours in the mesh that are part of the fiber
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


        void bfsSearch(const int &root, std::unordered_set<int> &visited, const TetMesh &tetMesh)
        {
            std::queue<int> q;
            q.push(root);

            visited.insert(root);

            // Obtain the next available componentId
            const int componentId = PreimageGraph::componentCount++;

            // All new components will have this as their root
            this->componentRoot[root] = componentId;

            // This is a new component now
            this->uniqueComponentIds.insert(componentId);

            //std::cout << "\nConnected component with root " << componentId << "...\n";
            //std::cerr << "\nComponent FAST INSIDE: ";
            //for (int c : this->getUniqueComponentsFase())
            //{
                //std::cerr << c << " ";
            //}

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
                        this->componentRoot[neighbourTriangleId] = componentId;
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
            std::unordered_set<int> visited;
            for (const auto &[triangleId, rootId] : this->componentRoot)
            {
                if (false == visited.contains(triangleId)) 
                {  
                    this->bfsSearch(triangleId, visited, tetMesh);
                }
            }
        }




        void updateConnectedComponents(const TetMesh &tetMesh, const std::pair<int, bool> &intersectingSegment, const PreimageGraph &preimageGraphPrevious)
        {
            this->componentRoot = preimageGraphPrevious.componentRoot;
            this->uniqueComponentIds = preimageGraphPrevious.uniqueComponentIds;

            const std::vector<int> &minusTriangles = tetMesh.getMinusTriangles(intersectingSegment.first, intersectingSegment.second);
            const std::vector<int> &plusTriangles = tetMesh.getPlusTriangles(intersectingSegment.first, intersectingSegment.second);

            for (auto &triangleId : minusTriangles)
            {
                auto it = componentRoot.find(triangleId);

                if (it == componentRoot.end())
                {
                    throw std::runtime_error("Minus triangle not found in preimage graph.");
                }

                // This connected component will no longer exist
                uniqueComponentIds.erase(it->second);

                // Remove it from the list
                componentRoot.erase(it);
            }

            // Add the plus triangles with a sentinel value
            for (auto &triangleId : plusTriangles)
            {
                this->componentRoot[triangleId] = triangleId;
            }

            // @TODO You can move this to the bfsSearch function
            std::unordered_set<int> visited;
            for (auto &triangleId : plusTriangles)
            {
                if (false == visited.contains(triangleId)) 
                {  
                    this->bfsSearch(triangleId, visited, tetMesh);
                }
            }
        }



        void updateConnectedComponentsEdge2(const TetMesh &tetMesh, const std::vector<std::vector<int>> &minusTriangles, const std::vector<std::vector<int>> &plusTriangles, const PreimageGraph &preimageGraphPrevious)
        {
            //std::cout << "-------------------------------------------------------- GOING INTO A NEW COMPONENT!!!\n";

            this->componentRoot = preimageGraphPrevious.componentRoot;
            this->uniqueComponentIds = preimageGraphPrevious.uniqueComponentIds;

            //std::cout << "\nHere are the INITIAL roots and components : " << std::endl;
            //for (const auto &[triangle, root] : this->componentRoot)
            //{
                //printf("triangle ID = %d, rootId = %d\n", triangle, root);
            //}

            //for (int i = 0 ; i < minusTriangles.size() ; i++)
            //{
                //printf("\nThe minus triangles are : ");
                //for (int t : minusTriangles[i])
                //{
                    //std::cout << t << " ";
                //}
                //printf("\nThe plus triangles are : ");
                //for (int t : plusTriangles[i])
                //{
                    //std::cout << t << " ";
                //}
            //}
            //std::cout << "\n\n";


            //std::cout << "\n\nNow recomputing properly...\n";

            for (int i = 0 ; i < minusTriangles.size() ; i++)
            {

                if (minusTriangles[i].empty())
                {
                    throw std::runtime_error("minusTriangles[i] is empty");
                }

                //std::cout << "\nHere are the roots and components : " << std::endl;
                //for (const auto &[triangle, root] : this->componentRoot)
                //{
                    //printf("triangle ID = %d, rootId = %d\n", triangle, root);
                //}


                //std::cout << "The first minus triangle is " << minusTriangles[i][0] << std::endl;

                const int rootId = this->componentRoot.at(minusTriangles[i][0]);

                //std::cout << "Here are the minus triangles : \n";
                for (auto &triangleId : minusTriangles[i])
                {

                    auto it = componentRoot.find(triangleId);

                    if (it == componentRoot.end())
                    {
                        throw std::runtime_error("Minus triangle now found in preimage graph.");
                    }

                    //printf("triangle ID = %d\n", triangleId);
                    // @TODO Double call
                    if (it->second != rootId)
                    {
                        throw std::runtime_error( "Not all triangles are removed from the same root!");
                    }

                    this->componentRoot.erase(it);

                }

                //std::cout << "Here are the plus triangles : \n";
                for (auto &triangleId : plusTriangles[i])
                {
                    //printf("triangle ID = %d\n", triangleId);
                    this->componentRoot[triangleId] = rootId;
                }
            }
        }

        void updateConnectedComponentsEdge3(TetMesh &tetMesh, const std::vector<std::pair<int, bool>> &intersectingEdges, const PreimageGraph &preimageGraphPrevious)
        {
            this->componentRoot = preimageGraphPrevious.componentRoot;
            this->uniqueComponentIds = preimageGraphPrevious.uniqueComponentIds;


            //std::cout << "\nHere are the INITIAL roots and components : " << std::endl;
            //for (const auto &[triangle, root] : this->componentRoot)
            //{
                //printf("triangle ID = %d, rootId = %d\n", triangle, root);
            //}



            //std::cout << "\n\nNow recomputing properly...\n";



            for (const auto &[edgeId, isDirectionLowerToUpper] : intersectingEdges)
            {
                const std::vector<int> &minusTriangles = tetMesh.getMinusTriangles(edgeId, isDirectionLowerToUpper);
                const std::vector<int> &plusTriangles = tetMesh.getPlusTriangles(edgeId, isDirectionLowerToUpper);

                //std::cout << edgeId << " - " << isDirectionLowerToUpper << std::endl;


                //for (int i = 0 ; i < minusTriangles.size() ; i++)
                //{
                    //printf("\nThe minus triangles are : ");
                    //for (int t : minusTriangles)
                    //{
                        //std::cout << t << " ";
                    //}
                    //printf("\nThe plus triangles are : ");
                    //for (int t : plusTriangles)
                    //{
                        //std::cout << t << " ";
                    //}
                //}
                //std::cout << "\n\n";


                //if (minusTriangles.empty() || plusTriangles.empty())
                //{
                    //throw std::runtime_error("MinusTriangles or plus triangles is empty");
                //}

                const int rootId = this->componentRoot.at(minusTriangles[0]);

                //std::cout << "Here are the minus triangles : \n";
                for (auto &triangleId : minusTriangles)
                {

                    auto it = componentRoot.find(triangleId);

                    if (it == componentRoot.end())
                    {
                        throw std::runtime_error("Minus triangle not found in preimage graph.");
                    }

                    //printf("triangle ID = %d\n", triangleId);
                    // @TODO Double call
                    if (it->second != rootId)
                    {
                        throw std::runtime_error( "Not all triangles are removed from the same root!");
                    }

                    this->componentRoot.erase(it);

                }

                //std::cout << "Here are the plus triangles : \n";
                for (auto &triangleId : plusTriangles)
                {
                    //printf("triangle ID = %d\n", triangleId);
                    this->componentRoot[triangleId] = rootId;
                }
            }

        }

};
