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

        // Which triangles belong to which root
        std::unordered_map<int, int> componentRoot;

        // For each component save a representatives so that we can do correspondence
        std::unordered_map<int, int> componentRepresentative;

        // All the unique components in this preimage graphs (effectively the unique values in componentRoot)
        //std::unordered_set<int> uniqueComponentIds;

        // The number of components found so far
        inline static int componentCount = 0; // declaration inside class

        PreimageGraph() { }

        // Initialize vertices (triangles), the edges are implicit
        PreimageGraph(const std::unordered_map<int, int> &_componentRoot)  : componentRoot(_componentRoot) { }


        // Depricated
        // @TODO Refactor with a better clear
        void clear()
        {
            this->componentRoot  = std::unordered_map<int, int>();
            this->componentRepresentative  = std::unordered_map<int, int>();
            //this->uniqueComponentIds  = std::unordered_set<int>();
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
            std::set<int> keys;

            for (const auto &pair : componentRepresentative)
            {
                keys.insert(pair.first);
            }

            return std::vector<int>(keys.begin(), keys.end());
        }

        void unitTestGetUniqueComponents()
        {
            //std::set<int> componentIdsA;

            //for (const auto &[triangleId, componentId] : this->componentRoot)
            //{
                //componentIdsA.insert(componentId);
            //}

            //std::set<int> componentIdsB;
            //for (const auto &componentId : this->uniqueComponentIds)
            //{
                //componentIdsB.insert(componentId);

            //}

            //if (componentIdsA != componentIdsB)
            //{
                //throw std::runtime_error( "Unique componentIds not computed correctly.");
            //}
        }

        const std::unordered_map<int, std::set<int>> groupComponents() const
        {
            std::unordered_map<int, std::set<int>> grouped;

            for (const auto &[key, value] : componentRoot)
            {
                grouped[value].insert(key);
            }

            return grouped;
        }

        const void printByRoot() const
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

            this->componentRepresentative[componentId] = root;

            // This is a new component now
            //this->uniqueComponentIds.insert(componentId);

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




        void updateComponentsSingular(const TetMesh &tetMesh, const std::pair<int, bool> &intersectingSegment)
        {
            //this->componentRoot = preimageGraphPrevious.componentRoot;
            //this->uniqueComponentIds = preimageGraphPrevious.uniqueComponentIds;

            const std::vector<int> &minusTriangles = tetMesh.getMinusTriangles(intersectingSegment.first, intersectingSegment.second);
            const std::vector<int> &plusTriangles = tetMesh.getPlusTriangles(intersectingSegment.first, intersectingSegment.second);


            for (auto &triangleId : minusTriangles)
            {
                auto it = componentRoot.find(triangleId);

                if (it == componentRoot.end())
                {
                    throw std::runtime_error("Minus triangle not found in preimage graph for singular.");
                }

                // This connected component will no longer exist
                //this->uniqueComponentIds.erase(it->second);

                // This component is over, no need to keep its representative
                this->componentRepresentative.erase(it->second);

                // Remove it from the list
                this->componentRoot.erase(it);
            }

            // Add the plus triangles with a sentinel value
            for (auto &triangleId : plusTriangles)
            {
                assert(false == this->componentRoot.contains(triangleId) && "Plus triangle is already in the preimage graph.");
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

        void updateComponentsRegular(TetMesh &tetMesh, const std::vector<std::pair<int, bool>> &intersectingEdges)
        {
            //this->componentRoot = preimageGraphPrevious.componentRoot;
            //this->uniqueComponentIds = preimageGraphPrevious.uniqueComponentIds;

            for (const auto &[edgeId, isDirectionLowerToUpper] : intersectingEdges)
            {
                const std::vector<int> &minusTriangles = tetMesh.getMinusTriangles(edgeId, isDirectionLowerToUpper);
                const std::vector<int> &plusTriangles = tetMesh.getPlusTriangles(edgeId, isDirectionLowerToUpper);

                assert(minusTriangles.size() > 0 && "Minus triangles of a regular edge are empty.");
                assert(plusTriangles.size() > 0 && "Plus triangles of a regular edge are empty.");

                //if (this->componentRoot.contains(minusTriangles[0]) == false)
                //{
                    //this->printByRoot();
                    //printf("\nMinus triangles: ");
                    //for (auto x : minusTriangles)
                    //{
                        //std::cout << x << " ";
                    //}
                    //printf("\n");
                    //printf("\nPlus triangles: ");
                    //for (auto x : plusTriangles)
                    //{
                        //std::cout << x << " ";
                    //}
                    //printf("\n");
                    //fflush(stdout);
                    //return;
                //}

                assert(this->componentRoot.contains(minusTriangles[0]) && "First minus triangle is not in the preimage graph.");

                const int componentId = this->componentRoot.at(minusTriangles[0]);

                for (auto &triangleId : minusTriangles)
                {
                    auto it = componentRoot.find(triangleId);

                    #ifndef NDEBUG
                    if (it == componentRoot.end())
                    {
                        //this->printByRoot();
                        //printf("\nMinus triangles: ");
                        //for (auto x : minusTriangles)
                        //{
                            //std::cout << x << " ";
                        //}
                        //printf("\n");
                        //printf("\nPlus triangles: ");
                        //for (auto x : plusTriangles)
                        //{
                            //std::cout << x << " ";
                        //}
                        //printf("\n");
                        //fflush(stdout);
                        throw std::runtime_error("Minus triangle not found in preimage graph.");
                    }

                    if (it->second != componentId)
                    {
                        throw std::runtime_error( "Not all triangles are removed from the same root!");
                    }
                    #endif

                    this->componentRoot.erase(it);

                }

                // Make sure the component representative is a triangle in the preimage graph (in case we deleted it with minus triangles)
                this->componentRepresentative[componentId] = plusTriangles[0];

                //std::cout << "Here are the plus triangles : \n";
                for (auto &triangleId : plusTriangles)
                {
                    assert(false == this->componentRoot.contains(triangleId) && "Plus triangle is already in the preimage graph.");
                    this->componentRoot[triangleId] = componentId;
                }
            }
        }

        const std::vector<std::pair<int, int>> establishCorrespondence(const TetMesh &tetMesh, const std::pair<int, bool> &intersectingSegment, const PreimageGraph &pg2) const
        {
            const std::vector<int> &minusTriangles = tetMesh.getMinusTriangles(intersectingSegment.first, intersectingSegment.second);
            const std::vector<int> &plusTriangles = tetMesh.getPlusTriangles(intersectingSegment.first, intersectingSegment.second);

            std::unordered_set<int> affectedComponents;

            //printf("Affected roots...");
            // These are all affected components, all other have a 1-1 correspondence
            for (const int &triangle : minusTriangles)
            {
                affectedComponents.insert(this->componentRoot.at(triangle));
                //printf("%d ", this->componentRoot.at(triangle));
            }

            //printf("Affected roots twin...");
            std::unordered_set<int> twinAffectedComponents;
            for (const int &triangle : plusTriangles)
            {
                twinAffectedComponents.insert(pg2.componentRoot.at(triangle));
                //printf("%d ", pg2.componentRoot.at(triangle));
            }


            //printf("\n\nOur preimage graph...\n");
            //this->printByRoot();
            //printf("\nTheir preimage graph...\n");
            //pg2.printByRoot();


            std::vector<std::pair<int, int>> componentCorrespondence;

            // 1-1 with either a twist or open at the boundary
            if (twinAffectedComponents.size() == 1 && affectedComponents.size() == 1)
            {
                //printf("--------------------------------------------------------------------- SPECIAL CASE HAPPENING!\n");
                const int affectedComponentId = *affectedComponents.begin();
                const int twinAffectedComponentId = *twinAffectedComponents.begin();

                componentCorrespondence.push_back({affectedComponentId, twinAffectedComponentId});
            }


            for (const auto &[componentId, representativeTriangleId] : this->componentRepresentative)
            {
                // If this component is not affected
                if (false == affectedComponents.contains(componentId))
                {
                    //printf("FOUND CORRESPONDENCE - Component %d corresponds to component %d\n", componentId, pg2.componentRoot.at(representativeTriangleId));
                    componentCorrespondence.push_back({componentId, pg2.componentRoot.at(representativeTriangleId)});
                }
            }

            //printf("\n\n");

            return componentCorrespondence;
        }


            
};
