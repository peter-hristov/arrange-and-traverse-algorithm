#pragma once

#include <cassert>
#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <unordered_map>

#include "./Hashing.h"

template<typename ElementType>
class DisjointSet {

    public:

        std::vector<int> parent;
        std::vector<int> rank;

        // Map a triangle to an ID (element ID)
        //std::map<std::set<int>, int> data;
        
        // Map from the data value to it's index in the data structure
        //std::unordered_map<ElementType, int, MyHash<ElementType>> data;
        std::unordered_map<ElementType, int, MyHash<ElementType>> data;
        //std::map<ElementType, int> data;


        DisjointSet() 
        {

        }

        DisjointSet(const std::set<ElementType> &preimageGraph) 
        {
            initialize(preimageGraph);
        }

        bool isEmpty()
        {
            return this->parent.size() == 0;
        }

        void clear()
        {
            this->parent = std::vector<int>();
            this->rank = std::vector<int>();
            this->data  = std::unordered_map<ElementType, int, MyHash<ElementType>>();
            //this->data  = std::map<ElementType, int>();
        }

        // Make sure everyone is poiting to the root
        void finalise()
        {
            for (int i = 0 ; i < parent.size() ; i++)
            {
                this->findIndex(i);
            }
        }

        std::vector<std::pair<int, ElementType>> getComponentRepresentatives()
        {
            std::set<int> uniqueComponetIds;
            std::vector<std::pair<int, ElementType>> componentRepresentatives;

            for (const auto &[element, index] : this->data)
            {
                const int componentId = this->findIndex(index);

                if (false == uniqueComponetIds.contains(componentId))
                {
                    componentRepresentatives.push_back({componentId, element});
                    uniqueComponetIds.insert(componentId);
                }
            }

            return componentRepresentatives;
        }

        //std::vector<int> getUniqueRepresentatives()
        //{

            //std::set<int> uniqueRoots;
            //std::vector<ElementType> representatives;

            //for (const auto &[key, value] : this->data)
            //{
                //const int root = this->findTriangle(key);

                //if (false == uniqueRoots.contains(root))
                //{
                    //representatives.push_back(key);
                    //uniqueRoots.insert(root);
                //}
            //}

            ////for(int i  = 0 ; i < parent.size() ; i++)
            ////{
                ////if (uniqueRoots.contains())
                ////uniqueRoots.insert(find(i));
            ////}

            //return representatives;
        //}

        //std::vector<int> getUniqueIndices()
        //{
            //std::set<int> uniqueRoots;

            //for(int i  = 0 ; i < parent.size() ; i++)
            //{
                //uniqueRoots.insert(findIndex(i));
            //}

            //return std::vector<int>(uniqueRoots.begin(), uniqueRoots.end());
        //}

        void initialize(const std::set<ElementType> &preimageGraph) 
        {
            int n = preimageGraph.size();

            rank.resize(n, 0);

            parent.resize(n);
            for (int i = 0; i < n; ++i) 
            {
                parent[i] = i;
            }

            // Map each triangle to an ID in the disjoint set
            int counter = 0;
            for(const ElementType triangle: preimageGraph)
            {
                data[triangle] = counter++;
            }

        }

        // Equals the number of distinct roots
        int countConnectedComponents()
        {
            std::set<int> roots;

            for(int i  = 0 ; i < parent.size() ; i++)
            {
                roots.insert(findIndex(i));
            }

            return roots.size();
        }

        void addElement(const ElementType &triangle1)
        {
            const int newElementId = parent.size();

            this->parent.push_back(newElementId);
            this->rank.push_back(0);
            this->data[triangle1] = newElementId;
        }

        // Interface for triangles
        int findElement(const ElementType &triangle)
        {
            assert(data.contains(triangle));
            return findIndex(data[triangle]);
        }

        void unionElements(const ElementType &triangle1, const ElementType &triangle2) 
        {
            assert(data.contains(triangle1));
            assert(data.contains(triangle2));
            unionIndices(data[triangle1], data[triangle2]);

        }

        bool connectedElements(const ElementType &triangle1, const ElementType &triangle2)
        {
            assert(data.contains(triangle1));
            assert(data.contains(triangle2));
            return connected(data[triangle1], data[triangle2]);
        }

        // Find with path compression
        int findIndex(const int &x) {
            if (parent[x] != x) {
                // Path compression
                parent[x] = findIndex(parent[x]);  
            }
            return parent[x];
        }

        // Union by rank
        void unionIndices(const int &x, const int &y) 
        {
            const int rootX = findIndex(x);
            const int rootY = findIndex(y);

            if (rootX != rootY) 
            {
                if (rank[rootX] > rank[rootY]) 
                {
                    parent[rootY] = rootX;
                } 
                else if (rank[rootX] < rank[rootY]) 
                {
                    parent[rootX] = rootY;
                } 
                else 
                {
                    parent[rootY] = rootX;
                    rank[rootX]++;
                }
            }
        }

        // Check if two nodes are connected
        bool connectedIndex(const int &x, const int &y) {
            return findIndex(x) == findIndex(y);
        }

        void print(const std::function<void(const ElementType&)> &printElement)
        {
            const int connectedComponents = this->countConnectedComponents();

            std::cout << "This preimage graph has " << connectedComponents << " connected components." << std::endl;


            std::map<int, std::vector<ElementType>> components;

            for (const auto &[element, index] : this->data)
            {
                const int componentId = this->findIndex(index);
                components[componentId].push_back(element);
            }


            for (const auto &[componentId, elements] : components)
            {
                std::cout << "Component " << componentId << " has elements:\n";
                for (const auto &element : elements)
                {
                    printElement(element);
                }
                std::cout << "\n\n";
            }
        }



        const std::unordered_map<int, std::set<int>> groupComponents()
        {
            std::unordered_map<int, std::set<int>> grouped;

            for (const auto &[key, value] : this->data)
            {
                grouped[this->findIndex(value)].insert(key);
            }

            return grouped;

        }


        void printByRoot()
        {
            std::unordered_map<int, std::set<int>> grouped = this->groupComponents();

            // Now grouped[v] contains all keys with that value v
            for (const auto &[value, keys] : grouped)
            {
                std::cout << "Component root " << value << ": ";
                for (int k : keys)
                {
                    std::cout << k << " ";
                }
                std::cout << "\n";
            }
        }
};


template<typename ElementType>
bool operator==(const DisjointSet<ElementType>& lhs, const DisjointSet<ElementType>& rhs)
{
    if (lhs.parent.size() != rhs.parent.size())
        return false;

    if (lhs.rank.size() != rhs.rank.size())
        return false;

    // Compare data maps
    if (lhs.data.size() != rhs.data.size())
        return false;

    for (const auto& [element, parentId] : lhs.data)
    {
        auto it = rhs.data.find(element);

        if (it == rhs.data.end() || it->second != parentId)
        {
            return false;
        }
    }

    return true;
}

template<typename ElementType>
bool operator!=(const DisjointSet<ElementType>& lhs, const DisjointSet<ElementType>& rhs)
{
    return !(lhs == rhs);
}
