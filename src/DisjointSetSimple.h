#pragma once

#include <cassert>
#include <iostream>
#include <vector>
#include <set>


class DisjointSetSimple {

    public:

        std::vector<int> parent;
        std::vector<int> rank;

        DisjointSetSimple() 
        {

        }

        DisjointSetSimple(const int numberOfElements)
        {
            rank.assign(numberOfElements, 0);

            parent.resize(numberOfElements);
            std::iota(parent.begin(), parent.end(), 0);
        }

        bool isEmpty()
        {
            return this->parent.size() == 0;
        }

        // Make sure everyone is poiting to the root
        void finalise()
        {
            for (int i = 0 ; i < parent.size() ; i++)
            {
                this->find(i);
            }
        }

        //std::vector<std::pair<int, ElementType>> getComponentRepresentatives()
        //{
            //std::set<int> uniqueComponetIds;
            //std::vector<std::pair<int, ElementType>> componentRepresentatives;

            //for (const auto &[element, index] : this->data)
            //{
                //const int componentId = this->findIndex(index);

                //if (false == uniqueComponetIds.contains(componentId))
                //{
                    //componentRepresentatives.push_back({componentId, element});
                    //uniqueComponetIds.insert(componentId);
                //}
            //}

            //return componentRepresentatives;
        //}

        int countComponents()
        {
            std::set<int> roots;

            for(int i  = 0 ; i < parent.size() ; i++)
            {
                roots.insert(find(i));
            }

            return roots.size();
        }

        void add(const int numberOfElements)
        {
            const size_t oldSize = parent.size();

            parent.resize(oldSize + numberOfElements);
            std::iota(parent.begin() + oldSize, parent.end(), static_cast<int>(oldSize));

            rank.insert(rank.end(), numberOfElements, 0);
        }

        // Find with path compression
        int find(const int &x) {
            if (parent[x] != x) {
                // Path compression
                parent[x] = find(parent[x]);  
            }
            return parent[x];
        }

        // Union by rank
        void unify(const int &x, const int &y) 
        {
            const int rootX = find(x);
            const int rootY = find(y);

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
        bool connected(const int &x, const int &y) 
        {
            return find(x) == find(y);
        }

        //void print(const std::function<void(const ElementType&)> &printElement)
        //{
            //const int connectedComponents = this->countConnectedComponents();

            //std::cout << "This preimage graph has " << connectedComponents << " connected components." << std::endl;


            //std::map<int, std::vector<ElementType>> components;

            //for (const auto &[element, index] : this->data)
            //{
                //const int componentId = this->findIndex(index);
                //components[componentId].push_back(element);
            //}


            //for (const auto &[componentId, elements] : components)
            //{
                //std::cout << "Component " << componentId << " has elements:\n";
                //for (const auto &element : elements)
                //{
                    //printElement(element);
                //}
                //std::cout << "\n\n";
            //}
        //}



        //const std::unordered_map<int, std::set<int>> groupComponents()
        //{
            //std::unordered_map<int, std::set<int>> grouped;

            //for (const auto &[key, value] : this->data)
            //{
                //grouped[this->findIndex(value)].insert(key);
            //}

            //return grouped;

        //}


        //void printByRoot()
        //{
            //std::unordered_map<int, std::set<int>> grouped = this->groupComponents();

            //// Now grouped[v] contains all keys with that value v
            //for (const auto &[value, keys] : grouped)
            //{
                //std::cout << "Component root " << value << ": ";
                //for (int k : keys)
                //{
                    //std::cout << k << " ";
                //}
                //std::cout << "\n";
            //}
        //}
};
