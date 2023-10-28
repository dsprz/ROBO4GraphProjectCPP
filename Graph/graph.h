#pragma once
#include <string>
#include <vector>
#include <map>
#include "vertex.h"
#include "edge.h"

using std::vector;
using std::string;
using std::map;
using std::pair;

struct Graph{
    private:
        double longitudeOfCenter_ = 0; //degrees
        double latitudeOfCenter_ = 0; //degrees
        
        string filename_;
        map<unsigned int, Vertex> allVertices_;                    //map<objectId, Vertex>
        map<pair<unsigned int, unsigned int>, Edge> allEdges_;    //map<<startId, endId>, Edge>

        void storeVerticesAndEdges();
        void makeConnections();
        vector<Vertex> reconstructPath(const map<Vertex, Vertex> &parents, 
                                        const Vertex &startVertex, 
                                        const Vertex &endVertex) const;
        
    public:
        Graph();
        Graph(string filename);
        double computeWeight(const Vertex &startVertex, const Vertex &endVertex);
        void bfs(const Vertex &startVertex, const Vertex &endVertex);
        void dijkstra(const Vertex &startVertex, const Vertex &endVertex);
        void astar(const Vertex &startVertex, const Vertex &endVertex);
        //void static bfs(uint32_t vstart);
};
