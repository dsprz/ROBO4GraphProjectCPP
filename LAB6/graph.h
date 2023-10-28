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
        //map<unsigned int, std::pair<double, double>> verticesMap_; //<vertexId, <longitude, latitude>>
        map<unsigned int, Vertex> allVertices_;                    //map<objectId, Vertex>
        map<pair<unsigned int, unsigned int>, Edge> allEdges_;    //map<<startId, endId>, Edge>
        void createEdge();
        Vertex createVertex(unsigned int objectId, double longitude, double latitude);
        //void createAllVertices();
        void createGraph();
        void storeVerticesAndEdges();
        void makeConnections();
        vector<Vertex> reconstructPath(map<unsigned int, unsigned int> &parents, 
                                        const Vertex &startVertex, 
                                        const Vertex &endVertex);
        

    public:
        Graph();
        Graph(string filename);
        map<unsigned int, map<unsigned int, double>> getEdgesMap() const;
        map<unsigned int, std::pair<double, double>> getVerticesMap() const;
        //vector<Vertex> getAllVertices() const;
        map<unsigned int, map<unsigned int, Edge>> getAllEdges() const;
        long double computeWeight(const Vertex &startVertex, const Vertex &endVertex);
        void bfs(const Vertex &startVertex, const Vertex &endVertex);
        void dijkstra(const Vertex &startVertex, const Vertex &endVertex);
        void astar(const Vertex &startVertex, const Vertex &endVertex);
        //void static bfs(uint32_t vstart);
};