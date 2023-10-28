#pragma once
//#include "vertex.h"

struct Edge{
    private:
        double shapelen_;
        unsigned int startVertexId_;
        unsigned int endVertexId_;
    public:
        Edge();
        Edge(unsigned int startVertexId, unsigned int endVertexId, double shapelen);
        //long double computeWeight(Vertex startVertex, Vertex endVertex);
        unsigned int getStartVertexId() const;
        unsigned int getEndVertexId() const;
        double getEdgeLength() const;
        bool operator==(const Edge &otherEdge) const;

};