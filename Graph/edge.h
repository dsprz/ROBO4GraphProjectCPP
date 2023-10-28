#pragma once
//#include "vertex.h"

struct Edge{
    private:
        double shapelen_;
        unsigned int startVertexId_;
        unsigned int endVertexId_;
    public:
        Edge();
        Edge(const unsigned int &startVertexId, 
            const unsigned int &endVertexId, 
            const double &shapelen);
        //long double computeWeight(Vertex startVertex, Vertex endVertex);
        unsigned int getStartVertexId() const;
        unsigned int getEndVertexId() const;
        double getEdgeLength() const;
        bool operator==(const Edge &otherEdge) const;

};
