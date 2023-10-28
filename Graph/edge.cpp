#include "edge.h"
#include <iostream>
#include <math.h>
#include "vertex.h"

Edge::Edge()
{

}
Edge::Edge(const unsigned int &startVertexId, 
    const unsigned int &endVertexId, const double &shapelen)
{   
    shapelen_ = shapelen;
    startVertexId_ = startVertexId;
    endVertexId_ = endVertexId;
}

bool Edge::operator==(const Edge &otherEdge) const
{
    return shapelen_ == otherEdge.getEdgeLength() 
            && startVertexId_ == otherEdge.getStartVertexId()
            && endVertexId_ == otherEdge.getEndVertexId();
}

unsigned int Edge::getStartVertexId() const
{
    return startVertexId_;
}

unsigned int Edge::getEndVertexId() const
{
    return endVertexId_;
}

double Edge::getEdgeLength() const
{
    if(startVertexId_ != endVertexId_)
    {
        return shapelen_;
    }
    return 0;
}

/*long double Edge::computeWeight(Vertex startVertex, Vertex endVertex)
{
    const int R0 = 6378137;

    double cartesianStartVertexX = Vertex::convertToCartesianX(startVertex.getCartesianX());
    double startVertexY = Vertex::convertToCartesianY(startVertex.getCartesianY());
    double endVertexX = Vertex::convertToCartesianX(endVertex.getCartesianX());
    double endVertexY = Vertex::convertToCartesianY(endVertex.getCartesianY());
    //euclidian distance between startVertex and endVertex
    //return sqrt(startVertex.get)
    return 0;
}*/
