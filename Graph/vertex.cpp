
#include <math.h>
#include <algorithm>
#include "vertex.h"
#include "graph.h"
#include "edge.h"
#include <map>
#include <iostream>
#include <set>

using std::set;
using std::map;
using std::cout;
using std::endl;

Vertex::Vertex()
{
    
}

Vertex::Vertex(unsigned int objectId, 
            double longitude, 
            double latitude)
{
    objectId_ = objectId;
    longitude_ = longitude;
    latitude_ = latitude;
}

bool Vertex::operator<(const Vertex &otherVertex) const
{
    return objectId_ < otherVertex.getObjectId();
}

bool Vertex::operator==(const Vertex &otherVertex) const
{
    return objectId_ == otherVertex.getObjectId();
}

bool Vertex::operator!=(const Vertex &otherVertex) const
{
    return objectId_ != otherVertex.getObjectId();
}

bool Vertex::strictWeightCompare(const Vertex &thisVertex, const Vertex &otherVertex)
{
    return thisVertex.getWeight() < otherVertex.getWeight(); 
}

bool Vertex::strictEstimateCompare(const Vertex &thisVertex, const Vertex &otherVertex)
{
    return thisVertex.getEstimate() < otherVertex.getEstimate();
}

double Vertex::getLongitude() const
{
    return this->longitude_;
}

double Vertex::getLatitude() const
{
    return this->latitude_;
}

unsigned int Vertex::getObjectId() const
{
    return this->objectId_;
}


double Vertex::convertToCartesianX(const double &meanLongitude, const double &meanLatitude) const
{
    return R0_*cos(angleToRad(meanLatitude))*angleToRad(longitude_ - meanLongitude);
}


double Vertex::convertToCartesianY(const double &meanLatitude) const
{
    double latitudeDifferenceInRad = angleToRad(latitude_ - meanLatitude);
    return R0_*log(
                tan((latitudeDifferenceInRad/2) + (M_PI/4))
                );
}


float Vertex::angleToRad(const double &angleInDegrees) const
{
    return angleInDegrees*M_PI/180;
}



/*Vertex Vertex::getVertexfromId(unsigned int objectId,
                                map<unsigned int, std::pair<double, double>> graphVerticesMap) const
{
    double longitude = graphVerticesMap[objectId].first;
    double latitude = graphVerticesMap[objectId].second;
    Vertex vertex(objectId, longitude, latitude);
    return vertex;
}*/

set<Vertex> &Vertex::getAdjacencyList()
{
    return adjacencyList_;
}

double Vertex::getWeight() const
{
    return weight_;
}

void Vertex::setWeight(const double &weight)
{
    weight_ = weight;
}

double Vertex::getEstimate() const
{
    return estimate_;
}

void Vertex::setEstimate(const double &estimate)
{
    estimate_ =  estimate;
}

void Vertex::printVertexData() const
{
    cout << "ObjectId = " << objectId_
         << "      " 
         << "Longitude = " << longitude_
         << "      "
         << "Latitude = " << latitude_
         <<endl;
}

/*void Vertex::visitAllVertices()
{
    // Basic BFS algorithm in pseudo C+ code
    Graph::bfs(uint32_t vstart) {
    container<uint32_t> active_queue;
    set<uint32_t> closed_set;
    // ID of the start vertex
    active_queue.push_end(vstart);
    do {
    // from the current vertex in the front of the queue
    // compute all vertices reachable in 1 step
    auto vcurrent = active_queue.pop_front();
    closed_set.add(vcurrent);
    for(vnext in adjacency_list of vcurrent) {
    if (vnext is in closed_set) {
    continue;
    }
    if (vnext is not already in active_queue) {
    active_queue.push_end(vnext);
    }
    }
    } while (active_queue.size() != 0)
    }
}*/
