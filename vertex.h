#pragma once
#include <vector>
#include <map>
#include "edge.h"
#include <set>
#include <float.h>

using std::set;
using std::map;
using std::vector;

struct Vertex{
    private:
        unsigned int objectId_;
        double longitude_; //deg
        double latitude_; //deg
        double cartesianX_;
        double cartesianY_;
        double weight_ = 0;
        double estimate_ = 0;
        set<Vertex> adjacencyList_;
        const static int R0_ = 6378137;

        float angleToRad(const double &angleInDegrees) const;

    public:
        Vertex();
        Vertex(unsigned int objectId, 
            double longitude, 
            double latitude);
        double getLongitude() const;
        double getLatitude() const;
        set<Vertex> &getAdjacencyList();
        double getWeight() const;
        void setWeight(const double &weight);
        unsigned int getObjectId() const;
        double getEstimate() const;
        void setEstimate(const double &estimate);
        void printVertexData() const;
        /*void setAdjacencyList(//unsigned int objectId,
                                map<unsigned int, Vertex> vertices,
                            //map<unsigned int, std::pair<double, double>> graphVerticesMap,
                           map<unsigned int, map<unsigned int, Edge>> edges);*/

        static double getLength(Vertex startVertex, Vertex endVertex);
        double convertToCartesianX(const double &meanLongitude, const double &meanLatitude) const;
        double convertToCartesianY(const double &meanLatitude) const;
        bool operator <(const Vertex &otherVertex) const;
        bool operator ==(const Vertex &otherVertex) const;
        bool operator !=(const Vertex &otherVertex) const;
        static bool strictWeightCompare(const Vertex &thisVertex, const Vertex &otherVertex);
        static bool strictEstimateCompare(const Vertex &thisVertex, const Vertex &otherVertex);
};
