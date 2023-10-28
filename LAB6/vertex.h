#pragma once
#include <vector>
#include <map>
#include "edge.h"
#include <set>

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
        long double weight_ = 0;
        long double estimate_ = 0;
        set<Vertex> adjacencyList_;
        const static int R0_ = 6378137;

        //void visitAllVertices();

        void setCartesianX(double cartesianX);
        void setCartesianY(double cartesianY);
        bool isAdjacentTo(Vertex vertex, map<std::pair<unsigned int, unsigned int>, double>) const;

        double computeMapLongitudeCenter() const; //Computes lamdaCenter
        double computeMapLatitudeCenter() const; //Computes phiCenter
        float angleToRad(double angleInDegrees) const;
        vector<unsigned int> extract_keys(map<unsigned int, Edge> const& input_map) const;

    public:
        Vertex();
        Vertex(unsigned int objectId, double longitude, double latitude);
        Vertex getVertexfromId(unsigned int objectId,
                            map<unsigned int, std::pair<double, double>> graphVerticesMap) const;
        /*double getCartesianX() const;
        double getCartesianY() const;*/
        double getLongitude() const;
        double getLatitude() const;
        set<Vertex> &getAdjacencyList();
        long double getWeight() const;
        void setWeight(long double weight);
        unsigned int getObjectId() const;
        long double getEstimate() const;
        void setEstimate(long double estimate);
        void printVertexData() const;
        /*void setAdjacencyList(//unsigned int objectId,
                                map<unsigned int, Vertex> vertices,
                            //map<unsigned int, std::pair<double, double>> graphVerticesMap,
                           map<unsigned int, map<unsigned int, Edge>> edges);*/

        static double getLength(Vertex startVertex, Vertex endVertex);
        double convertToCartesianX(const double &meanLongitude, const double &meanLatitude) const;
        double convertToCartesianY(const double &meanLongitude, const double &meanLatitude) const;
        bool operator <(const Vertex &otherVertex) const;
        bool operator >(const Vertex &otherVertex) const;
        bool operator ==(const Vertex &otherVertex) const;
        bool operator >=(const Vertex &otherVertex) const;
        bool operator <=(const Vertex &otherVertex) const;
        static bool strictWeightCompare(const Vertex &thisVertex, const Vertex &otherVertex);
        static bool strictEstimateCompare(const Vertex &thisVertex, const Vertex &otherVertex);
};
