#include <iostream>
#include <ostream>
#include <string>
#include <fstream>
#include "graph.h"
#include "vertex.h"
#include <map>
#include <vector>
#include <sstream>

using std::cout;
using std::endl;
using std::vector;
using std::string;

int main(int argc, char* argv[])
{

    std::string filename{argv[1]};
    Graph graph(filename);
    //cout << "lol" << endl;

    //int i = 0;
   // graph.getAllVertices()[1].getAdjacencyList()[0].printVertexData();
   //std::map<std::pair<int,int>, double> dmap;
   //dmap[{2,3}] = 3;
   //cout << dmap[{2,3}] << endl;
   //Vertex startVertex(193, -77.0121568832,38.9565958412);
   //Vertex endVertex(194, -77.0125834045,38.956595438);
   //cout << graph.getLength(startVertex, endVertex) << endl;
    return 0;
}
