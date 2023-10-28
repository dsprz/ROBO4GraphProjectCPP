#include "graph.h"
#include "vertex.h"
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <math.h>
#include <deque>
#include <set>
#include <algorithm>

using std::vector;
using std::string;
using std::deque;
using std::set;
using std::cout;
using std::endl;

Graph::Graph()
{

}


Graph::Graph(std::string filename)
{
    filename_ = filename;
    storeVerticesAndEdges();
    makeConnections();
    //std::cout << latitudeOfCenter_ << std::endl;
    //std::cout << longitudeOfCenter_ << std::endl;
    /*map<unsigned int, Vertex>::iterator it;
    for(it=allVertices_.begin(); it!=allVertices_.end(); ++it)
    {
        vector<Vertex> vert = it->second.getAdjacencyList();
        for(Vertex v : vert)
        {
            v.printVertexData();
        }
    }*/
    unsigned int vid = 86771;
    cout << "adjacent Vertices for Vertex " << vid <<" : { " << endl;; 
    for(const Vertex &v : allVertices_[vid].getAdjacencyList())
    {
        v.printVertexData();
    }
    cout << "} " << endl;
    //cout << computeWeight(allVertices_[86771], allVertices_[24774]) << endl;
    cout << "#####################################################################################" << endl;
    //cout << computeWeight(allVertices_[200], allVertices_[197]) << endl;    
    //bfs(allVertices_[86771], allVertices_[110636]);
    //dijkstra(allVertices_[86771], allVertices_[110636]);
    //astar(allVertices_[86771], allVertices_[110636]);
    astar(allVertices_[86771], allVertices_[110636]);
    //dijkstra(allVertices_[86771], allVertices_[110636]);


    //dijkstra(allVertices_[42159], allVertices_[287]);
    //astar(allVertices_[42159], allVertices_[287]);

    //cout << "From 17792 to end passing by 17796 : " << 
      //      computeWeight(allVertices_.at(17792), allVertices_.at(17796)) << endl; 
            //+ computeWeight(allVertices_.at(17796), allVertices_.at(110636)) << endl;

    //cout << "From 17792 to end passing by 16696 : " <<  
      //   computeWeight(allVertices_.at(17792), allVertices_.at(16696)) << endl; 
           // +  computeWeight(allVertices_.at(16696), allVertices_.at(110636))<< endl;
    
    //cout << allVertices_.at(16696).convertToCartesianX(longitudeOfCenter_, latitudeOfCenter_) << endl;

    /*cout << "Edge length from 17792 to end = " << computeWeight(allVertices_.at(17792), allVertices_.at(110636))<< endl;//allEdges_.at({17792, 110636}).getEdgeLength() << endl;
    cout << "Edge length from 16696 to end = " << computeWeight(allVertices_.at(16696), allVertices_.at(110636))<< endl;*/
    //cout << allVertices_[24774].getWeight() << endl;
    /*map<pair<unsigned int, unsigned int>, Edge>::iterator it;
    unsigned int startVertexId;
    unsigned int endVertexId;
    Edge edge;
    /*for(it = allEdges_.begin(); it != allEdges_.end(); ++it)
    {
        startVertexId = it->first.first;
        endVertexId = it->first.second;
        edge = it->second;

    }*/
    //createAllVertices();
    /*for(string s : edges_[0])
    {
        std::cout << s << std::endl;
    }*/
}

void Graph::storeVerticesAndEdges()
{
    /*
    Stores edges data in the edges_ vector
    Stores vertices data in the map verticesMap_ as follows : verticesMap_[vertexID] = pair(longitude, latitude)
    */
    vector<string> row;
    string line, word;
    std::ifstream file(filename_, std::ios::in);

    unsigned int vertexId;
    double vertexLongitude;
    double vertexLatitude;
    unsigned int startVertexId;
    unsigned int endVertexId;
    double length;

    if (file.is_open())
    {
        while (getline(file, line))
        {
            row.clear();
            std::stringstream str(line);
            while(getline(str, word, ','))
            {
                row.push_back(word);
            }

            if (row[0] == "V")
            {
                vertexId = stod(row[1]);
                vertexLongitude = stod(row[2]);
                vertexLatitude = stod(row[3]);
                longitudeOfCenter_+=vertexLongitude;
                latitudeOfCenter_+=vertexLatitude;
                //verticesMap_[vertexId] = std::make_pair(vertexLongitude, vertexLatitude);
        
                Vertex vertex(vertexId, vertexLongitude, vertexLatitude);
                allVertices_[vertexId] = vertex;
                //vertices_.push_back(row);
            }
            else if (row[0] == "E")
            {
                startVertexId = stod(row[1]);
                endVertexId = stod(row[2]);
                length = stod(row[3]);

                //edgesMap_[startVertexId][endVertexId] = length;
                Edge edge(startVertexId, endVertexId, length);
                allEdges_[{startVertexId,endVertexId}] = edge;
                //edgesMap_[endVertexId].insert(std::pair<unsigned int, double>(startVertexId, length));
 
            }
        }
        longitudeOfCenter_ /= allVertices_.size();
        latitudeOfCenter_ /= allVertices_.size();
    }
}

void Graph::makeConnections()
{
    map<pair<unsigned int, unsigned int>, Edge>::iterator it;

    unsigned int startVertexId;
    unsigned int endVertexId;
    for(it = allEdges_.begin(); it!=allEdges_.end(); ++it)
    {
        startVertexId = it->first.first;
        endVertexId = it->first.second;
        allVertices_.at(startVertexId).getAdjacencyList().insert(allVertices_.at(endVertexId));
        allVertices_.at(endVertexId).getAdjacencyList().insert(allVertices_.at(startVertexId));
    }
}

void Graph::bfs(const Vertex &startVertex, const Vertex &endVertex)
{
    //one way problem
    //vector<Vertex> activeQueue;
    deque<Vertex> activeQueue;
    set<Vertex> alreadyVisitedVertices;
    map<Vertex, Vertex> parents; //map<parent, child>, map[parent] = child, keep track of parents to reconstruct the path
    Vertex currentVertex;
    bool found = false;

    activeQueue.emplace_back(startVertex); 
    alreadyVisitedVertices.insert(startVertex);

    while(!activeQueue.empty() && !found)
    {
        //why doesn't currentVertex = activeQueue.front() work properly ?
        currentVertex = allVertices_.at(activeQueue.front().getObjectId());
        //cout << "currentVertex = " ;
        //currentVertex.printVertexData();
        activeQueue.pop_front();
        alreadyVisitedVertices.insert(currentVertex);

        for(Vertex adjacentVertex : currentVertex.getAdjacencyList())
        {
            //adjacentVertex.printVertexData();
            auto adjacentVertexId = adjacentVertex.getObjectId();
            auto currentVertexId = currentVertex.getObjectId();
            try
            {
                auto edgeStartId = allEdges_.at({currentVertexId, adjacentVertexId}).getStartVertexId();
                if (alreadyVisitedVertices.find(adjacentVertex) != alreadyVisitedVertices.end()
                || currentVertexId != edgeStartId)
                //closed set
                 //If Vertex already visited
                { 
                    continue;
                }
            }
            catch(const std::out_of_range& e)
            {
                continue;
            }

            if(std::find(activeQueue.begin(), activeQueue.end(), adjacentVertex) == activeQueue.end())
            {                
                activeQueue.push_back(adjacentVertex);
                parents[adjacentVertex] =  currentVertex;
                if(adjacentVertex == endVertex)
                {
                    found = true;
                    cout << "found ! " << endl;
                    break;
                }
                
            }
        }
        //cout << activeQueue.size() << endl;

    }

    /*while(!found)
    {
        for (Vertex nextVertex : activeQueue)
        {
            alreadyVisited[nextVertex.getObjectId()] = nextVertex;
            //nextVertex.printVertexData();
            activeQueue.pop_front();
            for(Vertex adjacentVertex : nextVertex.getAdjacencyList())
            {
                adjacentVertex.printVertexData();
                if (alreadyVisited.find(adjacentVertex.getObjectId()) == alreadyVisited.end())
                {
                    activeQueue.push_back(adjacentVertex);
                    cout << "je suis ajouté à activeQueue" << endl;

                    //parents[adjacentVertex.getObjectId()] = nextVertex.getObjectId();
                }

                else cout << "je suis déjà visité" << endl;

                if(adjacentVertex.getObjectId() == endVertexId)
                {
                    found = true;
                    cout << "found ! " << endl;
                    cout << adjacentVertex.getObjectId() << endl;
                    cout << endVertexId << endl;
                    break;
                }
                cout << activeQueue.size() << endl;

            }
        }
        cout << "end" << endl;*/

    
    unsigned int numberOfVertices = alreadyVisitedVertices.size();
    cout << "Number of vertices visited = " << numberOfVertices << endl;
    vector<Vertex> VerticesOnPath = reconstructPath(parents, startVertex, endVertex);
    int i = 1;
    double length = 0;
    //auto childId = startVertex.getObjectId();
    unsigned int childVertexId;
    unsigned int parentVertexId = startVertex.getObjectId();
    for(auto it = VerticesOnPath.rbegin(); it != VerticesOnPath.rend(); ++it)
    {   
        childVertexId = it->getObjectId();
        if (parentVertexId!=childVertexId)
        {
            length += allEdges_.at({parentVertexId, childVertexId}).getEdgeLength();
        }        
        parentVertexId = childVertexId;
        cout << "Vertex["<< i << "] = " << childVertexId << ", length = " << length <<endl;
        i++;  
    }
    cout << "Total vertices on path from start to end = " << VerticesOnPath.size() << endl;

}

void Graph::dijkstra(const Vertex &startVertex, const Vertex &endVertex)
{
    cout << "dijkstra" << endl;
    deque<Vertex> activeQueue;
    set<Vertex> alreadyVisitedVertices; //closed set
    map<Vertex, Vertex> parents;
    activeQueue.push_back(startVertex);
    //bool found = false;
    Vertex currentVertex;
    unsigned int currentVertexId;
    unsigned int adjacentVertexId;
    double sumOfWeights;
    
    while(!activeQueue.empty())
    {
        currentVertex = allVertices_.at(activeQueue.front().getObjectId());
        activeQueue.pop_front();
        //currentVertex.setWeight(FLT_MAX);

        if(currentVertex == endVertex)
            { cout << "coucou" << endl;
            break;}

        alreadyVisitedVertices.insert(currentVertex);

        for(Vertex v : currentVertex.getAdjacencyList())
        {
            //cout <<"size = " << adjacentVertex.getAdjacencyList().size() << endl;
            
            adjacentVertexId =  v.getObjectId();
            currentVertexId = currentVertex.getObjectId();
            Vertex &adjacentVertex = allVertices_.at(adjacentVertexId);

            try
            {
                auto edgeStartId = allEdges_.at({currentVertexId, adjacentVertexId}).getStartVertexId();
                if (alreadyVisitedVertices.find(adjacentVertex) != alreadyVisitedVertices.end()
                || currentVertexId != edgeStartId)
                //closed set
                 //If Vertex already visited
                { 
                    continue;
                }
            }
            catch(const std::out_of_range& e)
            {
                continue;
            }

            sumOfWeights = currentVertex.getWeight() + computeWeight(currentVertex, adjacentVertex);
            //sum of all Weights up until currentVertex; 

            if(std::find(activeQueue.begin(), activeQueue.end(), adjacentVertex) == activeQueue.end())
            // if adjacentVertex is not already in activeQueue
            {
                adjacentVertex.setWeight(sumOfWeights);
                activeQueue.emplace_back(adjacentVertex);
                parents[adjacentVertex] = currentVertex;
            }
            else if(sumOfWeights < adjacentVertex.getWeight())
            //if adjacentVertex is already in activeQueue and the sum of weights is lesser than its current weight
            {
                adjacentVertex.setWeight(sumOfWeights);
                parents[adjacentVertex] = currentVertex;
            }
        //cout << "weight : " << adjacentVertex.getWeight() << endl;
        }
        std::sort(activeQueue.begin(), activeQueue.end(), Vertex::strictWeightCompare);
    
    }
    unsigned int numberOfVertices = alreadyVisitedVertices.size();
    cout << "Number of vertices visited = " << numberOfVertices << endl;
    vector<Vertex> VerticesOnPath = reconstructPath(parents, startVertex, endVertex);

    int i = 1;
    double length = 0;
    //auto childId = startVertex.getObjectId();
    unsigned int childVertexId;
    unsigned int parentVertexId = startVertex.getObjectId();
    for(auto it = VerticesOnPath.rbegin(); it != VerticesOnPath.rend(); ++it)
    {   
        childVertexId = it->getObjectId();
        if (parentVertexId!=childVertexId)
        {
            length += allEdges_.at({parentVertexId, childVertexId}).getEdgeLength();
        }
        parentVertexId = childVertexId;
        cout << "Vertex["<< i << "] = " << childVertexId << ", length = " << length <<endl;
        i++;  
    }
    cout << "Total vertices on path from start to end = " << VerticesOnPath.size() << endl;
    cout << "length = " << length << endl;
}

void Graph::astar(const Vertex &startVertex, const Vertex &endVertex)
{
    cout << "astar" << endl;
    deque<Vertex> activeQueue;
    set<Vertex> alreadyVisitedVertices; //closed set
    map<Vertex, Vertex> parents;
    activeQueue.emplace_back(startVertex);
    unsigned int currentVertexId;
    unsigned int adjacentVertexId;

    
    while(!activeQueue.empty())
    {
        Vertex &currentVertex = allVertices_.at(activeQueue.front().getObjectId());
        activeQueue.pop_front();
        /*cout << "currentVertexId = " << currentVertex.getObjectId() <<
        " estimate = " << currentVertex.getEstimate()
        << " weight = " << currentVertex.getWeight()
        << endl;*/

        if(currentVertex == endVertex)
        { 
            break;
        }

        alreadyVisitedVertices.insert(currentVertex);

        for(Vertex v : currentVertex.getAdjacencyList())
        {

            //cout <<"size = " << adjacentVertex.getAdjacencyList().size() << endl;
            adjacentVertexId = v.getObjectId();
            currentVertexId = currentVertex.getObjectId();
            Vertex &adjacentVertex = allVertices_.at(adjacentVertexId);

            try
            {
                auto edgeStartId = allEdges_.at({currentVertexId, adjacentVertexId}).getStartVertexId();
                if (alreadyVisitedVertices.find(adjacentVertex) != alreadyVisitedVertices.end() || currentVertexId != edgeStartId)
                //closed set
                 //If Vertex already visited
                { 
                    continue;
                }
            }
            catch(const std::out_of_range& e)
            {
                continue;
            }

            auto g = currentVertex.getWeight() + allEdges_.at({currentVertexId, adjacentVertexId}).getEdgeLength();
            auto f = g + computeWeight(adjacentVertex, endVertex); //heuristic distance = computeWeight 

            if(std::find(activeQueue.begin(), activeQueue.end(), adjacentVertex) == activeQueue.end() || f < adjacentVertex.getEstimate())
            // if adjacentVertex is not already in activeQueue
            {
                adjacentVertex.setWeight(g);
                adjacentVertex.setEstimate(f);
                parents[adjacentVertex] = currentVertex;
                activeQueue.emplace_back(adjacentVertex);
            }
            /*else if(f < adjacentVertex.getEstimate()) 
            //if adjacentVertex is already in activeQueue
            {
                cout << "ALREADY IN QUEUE" << endl;
                cout << "EMPLACING BACK : " << adjacentVertexId << endl;
                cout << "OLD WEIGHT : " << adjacentVertex.getWeight() << " "
                 << "OLD ESTIMATE : " << adjacentVertex.getEstimate() << endl;
                adjacentVertex.setWeight(g);
                adjacentVertex.setEstimate(f);
                cout << "NEW WEIGHT : " << adjacentVertex.getWeight() << " " << 
                "NEW ESTIMATE : " << adjacentVertex.getEstimate() << endl;
                parents[adjacentVertex] = currentVertex;
                activeQueue.emplace_back(adjacentVertex);

            }*/
        }
        std::sort(activeQueue.begin(), activeQueue.end(), Vertex::strictEstimateCompare);
    }

    unsigned int numberOfVertices = alreadyVisitedVertices.size();
    cout << "Number of vertices visited = " << numberOfVertices << endl;
    vector<Vertex> VerticesOnPath = reconstructPath(parents, startVertex, endVertex);

    int i = 1;
    double length = 0;
    //auto childId = startVertex.getObjectId();
    unsigned int childVertexId;
    unsigned int parentVertexId = startVertex.getObjectId();
    for(auto it = VerticesOnPath.rbegin(); it != VerticesOnPath.rend(); ++it)
    {   
        childVertexId = it->getObjectId();
        if(parentVertexId != childVertexId)
        {
            length += allEdges_.at({parentVertexId, childVertexId}).getEdgeLength();
        }
        parentVertexId = childVertexId;
        cout << "Vertex["<< i << "] = " << childVertexId << ", length = " << length << endl;
        i++;  
    }
    cout << "Total vertices on path from start to end = " << VerticesOnPath.size() << endl;
    cout << "length = " << length << endl;

}

vector<Vertex> Graph::reconstructPath(const map<Vertex, Vertex> &parents, 
                                    const Vertex &startVertex, 
                                    const Vertex &endVertex) const
{
    //Reconstructs the path from end to start

    unsigned int startVertexId = startVertex.getObjectId();
    unsigned int endVertexId =  endVertex.getObjectId();

    vector<Vertex> vertices;
    vertices.emplace_back(endVertex);
    Vertex parent = parents.at(endVertex);
    vertices.emplace_back(allVertices_.at(parent.getObjectId()));
    Vertex oldParent;
    while(parent != startVertex)
    {
        oldParent = parent;
        parent = parents.at(parent);
        vertices.emplace_back(allVertices_.at(parent.getObjectId()));

    }
    return vertices;  
}

double Graph::computeWeight(const Vertex &startVertex, const Vertex &endVertex)
{
    /*
    Returns the euclidian distance between 2 vertices if map<startVertexId, endVertexId> is not found.
    Otherwise, returns the edge length already written in the CSV
    */
    unsigned int startVertexId = startVertex.getObjectId();
    unsigned int endVertexId = endVertex.getObjectId();
    if (allEdges_.find({startVertexId, endVertexId}) == allEdges_.end())
    {
        double cartesianStartVertexX = startVertex.convertToCartesianX(longitudeOfCenter_, latitudeOfCenter_);
        double cartesianStartVertexY = startVertex.convertToCartesianY(latitudeOfCenter_);
        double cartesianEndVertexX = endVertex.convertToCartesianX(longitudeOfCenter_, latitudeOfCenter_);
        double cartesianEndVertexY = endVertex.convertToCartesianY(latitudeOfCenter_);

        //cout << "I am going to return the formula with the square root " << endl;
        return sqrt(pow(cartesianStartVertexX - cartesianEndVertexX, 2) 
                + pow(cartesianStartVertexY - cartesianEndVertexY, 2)
                    ); 
    }

    return allEdges_.at({startVertexId, endVertexId}).getEdgeLength();    
}


/*Vertex Graph::createVertex(unsigned int objectId, double longitude, double latitude)
{
    Vertex vertex(objectId, longitude, latitude);
    return vertex;
}*/

/*void Graph::createAllVertices()
{
    for(auto const& v : verticesMap_)
    {
        unsigned int objectId = v.first;
        double longitude = v.second.first;
        double latitude = v.second.second;
        Vertex vertex(objectId, longitude, latitude);
        vertex.setAdjacencyList(verticesMap_, edgesMap_); //BEAUCOUP TROP LONGUE
        allVertices_.emplace_back(vertex);
    }
}
*/

/*double Graph::getLength(Vertex startVertex, Vertex endVertex)
{
    //std::make_pair(startVertex.getObjectId(), endVertex.getObjectId());
    return edgesMap_.at({startVertex.getObjectId(), endVertex.getObjectId()});

    //return edgesMap_[{startVertex.getObjectId(), endVertex.getObjectId()}];
}
*/
/*map<unsigned int, map<unsigned int, double>> Graph::getEdgesMap() const
{
    return edgesMap_;
}

map<unsigned int, std::pair<double, double>> Graph::getVerticesMap() const
{
    return verticesMap_;
}

map<unsigned int, map<unsigned int, Edge>> Graph::getAllEdges() const
{
    return allEdges_;
}
*/
/*vector<Vertex> Graph::getAllVertices() const
{
    return allVertices_;
}
*/
/*vector<vector<string>> Graph::getAllEdges() const
{
    return edges_;
}*/

/*static void Graph::bfs(uint32_t startVertexId, uint32_t endVertexId)
{
    container<uint32_t> active_queue;
    set<uint32_t> closed_set;
    // ID of the start vertex
    active_queue.push_end(startVertexId);
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
}*/

    // Basic BFS algorithm in pseudo C+ code
    //Graph::bfs(uint32_t vstart) {
//Implement a variant of the above BFS algorithm with a start and end vertices. 
//The do loop must stop when the end vertex is found. You have to find the best STL container for the active queue.
    //}
