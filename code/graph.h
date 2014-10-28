#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include <cmath>
#include <iostream>
#include <utility>      // std::pair, std::make_pair
#include <unordered_map>
#include <queue>

#include "dbscan.h"
#include "node.h"

#include <lemon/list_graph.h>
//#include <lemon/concepts/graph.h>
#include <lemon/concepts/digraph.h>
#include <lemon/maps.h>
#include <lemon/path.h>

#include <lemon/lgf_reader.h>
#include <lemon/lgf_writer.h>

// build structure edges based on combined distances
//       attribute edges based on clustering results of individual clustering

// Structure Node -> 0,0
// Attribute Node -> ai, jth value
//

// Prototype
// GraphHandler
class GraphHandler {
public:
    GraphHandler(unsigned int elements);
    virtual ~GraphHandler();
    
    virtual void addEdge(unsigned int i, unsigned int j) = 0;
    virtual void addEdge(unsigned int i, unsigned int j, double weight) = 0;
    virtual void addNodeAttribute(unsigned int node, unsigned int attribute) = 0;
    virtual void addAttributeLimit(unsigned int limit) = 0;

    virtual void saveGraph(string location) = 0;
    virtual void loadGraph(string location) = 0;
protected:
    unsigned int node_count;
};

//
// LemonGraph
//
class LemonGraph : public GraphHandler {
public:
    typedef lemon::ListDigraph Graph;
    typedef Graph::NodeMap<vector<unsigned int>> NodeAttributeMap;
    typedef Graph::NodeMap<pair<unsigned int, unsigned int>> NodeMap_pair;
    typedef Graph::ArcMap<double> ArcMap;
    typedef lemon::OutDegMap<Graph> OutDegMap;
    typedef lemon::SimplePath<Graph> Path;
    
    
    typedef Graph::OutArcIt OutArcIt;
    typedef Graph::NodeIt NodeIt;
    typedef Graph::Node Node;
    typedef Graph::Arc Arc;
    typedef lemon::IdMap<Graph, Node> IdMap;
    
    LemonGraph(unsigned int node_count);
    ~LemonGraph();
    
    void addEdge(unsigned int i, unsigned int j);
    void addEdge(unsigned int i, unsigned int j, double weight);

    void addNodeAttribute(unsigned int node, unsigned int attribute);
    void addAttributeLimit(unsigned int limit);
    
    void saveGraph(string location);
    void loadGraph(string location);
    
    void initAugmentedGraph();

    void updateAll(vector<vector<unsigned int>> &clustering, vector<unsigned int> &medians);
    
    double distance(unsigned int source, unsigned int target, unsigned int length);
    
    bool getNeighbors(unsigned int pt, vector<unsigned int>& rgpNodesFound, vector<double>* distances, double eps, unsigned int minPts);
    
    unsigned int bfs(unsigned int source, unsigned int length);
    
private:
    double p(Node i, Node j);
    unsigned int vote(unsigned int attr, Node i, Node j);
    void updateEdges();
    void updateWeights(vector<vector<unsigned int>> &clustering, vector<unsigned int> &medians);
    double distance_(Arc a_now, Node target, unsigned int length, double atm, vector<double> &reset, set<Node> &visited, unordered_map<unsigned int, double> &hashed);
    
    
    Graph g;
    //vector with the nodes, stru + attr
    vector<Node> nodes; 
    //different attributes of a node in a vector
    NodeAttributeMap attributes;
    
    //
    //Algo
    //
    NodeMap_pair node_repr; // returns the attribute-group a node represents, eg every structure node is 0, every attribute node from attribute 1 has 1, ...
    
    ArcMap arc_weights; // edge weights
    
    OutDegMap node_degree; // how much outgoing arcs for a given node
    IdMap node_ids;
    
    vector<double> w;
    vector<unsigned int> attrLimits;
    
    
    vector<bool> visited;
    
    vector<pair<Node, unsigned int>> distance_vec;
    queue<pair<Node, unsigned int>> todo;
    
};



//
// SimpleGraph
//

class SimpleGraph : public GraphHandler {
public:
    SimpleGraph(unsigned int elements);
    ~SimpleGraph();
    
    void addEdge(unsigned int i, unsigned int j);
    void addEdge(unsigned int i, unsigned int j, double weight);
    void addNodeAttribute(unsigned int node, unsigned int attribute);
    void addAttributeLimit(unsigned int limit);

    void saveGraph(string location);
    void loadGraph(string location);
private:
    vector<vector<pair<unsigned int, double>>> edges;
};

//
// CreateGraph
//
class Graph {
public:
    Graph(CombinedMetric* metrics, GraphHandler &handler);
};

class Graph2 {
public:
    Graph2(DB2& db, vector<vector<pair<unsigned int, float>>> *scan, double w_1, double c, unsigned int dist_limit, unsigned int w_threshold, bool compressed);
};
#endif
