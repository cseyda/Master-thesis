#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include <cmath>
#include "dbscan.h"

#include <lemon/list_graph.h>
//#include <lemon/concepts/graph.h>
#include <lemon/concepts/digraph.h>
#include <lemon/maps.h>
#include <lemon/path.h>

// build structure edges based on combined distances
//       attribute edges based on clustering results of individual clustering

// Structure Node -> 0,0
// Attribute Node -> ai, jth value
//


class GraphHandler {
public:
    typedef lemon::ListDigraph Graph;
    //typedef lemon::concepts::Digraph Digraph;
    typedef Graph::NodeMap<vector<unsigned int>> NodeMap;
    typedef Graph::NodeMap<unsigned int> NodeMap_uint;
    typedef Graph::ArcMap<double> ArcMap;
    typedef lemon::OutDegMap<Graph> OutDegMap;
    typedef lemon::SimplePath<Graph> Path;
    
    typedef Graph::OutArcIt OutArcIt;
    typedef Graph::Node Node;
    
    typedef lemon::IterableBoolMap<Graph, Node> BoolMap;
    
    GraphHandler(unsigned int elements);
    ~GraphHandler();
    
    
    void addStructureEdge(unsigned int i, unsigned int j);
    void addAttributeEdge(unsigned int i, unsigned int j);
    
    double distance(unsigned int source, unsigned int target, unsigned int length);
    
    // number: count of attribute nodes to add
    // return: id of first added
    unsigned int newAttributes(unsigned int count);
    
    //init all Edges
    void initEdges();
    
    // update w weights
    void updateWeights(vector<vector<unsigned int>> &clustering, vector<unsigned int> &medians);
    
    //update Attribute Edges
    void updateEdges();
    
    //updateAll
    void updateAll(vector<vector<unsigned int>> &clustering, vector<unsigned int> &medians);
    
private:
    double p(Node i, Node j);
    
    unsigned int vote(unsigned int attr, unsigned int i, unsigned int j);
    
    double distance_(Node source, Node target, unsigned int length, Path& atm, vector<Path> &paths);
    
    unsigned int attribute_nodes;
    unsigned int structure_nodes;
    vector<vector<unsigned int>*> attributes;
    vector<double> w;
    
    double w_sum; // sum of w without w_0
    
    Graph g;
    NodeMap labels; //different attributes of a node
    
    NodeMap_uint w_node; // returns the attribute-group a node represents, eg every structure node is 0, every attribute node from attribute 1 has 1, ...
    
    ArcMap w_arc; // edge weights
    
    OutDegMap outDeg; // how much outgoing arcs for a given node
    
    vector<Node> nodes; //vector with the nodes, stru + attr
    
    BoolMap normalNode; // is a node a structure node?
};


class Graph {
public:
    Graph(vector<Node*> &vecNodes, CombinedMetric* metrics, GraphHandler &handler);
};

#endif
