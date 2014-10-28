#include <vector>
#include <cmath>
#include "graph.h"
#include "dbscan.h"
#include "node.h"
#include "metric.h"

#include <lemon/list_graph.h>
//#include <lemon/concepts/graph.h>
#include <lemon/concepts/digraph.h>
#include <lemon/maps.h>
#include <lemon/path.h>

Graph::Graph(vector<Node*> &vecNodes, CombinedMetric* metrics, GraphHandler &handler){
    //structure nodes
    vector<unsigned int> nodesFound;
    for(unsigned int i=0; i<metrics->getCount(); i++){
        //Node* pNode = vecNodes[i];
        
        metrics->Met::getNeighbors(i, nodesFound); //??? try other
        
        for(unsigned int j=0; j<nodesFound.size(); j++){
            handler.addStructureEdge(i, nodesFound[j]);
        }
        nodesFound.clear();
    }
    
    //attribute nodes
    vector<Met*>* metrics_vec = metrics->metrics();
    for (unsigned int i=0; i<metrics_vec->size(); i++){
        
        DBSCAN scan(*(metrics_vec->at(i)));
        vector<vector<unsigned int>> clustering;
        unsigned int clusters = scan.run(clustering);
        
        unsigned int attr = handler.newAttributes(clusters++);
        
        for(unsigned int cluster=0; cluster<clustering.size(); cluster++){
            for(unsigned int i=0; i<clustering[cluster].size(); i++){
                unsigned int node = clustering[cluster][i];
                handler.addAttributeEdge(node, attr + cluster);
            }
        }
        
        scan.clear_clustering();
    }
};







GraphHandler::GraphHandler(unsigned int elements) :
        attribute_nodes(0),
        structure_nodes(elements),
        labels(g),
        w_node(g),
        w_arc(g),
        outDeg(g),
        nodes(structure_nodes),
        normalNode(g, true){
    
    g.reserveNode(structure_nodes);
    for (unsigned int i=0; i<structure_nodes; i++){
        nodes[i] = g.addNode();
        w_node[nodes[i]] = 0;
    }
};
    
GraphHandler::~GraphHandler(){};
    
    
void GraphHandler::addStructureEdge(unsigned int i, unsigned int j){
    g.addArc(nodes[i], nodes[j]);
    g.addArc(nodes[j], nodes[i]);
};

void GraphHandler::addAttributeEdge(unsigned int i, unsigned int j){
    if ((i < structure_nodes) && (j >=structure_nodes)){
        labels[nodes[i]].push_back(j);
        g.addArc(nodes[i], nodes[j]);
        g.addArc(nodes[j], nodes[i]);
    }
};

double GraphHandler::distance(unsigned int source, unsigned int target, unsigned int length){
    Path atm;
    vector<Path> paths;
    
    double sum = 0.0;
    
    for (OutArcIt a(g, nodes[source]); a!=lemon::INVALID; ++a){
        Node n = g.target(a);
        atm.addBack(a);
        
        sum += distance_(n, nodes[target], length-1, atm, paths);
    }
    
    return sum;
};

// number: count of attribute nodes to add
// return: id of first added
unsigned int GraphHandler::newAttributes(unsigned int count){
    attributes.push_back(new vector<unsigned int>);
    for (unsigned int i=0; i<count; i++){
        Node n = g.addNode();
        nodes.push_back(n);
        normalNode[n] = false;
        w_node[n] = attributes.size();
        attributes[attributes.size()-1]->push_back(attribute_nodes+i);
    }
    attribute_nodes+=count;
    
    return attribute_nodes-count;
};

//init all Edges
void GraphHandler::initEdges(){
    for(unsigned int i = 0; i < attributes.size(); i++){
        w.push_back(1.0);
    }
    
    updateEdges();
};

// update w weights
void GraphHandler::updateWeights(vector<vector<unsigned int>> &clustering, vector<unsigned int> &medians){
    if (clustering.size() != medians.size()){
        cout << "error\n";
    }
    
    vector<unsigned int> votes(w.size(), 0);
    unsigned int votes_sum = 0;
    for (unsigned int i=0; i<w.size(); i++){
        for (unsigned int j=0; j<clustering.size(); j++){
            for (unsigned int k=0; k<clustering[j].size(); k++){
                votes[i] += vote(i, clustering[j][k], medians[j]);
            }
        }
        votes_sum += votes[i];
    }
    
    w_sum = 0;
    for (unsigned int i=0; i<w.size(); i++){
        w[i] = 0.5 * (w[i] + (double (votes[i] * w.size())) / double (votes_sum));
        w_sum += w[i];
    }
    
};

//update Attribute Edges
void GraphHandler::updateEdges(){
    for (unsigned int i=0; i<nodes.size(); i++){
        Node source = nodes[i];
        for (OutArcIt a(g, source); a!=lemon::INVALID; ++a){
            Node target = g.target(a);
            
            w_arc[a] = p(source, target);
        }
    }
};

//updateAll
void GraphHandler::updateAll(vector<vector<unsigned int>> &clustering, vector<unsigned int> &medians){
    updateWeights(clustering, medians);
    updateEdges();
};

//private:
double GraphHandler::p(Node i, Node j){
    //double p = 0.0;
    double w = w_node[j];
    
    double w_s = 0.0;
    if (normalNode[i]){
        w_s = w_sum;
    }
    
    return w / (outDeg[i] + w_s);
}

unsigned int GraphHandler::vote(unsigned int attr, unsigned int i, unsigned int j){
    unsigned int vote = 0;
    if (labels[nodes[i]][attr] == labels[nodes[j]][attr]){
        vote = 1;
    }
    return vote;
};



double GraphHandler::distance_(Node source, Node target, unsigned int length, Path& atm, vector<Path> &paths){
    //reached walking distance
    if (length == 0){
        return 0.0;
    }
    
    double c = 0.2;
    double sum = 0.0;
    
    // target found
    if (source == target){
        double w_path = 1.0;
        
        for(unsigned int i=0; i< (unsigned int) atm.length(); i++){
            w_path += w_arc[atm.nth(i)];
        }
        double reset = c * pow((1-c), length);
        w_path *= reset;
        
        sum += w_path;
    }
    
    //recursion
    for (OutArcIt a(g, source); a!=lemon::INVALID; ++a){
        Node n = g.target(a);
        atm.addBack(a);
        
        sum += distance_(n, target, length-1, atm, paths);
    }
    
    atm.eraseBack();
    
    return sum;
};
