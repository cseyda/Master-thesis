#ifndef DBSCAN_H
#define DBSCAN_H

#include <vector>
#include <string>

#include "node.h"
#include "metric.h"
#include "clustering.h"

using namespace std;

class DBSCAN {
public:
    DBSCAN(Met &met);
    
    //return a vector representing the different clusters
    // cluster 0 is NOISE
    unsigned int run(Clustering &clustering);

private:
    Met *met;

    const unsigned int node_count;
    
    vector<unsigned int> cluster;
    vector<bool> visited;
    vector<unsigned int> Neighbourhood; //used in expand
    
    void expand(vector<unsigned int>& rgp, unsigned int nCluster);
    
    vector<vector<unsigned int>> indices;
    vector<vector<double>> dists;
    
    vector<unsigned int> slice;
    vector<bool> enoughPts;
};
#endif
