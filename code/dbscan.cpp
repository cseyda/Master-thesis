#include <vector>
#include <string>

//#include "node.h"
#include "dbscan.h"
#include "clustering.h"

using namespace std;

DBSCAN::DBSCAN(Met &metric) :
        met(&metric), 
        node_count(met->getCount()),
        cluster(node_count, 0),
        visited(node_count, false) {
};

unsigned int DBSCAN::run(Clustering &c){
    unsigned int clusterID = 1;
    vector<unsigned int> neighbours;
    double cout_limit = 0.0;
    
    cout << "db node_count: " << node_count << endl;
    for(unsigned int nodeID = 0; nodeID < node_count; nodeID++){
        
        if (double(nodeID) / node_count >= cout_limit){
            cout << double(nodeID) / node_count * 100 << " %" << endl;
            cout_limit+=.05;
        }
        
        //cout << "h" << endl;
        if (!visited[nodeID]){
            visited[nodeID] = true;
            neighbours.clear();
            //cout << "i" << endl;
            if (met->getNeighbors(nodeID, neighbours)){
                //cout << "j" << endl;
                cluster[nodeID] = clusterID;
                //visited[nodeID] = true;
                expand(neighbours, clusterID);
                clusterID++;
            }
        }
    }
    
    c.resize(clusterID);
    
    met->expand(c, cluster);
    
    return clusterID;
};

void DBSCAN::expand(vector<unsigned int>& neighbours, const unsigned int clusterID){
    Neighbourhood.clear();
    
    for (auto nodeID : neighbours){
        if (!visited[nodeID]){
            cluster[nodeID] = clusterID;
            visited[nodeID] = true;
            
            Neighbourhood.push_back(nodeID);
        }
    }
    
    //single ask
    /*
    unsigned int nodeID, nodeID2;
    
    for (unsigned int i=0; i<Neighbourhood.size(); i++){
        neighbours.clear();
        nodeID = Neighbourhood[i];
        //cout << nodeID << endl;
        if (met->getNeighbors(nodeID, neighbours)){
            for (unsigned int j=0; j<neighbours.size(); j++){
                nodeID2 = neighbours[j];
                
                // just set true?
                if (!visited[nodeID2]){
                    visited[nodeID2] = true;
                }
                
                if (cluster[nodeID2] == 0){
                    cluster[nodeID2] = clusterID;
                    
                    Neighbourhood.push_back(nodeID2); //list implementation?
                }
            }
        }
    }*/
    
    //multi
    
    unsigned int slice_size = 10;
    
    while (!Neighbourhood.empty()){
        
        slice.clear();
        indices.clear();
        dists.clear();
        enoughPts.clear();
        
        while (!Neighbourhood.empty() && (slice.size()!=slice_size)){
            slice.push_back(Neighbourhood.back());
            Neighbourhood.pop_back();
        }
        
        indices.resize(slice.size());
        dists.resize(slice.size());
        enoughPts.resize(slice.size());
        
        //cout << "d" << endl;
        met->getNeighbors(slice, indices, dists, enoughPts, met->getEps());
        
        //cout << "e" << endl;
        for(unsigned int l=0; l<slice.size(); l++){
            if (enoughPts[l]){
                for (unsigned int nodeID2 : indices[l]){
                    // just set true?
                    if (!visited[nodeID2]){
                        visited[nodeID2] = true;
                    }
                    
                    if (cluster[nodeID2] == 0){
                        cluster[nodeID2] = clusterID;
                        
                        Neighbourhood.push_back(nodeID2);
                    }
                }
            }
        }
        //cout << "f" << endl;
    }
};
