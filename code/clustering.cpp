#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <ctime>        // std::time
#include <cstdlib>      // std::rand, std::srand
#include <algorithm>    // std::shuffle
#include "clustering.h"

using namespace std;

Clustering::Clustering(){};
Clustering::Clustering(const string filename){
    read_file(filename, 1, 0);
};

Clustering::~Clustering(){};

//unsigned int Clustering::count(){};
//unsigned int Clustering::elements(){};

void Clustering::resize(unsigned int clusters){
    clustering.resize(clusters);
};

void Clustering::add(unsigned int clusterID, unsigned int nodeID){
    clustering[clusterID].push_back(nodeID);
};
    
bool Clustering::read_file(const string filename, unsigned int sample_make, unsigned int sample_skip){
    ifstream is;
    is.open(filename);
    
    if (is.is_open()){
        //vector<unsigned int> sample_count;
        string line;
        unsigned int node;
        unsigned int cluster;

        while(getline(is, line)){
            if (line != "") {
                istringstream iss(line);
                //index cluster
                while (iss >> node >> cluster) {
                    //if (sample_count.size() < (cluster+1)){
                    //    sample_count.resize(cluster+1);
                    //}
                    
                    //if (!(sample_count[cluster] > sample_make)){
                        
                        if(clustering.size() < (cluster+1)){
                            clustering.resize(cluster+1);
                        }
                        clustering[cluster].push_back(node);
                    
                    //}
                    //sample_count[cluster] = (sample_count[cluster]+1)%(sample_make+sample_skip);
                }
            }
        }
        is.close();
        /*
        unsigned int min_size = 8;
        // remove clusters too small
        //cout << "Clusterings " << clustering.size() << endl;
        for (auto it=clustering.begin(); it!=clustering.end(); it++){
            if (it->size() < min_size){
                //cout << "Erased " << it->size() << " nodes\n";
                clustering.erase(it);
                it--;
            }
        }
        //cout << "Clusterings " << clustering.size() << endl;
        */
        return true;
    } else {
        return false;
    }
};


bool Clustering::write_file(const string filename, bool noise){
    ofstream myfile;
    myfile.open(filename);
    
    unsigned int clusterID = 0;
    
    if (noise) {
        clusterID = 1;
    }
    
    for(; clusterID < clustering.size(); clusterID++){
        for(unsigned int nodeID : clustering[clusterID]){
            myfile << nodeID << " " << clusterID << endl;
            //vid, cid = line.split(",")#vertex, component ID
        }
    }
    
    myfile.close();
    
    return true;
};

unsigned int dr(unsigned int val, unsigned int scale) {
    double mult = val / (double)scale;
    double trinum = (sqrt(8.0 * mult + 1.0) - 1.0) / 2.0;
    return (unsigned int) trinum * scale;
};

void sampling(Clustering& c, unsigned int min_size, unsigned int scale){
    //if (skip_first != 0){
    std::srand( unsigned(std::time(0)));
    
    for (auto &cluster : c.clustering){
        if (cluster.size() < min_size){
            cluster.clear();
        } else {
            std::random_shuffle(cluster.begin(), cluster.end());
        }
    }
    
    //remove now empty clusters
    c.clustering.erase(std::remove_if(c.clustering.begin(), c.clustering.end(), [](vector<unsigned int>& cl){return cl.empty();}), c.clustering.end());
    
    //remove points based on DR
    for (auto &cluster : c.clustering){
        unsigned int i = dr(cluster.size(), min_size);
        
        if (i<cluster.size()){
            cluster.erase(cluster.begin()+i, cluster.end());
        }
    }
};
