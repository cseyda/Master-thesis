#ifndef CLUSTERING_H
#define CLUSTERING_H

#include <vector>
#include <string>

using namespace std;

class Clustering{
public:
    Clustering();
    Clustering(const string filename);
    ~Clustering();
    
    //unsigned int count();
    //unsigned int elements();
    
    bool read_file(const string filename, unsigned int sample_make, unsigned int sample_skip);
    bool write_file(const string filename, bool noise);
    
    void resize(unsigned int clusters);
    void add(unsigned int clusterID, unsigned int nodeID);
    void setNoise(bool noise);
//private:
    // clustering[0] contains NOISE
    vector<vector<unsigned int>> clustering;
//private:
  //  bool noise;
};

unsigned int dr(unsigned int val, unsigned int scale);

void sampling(Clustering& c, unsigned int min_size, unsigned int scale);

#endif
