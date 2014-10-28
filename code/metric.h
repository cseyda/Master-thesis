#ifndef METRIC_H
#define METRIC_H

#include <vector>
#include <string>
//#include <tuple>
#include <utility>

#include "flann/flann.hpp"

#include "my_dists.h"
#include "simple_index.h"
#include "mongo.h"
#include "clustering.h"

using namespace std;

class Met {
public:
    Met(string name, string collection, double eps, unsigned int minPts);
    virtual ~Met();
    
    string getName();
    string getCollection();
    
    double getEps();
    void setEps(double eps);
    
    unsigned int getPoints();
    
    unsigned int getWeight();
    void setWeight(unsigned int weight);
    
    virtual bool getNeighbors(unsigned int pt, vector<unsigned int>& rgpNodesFound, vector<double>* distances, double eps, unsigned int minPts) = 0;
    
    virtual void getNeighbors(vector<unsigned int>& pts, vector<vector<unsigned int>>& indices, vector<vector<double>>& dists, vector<bool>& enoughPts, double eps) = 0;
    
    bool getNeighbors(unsigned int pt, vector<unsigned int>& rgpNodesFound);
    bool getNeighbors(unsigned int pt, vector<unsigned int>& rgpNodesFound, vector<double>& distances);
    
    virtual double distance(unsigned int a, unsigned int b) = 0;
    //virtual double distance(unsigned int a, unsigned int *b){return 0.0;};
    //virtual double distance(unsigned int *a, unsigned int *b){return 0.0;};

    virtual double getMin() = 0;
    virtual double getMax() = 0;
    
    virtual unsigned int getCount() = 0;
    
    virtual void expand(Clustering& c, vector<unsigned int>& cluster);
    
protected:
    unsigned int _weight;
    double       _eps;
    unsigned int _minPts;

private:
    const string _name;
    const string _collection;
};


//
// FLANN
//
class LocF : public Met {
public:
    typedef double DistType;
    typedef flann::L2_Simple<DistType> Dist;
    
    LocF(DB2& db, double eps, unsigned int minPts);
    LocF(DB2& db, double eps, unsigned int minPts, vector<unsigned int>& pts);
    
    ~LocF();
    
    void getNeighbors(vector<unsigned int>& pts, vector<vector<unsigned int>>& indices, vector<vector<double>>& dists, vector<bool>& enoughPts, double eps);
    
    bool getNeighbors(unsigned int pt, vector<unsigned int>& rgpNodesFound, vector<double>* distances, double eps, unsigned int minPts);

    double distance(unsigned int a, unsigned int b);
    double getMin();
    double getMax();
    unsigned int getCount();
    
private:
    flann::Index<Dist> tree;
    DB2* db;
    const unsigned int nn;
    flann::Matrix<DistType> dataset;
    Dist _distance;
    
    DistType* query_space;
    double min_lon;
    double max_lon;
    
    double min_lat;
    double max_lat;
    
    flann::SearchParams params;
    //std::vector< std::vector<int> > indices;
    //std::vector<std::vector<double> > dists;

};



//
// Jaccard Index
//
class JaccardIndex : public Met {
public:
    typedef unsigned int DistType;//char DistType;
    typedef Jaccard<DistType> Dist;

    JaccardIndex(DB2& db, double eps, unsigned int minPts, string kind, bool use_index=true);
    
    ~JaccardIndex();
    
    bool getNeighbors(unsigned int pt, vector<unsigned int>& rgpNodesFound, vector<double>* distances, double eps, unsigned int minPts);
    
    void getNeighbors(vector<unsigned int>& pts, vector<vector<unsigned int>>& indices, vector<vector<double>>& dists, vector<bool>& enoughPts, double eps);
    
    double distance(unsigned int a, unsigned int b);
    double distance(unsigned int a, unsigned int *b);
    double distance(unsigned int *a, unsigned int *b);
    
    double getMin();
    double getMax();
    unsigned int getCount();
    
private:
    DB2* db;
    flann::Matrix<DistType> dataset;
    Dist _distance;
    unsigned int nn; //number of elements per vector
    unsigned int bigram_ids; //number of bigrams
    SimpleIndex<Dist, flann::Matrix<DistType>, unsigned char> index;
    string kind;
};

//
// combined metric
//
class CombinedMetric : public Met {
public:
    CombinedMetric(DB2& db, vector<Met*> *metrics, double eps, unsigned int minPts);
    
    ~CombinedMetric();
    
    bool getNeighbors(unsigned int pt, vector<unsigned int>& rgpNodesFound, vector<double>* distances, double eps, unsigned int minPts);
    
    void getNeighbors(vector<unsigned int>& pts, vector<vector<unsigned int>>& indices, vector<vector<double>>& dists, vector<bool>& enoughPts, double eps);
    
    double distance(unsigned int a, unsigned int b);
    double getMin();
    double getMax();
    unsigned int getCount();
    
    vector<Met*>* metrics();
    
    void expand(Clustering& c, vector<unsigned int>& cluster);
    
private:
    vector<Met*>* _metrics;
    vector<double> _factor;
    DB2* db;
    vector<unsigned int> id;
    vector<vector<unsigned int>> trans_map;
};


//
// Graph
//
class GraphMetric : public Met {
public:

    GraphMetric(DB2& db, double eps, unsigned int minPts, double w_1, double c, unsigned int dist_limit, unsigned int w_threshold, bool compressed);
    
    ~GraphMetric();
    
    bool getNeighbors(unsigned int pt, vector<unsigned int>& rgpNodesFound, vector<double>* distances, double eps, unsigned int minPts);
    
    void getNeighbors(vector<unsigned int>& pts, vector<vector<unsigned int>>& indices, vector<vector<double>>& dists, vector<bool>& enoughPts, double eps);
    
    double distance(unsigned int a, unsigned int b);
    double getMin();
    double getMax();
    unsigned int getCount();
    
    void expand(Clustering& c, vector<unsigned int>& cluster);
    
private:
    vector<vector<pair<unsigned int, float>>> scan;
    bool loadIndex(string index_name, unsigned int elements);
    bool saveIndex(string index_name);
    
    vector<unsigned int> id;
    vector<vector<unsigned int>> trans_map;
    vector<vector<unsigned int>> tmp_ind;
};

#endif
