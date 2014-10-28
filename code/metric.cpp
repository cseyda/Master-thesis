#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <set>
#include <map>
#include <string>
//#include <array>
//#include <tuple>
//#include <mutex>
#include <utility>
#include <algorithm>
#include <iterator> //distance
#include <numeric>      // std::accumulate

#include "dbclient.h"
#include "flann/flann.hpp"
#include <omp.h>
#include <stdint.h>

#include "my_dists.h"
//#include "simple_index.h"
#include "metric.h"
#include "graph.h"
#include "clustering.h"

#include <algorithm>
#include <numeric>
#include <limits>

#include <chrono>


typedef std::chrono::high_resolution_clock Clock;
typedef std::chrono::milliseconds milliseconds;

using namespace std;

static const double DEG_TO_RAD = 0.017453292519943295769236907684886;
static const double EARTH_RADIUS_IN_METERS = 6372797.560856;

inline double meter2deg(double eps){
    return eps/(EARTH_RADIUS_IN_METERS * DEG_TO_RAD);
};

inline double deg2meter(double eps){
    return eps * (EARTH_RADIUS_IN_METERS * DEG_TO_RAD);
};

//
// MET
//
Met::Met(string name, string collection, double eps, unsigned int minPts) :
    _weight(1),
    _eps(eps),
    _minPts(minPts),
    _name(name), 
    _collection(collection) {};
Met::~Met(){};

string Met::getName() {return _name;};
string Met::getCollection() {return _collection;};

double Met::getEps() {return _eps;};
void Met::setEps(double eps){_eps = eps;};

unsigned int Met::getPoints(){return _minPts;};

unsigned int Met::getWeight(){return _weight;};
void Met::setWeight(unsigned int weight){_weight = weight;};

bool Met::getNeighbors(unsigned int pt, vector<unsigned int>& rgpNodesFound){
    return getNeighbors(pt, rgpNodesFound, 0, _eps, _minPts);
};

bool Met::getNeighbors(unsigned int pt, vector<unsigned int>& rgpNodesFound, vector<double>& distances){
        return getNeighbors(pt, rgpNodesFound, &distances, _eps, _minPts);
};

void Met::expand(Clustering& c, vector<unsigned int>& cluster){
    for(unsigned int nodeID=0; nodeID < cluster.size(); nodeID++){
        c.add(cluster[nodeID], nodeID);
    }
};

//
// FLANN
//

//typedef double DistType;
//typedef flann::L2_Simple<DistType> Dist;
LocF::LocF(DB2& db, double eps, unsigned int minPts) :
        Met("LocF", db.getCollection(), meter2deg(eps), minPts), 
        tree(flann::KDTreeSingleIndexParams()),
        db(&db),
        nn(2),
        dataset(new DistType[db.size()*nn], db.size(), nn),
        params(32)
        {
    params.cores = 4;
    params.sorted = true;
    
    bool first = true;
    
    unsigned int db_count = db.size();
    double longitude, latitude;
    vector<string> words;
    
    for (unsigned int i=0; i<db_count; i++){
        db.get(i, longitude, latitude, words);
        
        dataset[i][0] = longitude;
        dataset[i][1] = latitude;
        
        if (first){
            min_lon = longitude;
            max_lon = longitude;
    
            min_lat = latitude;
            max_lat = latitude;
            
            first = false;
        } else {
            if (max_lon < longitude) {
                max_lon = longitude;
            } else if (min_lon > longitude){
                min_lon = longitude;
            }
            
            if (max_lat < latitude) {
                max_lat = latitude;
            } else if (min_lat > latitude){
                min_lat = latitude;
            }
        }
    }
    
    tree.buildIndex(dataset);
    
    
    
    query_space = new DistType[200*nn];
    //delete[] dataset.ptr();
};

LocF::LocF(DB2& db, double eps, unsigned int minPts, vector<unsigned int>& pts) :
        Met("LocF", db.getCollection(), meter2deg(eps), minPts), 
        tree(flann::KDTreeSingleIndexParams()),
        db(&db),
        nn(2),
        dataset(new DistType[pts.size()*nn], pts.size(), nn),
        params(32)
        {
    params.cores = 4;
    params.sorted = true;
    
    bool first = true;
    
    //unsigned int db_count = db.size();
    double longitude, latitude;
    vector<string> words;
    
    for (unsigned int i=0; i<pts.size(); i++){
        db.get(pts[i], longitude, latitude, words);
        
        dataset[i][0] = longitude;
        dataset[i][1] = latitude;
        
        if (first){
            min_lon = longitude;
            max_lon = longitude;
    
            min_lat = latitude;
            max_lat = latitude;
            
            first = false;
        } else {
            if (max_lon < longitude) {
                max_lon = longitude;
            } else if (min_lon > longitude){
                min_lon = longitude;
            }
            
            if (max_lat < latitude) {
                max_lat = latitude;
            } else if (min_lat > latitude){
                min_lat = latitude;
            }
        }
    }
    
    tree.buildIndex(dataset);
    
    
    
    query_space = new DistType[200*nn];
    //delete[] dataset.ptr();
};

LocF::~LocF(){};

void LocF::getNeighbors(vector<unsigned int>& pts, vector<vector<unsigned int>>& indices, vector<vector<double>>& dists, vector<bool>& enoughPts, double eps){
    //flann::SearchParams params;
    //params.cores = 4;
    //params.checks = 32;
    flann::Matrix<DistType> queries(query_space, pts.size(), nn);
    //flann::Matrix<DistType> queries(new DistType[pts.size()*nn], pts.size(), nn);

    // set pts
    for(unsigned int i=0; i<pts.size(); i++){
        queries[i][0] = dataset[pts[i]][0];
        queries[i][1] = dataset[pts[i]][1];
    }
    vector<vector<int>> tmp_ind;
    
    //unsigned int number_near = 
    tree.radiusSearch(
        /*const Matrix<ElementType>&*/ queries,
        /*std::vector< std::vector<int> >&*/ tmp_ind,
        /*std::vector<std::vector<DistanceType> >&*/ dists,
        /*float radius*/ eps,
        /*const SearchParams&*/ params);
    
    
    for(unsigned int i=0; i<indices.size(); i++){
        indices[i].assign(tmp_ind[i].begin(), tmp_ind[i].end());
        enoughPts[i] = (indices[i].size() >= getPoints());
    }
}

bool LocF::getNeighbors(unsigned int pt, vector<unsigned int>& rgpNodesFound, vector<double>* distances, double eps, unsigned int minPts){
    std::vector< std::vector<int> > indices;
    std::vector<std::vector<double> > dists;
    //indices.clear();
    //dists.clear();
    
    flann::Matrix<DistType> queries(dataset[pt], 1, nn);
    
    unsigned int number_near = tree.radiusSearch(
        /*const Matrix<ElementType>&*/ queries,
        /*std::vector< std::vector<int> >&*/ indices,
        /*std::vector<std::vector<DistanceType> >&*/ dists,
        /*float radius*/ eps,
        /*const SearchParams&*/ params/*flann::SearchParams(-1)*/);
    
    // remove search point
    /*
    unsigned int found_at = 0;//? uninitialized
    for (unsigned int i=0; i<indices[0].size(); i++){
        if(pt == (unsigned int) indices[0][i]){
            number_near--;
            found_at = i;
        } else {
            rgpNodesFound.push_back(indices[0][i]);
        }
    }*/
    
    rgpNodesFound.clear();
    rgpNodesFound.assign(indices[0].begin(), indices[0].end());
    
    /*if (!std::is_sorted(dists[0].begin(),dists[0].end())){
        cout << "Not sorted" << endl;
        for(auto i:dists[0]){
            cout << i << " ";
        }
        cout << endl;
    }*/
        
    if (distances != 0){
        //dists[0][found_at] = dists[0][dists[0].size()-1];
        //dists[0].erase(dists[0].end()-1);

        *(distances) = move(dists[0]);
        
        //transform(distances->begin(),distances->end(),distances->begin(),
        //    bind2nd(std::multiplies<double>(),(EARTH_RADIUS_IN_METERS * DEG_TO_RAD)));
    }

    return number_near >= minPts;
};

double LocF::distance(unsigned int a, unsigned int b){
    //return _distance(dataset[a], dataset[b], nn) * (EARTH_RADIUS_IN_METERS * DEG_TO_RAD);
    return _distance(dataset[a], dataset[b], nn);
};

double LocF::getMin(){
    return 0.0;
};
double LocF::getMax(){
    /*
    double maxes[4] = {0};
    unsigned int elements = vecNodes->size();
    
    #pragma omp parallel default(shared) num_threads(4)
    {
    int tid = omp_get_thread_num();
    #pragma omp for schedule(dynamic)
    for(unsigned int i=0; i<elements-1; i++){
        for(unsigned int j=i+1; j<elements; j++){
            double dist = distance(i, j);
            if (dist > maxes[tid]){
                maxes[tid] = dist;
            }
        }
    }
    }
    maxes[0] = *max_element(maxes,maxes+3);
    return maxes[0];
    */

    double min[2] = {min_lat, min_lon};
    double max[2] = {max_lat, max_lon};
    //cout << min[0] << " " << min[1] << " " << max[0] << " " << max[1] << endl;

    double dist = _distance(min, max, nn);
    //cout << dist << endl;

    return dist;
};

unsigned int LocF::getCount(){
    return db->size();
};



void bigram_map(vector<string> &bigram_vec, map<string, unsigned int> &bigrams, unsigned int &nn){
    if (bigram_vec.empty()){
        unsigned int chars = 42;
        string character[42] = {"a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z", " ", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "_", ".", ",", "!", "?"}; 
        
        nn = chars*chars;
        
        bigram_vec.resize(nn);
        for (unsigned int i=0; i<chars; i++){
            for (unsigned int j=0; j<chars; j++){
                bigram_vec[i*chars+j] = character[i] + character[j];
            }
        }
    } else {
        nn = bigram_vec.size();
    }

    // set bigram map for vector creation
    unsigned int j=0;
    for (auto it = bigram_vec.begin();
            it != bigram_vec.end(); 
            it++, j++){
        bigrams[*it] = j;
    }
};


/*
template <class DistType>
void bigram_vector(vector<string> &bigram_vec, flann::Matrix<DistType> &dataset, vector<Node*> &db, unsigned int &nn){
    //create bigram string, if not present
    
    map<string, unsigned int> bigrams;
    
    bigram_map(bigram_vec, bigrams, nn);
    
    dataset = flann::Matrix<DistType>(new DistType[db.size()*nn], db.size(), nn);
    for (unsigned int i=0; i<db.size(); i++){
        for (unsigned int j=0; j<nn; j++){
            dataset[i][j] = 0;
        }
    }
    
    unsigned int c = 0;
    for (vector<Node*>::const_iterator it = vecNodes.begin(); 
            it != vecNodes.end(); 
            it++, c++){
            
        Node* pt = *it;
        // extract character bigrams from string
        for(unsigned int l = 0; l < (pt->tag.length() - 1); l++) {
            unsigned int index = bigrams[pt->tag.substr(l, 2)];
            dataset[c][index]+=1;
        }
    }
};*/





template <class DistType>
unsigned int bigram_list(vector<string> &bigram_vec, flann::Matrix<DistType> &dataset, vector<Node*> &vecNodes, unsigned int &nn){
    //create bigram string, if not present
    
    map<string, unsigned int> bigrams;
    
    bigram_map(bigram_vec, bigrams, nn);

    //heuristic: max bigram-len => ignore doubles
    unsigned int list_len = 0;
    for (vector<Node*>::const_iterator it = vecNodes.begin(); 
            it != vecNodes.end(); 
            it++){
        Node* pt = *it;
        
        if (pt->tag.length() > list_len){
            list_len = pt->tag.length();
        }
    }
    
    dataset = flann::Matrix<DistType>(new DistType[vecNodes.size()*list_len], vecNodes.size(), list_len);
    for (unsigned int i=0; i<vecNodes.size(); i++){
        for (unsigned int j=0; j<list_len; j++){
            dataset[i][j] = 0;
        }
    }
    
    set<unsigned int> tmp;
    unsigned int c = 0;
    for (vector<Node*>::const_iterator it = vecNodes.begin(); 
            it != vecNodes.end(); 
            it++, c++){
        tmp.clear();
        Node* pt = *it;
        
        // extract character bigrams from string
        for(unsigned int l = 0; l < (pt->tag.length() - 1); l++) {
            unsigned int index = bigrams[pt->tag.substr(l, 2)];
            //dataset[c][index]+=1;
            tmp.insert(index+1);
        }
        
        unsigned int pos = 0;
        for (auto it_set=tmp.begin(); it_set != tmp.end(); it_set++){
            dataset[c][pos] = *it_set;
            pos++;
        }
    }
    
    return list_len;
};



template <class DistType>
unsigned int jac_list(vector<vector<unsigned int>> &jac_vec, flann::Matrix<DistType> &dataset, DB2& db, unsigned int &nn){
    
    unsigned int list_len = 0;
    for (auto vec : jac_vec){
        if (vec.size() > list_len){
            list_len = vec.size();
        }
    }
    
    dataset = flann::Matrix<DistType>(new DistType[db.size()*list_len], db.size(), list_len);
    for (unsigned int i=0; i<db.size(); i++){
        for (unsigned int j=0; j<list_len; j++){
            dataset[i][j] = 0;
        }
    }
    
    for (unsigned int i=0; i<db.size(); i++){
        for(unsigned int j=0; j<jac_vec[i].size(); j++){
            dataset[i][j] = jac_vec[i][j];
        }
    }
    
    return list_len;
};
//
// Simple Index
//

//typedef unsigned char DistType;
//typedef Jaccard<DistType> Dist;

JaccardIndex::JaccardIndex(DB2& db, double eps, unsigned int minPts, string kind, bool use_index) :
        Met("JaccardIndex", db.getCollection(), eps, minPts),
        db(&db),
        kind(kind) {
    
    //bigram
    //words
    if ("bigram" != kind){
        kind = "words";
    }
    
    //unsigned int tmp = bigram_list<DistType>(bigram_vec, dataset, vecNodes, nn);
    vector<vector<unsigned int>> jac_vec;
    
    db.jaccard(jac_vec, kind);
    
    unsigned int tmp = jac_list<DistType>(jac_vec, dataset, db, nn);
    bigram_ids=nn;
    nn=tmp;
    
    string index_name = getName() + "_" + getCollection() + "_" + kind;
    
    if (use_index){
        if(!index.loadIndex(index_name, db.size(), 100.0)){
            cout << "no tree found\n";
            index.buildIndex(
                    //double eps 
                    0.7, 
                    //unsigned int vector_len 
                    nn, 
                    //unsigned int elements 
                    db.size(), 
                    //&Matrixtype 
                    dataset, 
                    //Dist ancetype 
                    _distance, 
                    //double factor (scaling) 
                    100.0);
            index.saveIndex(index_name);
        } else {cout << "tree loaded\n";}
    } else {cout << "no index used\n";}
};

JaccardIndex::~JaccardIndex(){};

bool JaccardIndex::getNeighbors(unsigned int pt, vector<unsigned int>& rgpNodesFound, vector<double>* distances, double eps, unsigned int minPts){
    //vector<unsigned int> indices;
    vector<double> dists;
    
    unsigned int number_near = index.search(pt, eps, 
        //vector<unsigned int> &
        rgpNodesFound, 
        //vector<double> 
        dists);
    if (distances != 0){
        distances->insert(distances->begin(), dists.begin(), dists.end());
    }
    return (number_near >= minPts);
};

void JaccardIndex::getNeighbors(vector<unsigned int>& pts, vector<vector<unsigned int>>& indices, vector<vector<double>>& dists, vector<bool>& enoughPts, double eps){
    //vector<unsigned int> tmp_ind;
    for(unsigned int i=0; i<pts.size(); i++){
        //tmp_ind.clear();
        enoughPts[i] = getNeighbors(pts[i], indices[i], 0, eps, getPoints());
        //indices[i].assign(tmp_ind.begin(), tmp_ind.end());
    }
};

double JaccardIndex::distance(unsigned int a, unsigned int b){
    return _distance(dataset[a], dataset[b], nn);
};

double JaccardIndex::distance(unsigned int a, unsigned int *b){
    return _distance(dataset[a], b, nn);
};
double JaccardIndex::distance(unsigned int *a, unsigned int *b){
    return _distance(a, b, nn);
};
    
    
double JaccardIndex::getMin(){
    return 0.0;
};
double JaccardIndex::getMax(){
    return 1.0;
};

unsigned int JaccardIndex::getCount(){
    return db->size();
};


//
// COMBINED
//
CombinedMetric::CombinedMetric(DB2& db, vector<Met*> *metrics, double eps, unsigned int minPts) :
        Met("Combined", db.getCollection(), eps, minPts),
        _metrics(metrics),
        _factor(metrics->size()),
        db(&db){

    unsigned int weight_sum = 0;
    for (auto it = _metrics->begin(); it != _metrics->end(); it++){
        weight_sum += (*it)->getWeight();
    };
    
    unsigned int i=0;
    for (auto it = _metrics->begin(); 
            it != _metrics->end(); 
            it++, i++){
        _factor[i] = (double) (*it)->getWeight() / (weight_sum * (*it)->getEps());
        cout << _factor[i] << " " << (*it)->getEps() << endl;
    }
    
    //compressed
    try {
        vector<unsigned int> edges;
        unsigned int id_count = 0;
        
        cout << "graph" << endl;
        //get graph weighted nodes
        db.graph(id, edges, id_count, 0);
        
        trans_map.resize(id_count);
        
        for (unsigned int i=0; i<id.size(); i++){
            trans_map[id[i]].push_back(i);
        }
        
        cout << "alter" << endl;
        //alter based on text distances
        cout << "trans_map: " << trans_map.size() << endl;
        for(unsigned int i=0; i<trans_map.size(); i++){
            //used for comparison
            unsigned int first = trans_map[i][0];
            
            // every element before bound is near enough to stay
            auto bound = std::partition(trans_map[i].begin()+1, 
                    trans_map[i].end(),
                    [&](unsigned int id){return _metrics->at(1)->distance(id, first) < 0.1;});
            if (bound != trans_map[i].end()){
                trans_map.push_back(vector<unsigned int>(bound, trans_map[i].end()));
                trans_map[i].erase(bound, trans_map[i].end());
            }
        }
        cout << "trans_map: " << trans_map.size() << endl;
        for(unsigned int i=0; i<trans_map.size(); i++){
            for(unsigned int j=0; j<trans_map[i].size(); j++){
                id[trans_map[i][j]] = i;
            }
        }
        
        vector<unsigned int> tmp_pts;
        for(auto pts : trans_map){
            tmp_pts.push_back(pts[0]);
        }
        
        cout << "new LocF" << endl;
        //new LocF, using the new trans ids
        delete _metrics->at(0);
        _metrics->at(0) = new LocF(db, 10000.0, 2, tmp_pts);
        
    } catch(...) {
        id.resize(db.size());
        trans_map.resize(db.size());
        
        for (unsigned int i=0; i<id.size(); i++){
            id[i]=i;
            trans_map[i].push_back(i);
        }
    }
    
    //
    //
    
    
    // factor: multiply to normalize
    // weight: indicates importance
};

CombinedMetric::~CombinedMetric(){};

bool CombinedMetric::getNeighbors(unsigned int pt, vector<unsigned int>& rgpNodesFound, vector<double>* distances, double eps, unsigned int minPts){
    vector<double> my_dists;
    
    //get initial points
    _metrics->at(0)->getNeighbors(pt, rgpNodesFound, &my_dists, eps / _factor[0], 2);
    //find first not 0 distance
    unsigned int pos_not_zero = std::distance(my_dists.begin(),
        std::find_if(my_dists.begin(),
                    my_dists.end(), 
                    [](double val){return val >= 0.00009;}));
    //0.00009 ~ 10m exact
    // BS!! ??
    
    //normalize distances
    //transform(my_dists.begin()+pos_not_zero,my_dists.end(),my_dists.begin(),
        //bind2nd(std::multiplies<double>(),_factor[0]));
    
    
    // generic version
    /*
    #pragma omp parallel default(shared) num_threads(4)
    {
    #pragma omp for schedule(dynamic)
    for (unsigned int j=pos_not_zero; j<rgpNodesFound.size(); j++){
        my_dists[j] *= _factor[0];
        
        for (unsigned int i=1; i<_metrics->size(); i++){
            my_dists[j] += _factor[i] * _metrics->at(i)->distance(pt, rgpNodesFound[j]);
        }
    }
    }
    */
    
    //#pragma omp parallel default(shared) num_threads(4)
    //{
    //#pragma omp for schedule(dynamic)
    
    //cout << "g " << pt << " g" << endl;
    if (pt == 3074){
    cout << _factor[0] << " " << my_dists[1] << "  " << _factor[1] << " " << _metrics->at(1)->distance(trans_map[pt][0], trans_map[rgpNodesFound[1]][0]) << endl;
    
    cout << _factor[0] * my_dists[1] + _factor[1] * _metrics->at(1)->distance(trans_map[pt][0], trans_map[rgpNodesFound[1]][0]) << endl;
    }
    
    for (unsigned int j=pos_not_zero; j<rgpNodesFound.size(); j++){
        my_dists[j] = _factor[0] * my_dists[j]
                    + _factor[1] * _metrics->at(1)->distance(trans_map[pt][0], trans_map[rgpNodesFound[j]][0]);
    }
    //}
    
    
    
    //cout << "c" << endl;
    
    //remove elements greater distance than eps
    //removing using stl
    /*rgpNodesFound.erase(
        remove_if(rgpNodesFound.begin()+pos_not_zero, rgpNodesFound.end(),
        [&](double val) { return eps < val; }),
        rgpNodesFound.end());
    */
    unsigned int b_eps = 0;
    for (unsigned int i=pos_not_zero; i<my_dists.size()-b_eps; i++){
        if (my_dists[i] > eps){
            my_dists[i] = my_dists[my_dists.size()-b_eps-1];
            rgpNodesFound[i] = rgpNodesFound[my_dists.size()-b_eps-1];
            
            b_eps++;
            i--;
        }
    }
    rgpNodesFound.erase(rgpNodesFound.end()-b_eps, rgpNodesFound.end());
    
    
    
    //update found
    unsigned int found = std::accumulate(rgpNodesFound.begin(), rgpNodesFound.end(), 0, [&](unsigned int x, unsigned int id){return x + trans_map[id].size();});
    
    return (found >= minPts);
};

void CombinedMetric::getNeighbors(vector<unsigned int>& pts, vector<vector<unsigned int>>& indices, vector<vector<double>>& dists, vector<bool>& enoughPts, double eps){
    #pragma omp parallel default(shared) num_threads(1)
    {
    #pragma omp for schedule(dynamic)
    
    for(unsigned int i=0; i<pts.size(); i++){
        enoughPts[i] = getNeighbors(pts[i], indices[i], 0, eps, getPoints());
    }
    
    }
};

//not compressed ids
double CombinedMetric::distance(unsigned int a, unsigned int b) {
    double dist = 0;
    
    unsigned int i = 0;
    for(auto it = _metrics->begin(); 
            it != _metrics->end(); 
            it++, i++){
        dist += _factor[i] * (*it)->distance(a, b);
    }
    return dist;
};

double CombinedMetric::getMin(){
    return 0.0;
};
double CombinedMetric::getMax(){
    double dist = 0;
    
    unsigned int i = 0;
    for(auto it = _metrics->begin(); 
            it != _metrics->end(); 
            it++, i++){
        dist += _factor[i] * (*it)->getMax();
    }
    return dist;
};

vector<Met*>* CombinedMetric::metrics(){
    return _metrics;
};

unsigned int CombinedMetric::getCount(){
    //return _metrics->at(0)->getCount();
    return trans_map.size();
};

void CombinedMetric::expand(Clustering& c, vector<unsigned int>& cluster){
    for(unsigned int nodeID=0; nodeID < cluster.size(); nodeID++){
        for(unsigned int mapped : trans_map[nodeID]){
            c.add(cluster[nodeID], mapped);
        }
    }
};


//
// Graph
//
GraphMetric::GraphMetric(DB2& db, double eps, unsigned int minPts, double w_1, double c, unsigned int dist_limit, unsigned int w_threshold, bool compressed):
        Met("Graph", db.getCollection(), eps, minPts),
        //scan(db.size()),
        tmp_ind(4)
{
    string index_name = getName() + "_" + getCollection() + "_" + to_string(c) + "_" + to_string(w_1) + "_" + to_string(dist_limit) + "_" + to_string(w_threshold);
    
    if (compressed){
        index_name += "_compressed";
        
        vector<unsigned int> edges;
        unsigned int id_count = 0;
        
        db.graph(id, edges, id_count, 0);
        
        trans_map.resize(id_count);
        scan.resize(id_count);
    } else {
        id.resize(db.size());
        
        trans_map.resize(db.size());
        scan.resize(db.size());
        
        for (unsigned int i=0; i<id.size(); i++){
            id[i]=i;
        }
    }
    
    cout << id.size() << " " << scan.size() << endl;
    
    for (unsigned int i=0; i<id.size(); i++){
        trans_map[id[i]].push_back(i);
    }
    
    if(!loadIndex(index_name, scan.size())){
        cout << "no tree found\n";
        Graph2 graph(db, &scan, w_1, c, dist_limit, w_threshold, compressed);
        saveIndex(index_name);
    } else {cout << "tree loaded\n";}
};

GraphMetric::~GraphMetric(){};

struct sort_pred {
    bool operator()(const std::pair<unsigned int,float> &left, const std::pair<unsigned int,float> &right) {
        //cout << left.second << " " << right.second << endl;
        return left.second < right.second;
    }
};
        
bool GraphMetric::loadIndex(string index_name, unsigned int elements){
    //scan.resize(elements);
    
    ifstream is;
    is.open(index_name);
    
    if (is.is_open()){
        //unsigned int z_counter = 0;
        //unsigned int else_counter = 0;
        
        string line;
        unsigned int j;
        //Distancetype dist;
        float dist;
        
        unsigned int i = 0;
        while(getline(is, line)){
            if (line != "") {
                istringstream iss(line);
                
                while (iss >> j >> dist) {
                    scan[i].push_back(make_pair(j, dist));
                    scan[j].push_back(make_pair(i, dist));
                }
            }
            i++;
        }
        is.close();
        
        //cout << "zero: " << z_counter << "\nelse: " << else_counter << endl;
        //http://stackoverflow.com/questions/279854/how-do-i-sort-a-vector-of-pairs-based-on-the-second-element-of-the-pair
        for (auto &v : scan){
            std::sort(v.begin(), v.end(), sort_pred());
        }
        
        return true;
    } else {
        return false;
    }
};

bool GraphMetric::saveIndex(string index_name){
    ofstream is;
    is.open(index_name);
    
    if (is.is_open()){
        //write
        unsigned int i = 0;
        for (auto &it_in : scan){
            /*
            for (auto &it : it_in){
                if (it.first > i){
                    is << it.first << " " << it.second << " ";
                }
            }*/
            //std::copy_if (it_in.begin(), it_in.end(), std::ostream_iterator(is, ), [i](pair<unsigned int, float> &p){return p.first > i;});
            
            std::transform(it_in.begin(), it_in.end(), 
                    std::ostream_iterator<string>(is, " "), 
                    [&](pair<unsigned int, float> &p) -> string {return (p.first > i) ? to_string(p.first)+" "+to_string(p.second) : "";}
                );
            
            is << "\n";
            i++;
        }
        is.close();
        return true;
    } else {
        return false;
    }
};

bool GraphMetric::getNeighbors(unsigned int pt, vector<unsigned int>& rgpNodesFound, vector<double>* distances, double eps, unsigned int minPts){
    unsigned int found = trans_map[pt].size();
    /*
    for (auto it = scan[pt].begin(); it != scan[pt].end(); it++){
        double dist = it->second;

        if (dist <= eps){
            found += trans_map[it->first].size();
            rgpNodesFound.push_back(it->first);
        } else {
            break;
        }
    }
    */
    
    //get range
    auto end_iter = std::find_if(scan[pt].begin(), scan[pt].end(), [eps](pair<unsigned int, float> &p){return p.second >= eps;});
    
    //copy elements
    std::transform(scan[pt].begin(), end_iter, back_inserter(rgpNodesFound), [](pair<unsigned int, float> &p){return p.first;});
    
    //update found
    found = std::accumulate(rgpNodesFound.begin(), rgpNodesFound.end(), 0, [&](unsigned int x, unsigned int id){return x + trans_map[id].size();});
    
    return (found >= minPts);
};

void GraphMetric::getNeighbors(vector<unsigned int>& pts, vector<vector<unsigned int>>& indices, vector<vector<double>>& dists, vector<bool>& enoughPts, double eps){
    #pragma omp parallel default(shared) num_threads(1)
    {
    //unsigned int tid = omp_get_thread_num();
    
    #pragma omp for schedule(dynamic)
    for(unsigned int i=0; i<pts.size(); i++){
        enoughPts[i] = getNeighbors(pts[i], indices[i], 0, eps, getPoints());
    }
    }
};

void GraphMetric::expand(Clustering& c, vector<unsigned int>& cluster){
    for(unsigned int nodeID=0; nodeID < cluster.size(); nodeID++){
        for(unsigned int mapped : trans_map[nodeID]){
            c.add(cluster[nodeID], mapped);
        }
    }
};

double GraphMetric::distance(unsigned int a, unsigned int b){
    /*
    for (auto it = scan[a].begin(); it != scan[a].end(); it++){
        if (it->first == b){
            return it->second;
        }
    }*/
    auto iter = std::find_if(scan[a].begin(), scan[a].end(), [b](pair<unsigned int, float> &p){return p.first == b;});
    
    /*if (iter != scan[a].end()){
        return iter->second;
    } else {
        return getMax();
    }*/
    return (iter != scan[a].end()) ? iter->second : getMax();
};

inline double GraphMetric::getMin(){
    return 0.0;
};

inline double GraphMetric::getMax(){
    return 1.0;
};

unsigned int GraphMetric::getCount(){
    return scan.size();
};
