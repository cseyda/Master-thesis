//#include <sstream>
#include <string>
#include <vector>
#include <omp.h>
#include <algorithm>
#include <iostream>
#include <mutex>

#include <cmath>

#include "metric.h"
#include "eval.h"
#include "hist.h"
#include "clustering.h"

using namespace std;

double max_(Met* metric, unsigned int elements){
    double maxes[4] = {0};
    
    #pragma omp parallel default(shared) num_threads(4)
    {
    int tid = omp_get_thread_num();
    #pragma omp for schedule(dynamic)
    for(unsigned int i=0; i<elements-1; i++){
        for(unsigned int j=i+1; j<elements; j++){
            double dist = metric->distance(i, j);
            if (dist >= maxes[tid]){
                maxes[tid] = dist;
            }
        }
    }
    }
    return *max_element(maxes,maxes+3);
};

Eval::Eval(Met* metric, Clustering &c):
        metric(metric), clustering(c.clustering) {
    
};

Eval::~Eval(){};


void Eval::setMetric(Met* new_metric){
    metric = new_metric;
};

void Eval::all(double &dunn, double &db, double &c, double &sw, bool count){
    double metric_max = metric->getMax();
    double metric_min = metric->getMin();
    unsigned int elements = metric->getCount();
    
    unsigned int clusters = clustering.size();
    cout << "Clusters: " << clusters << " Nodes: " << elements << endl;
    
    vector<unsigned int> id(elements, 0);
    unsigned int clustered_nodes=0;
    for(auto it:clustering){
        for(unsigned int nid : it){
            id[nid] = clustered_nodes;
            clustered_nodes++;
        }
    }
    elements = clustered_nodes;
    
    // DB
    vector<double> db_sqr(elements, 0.0);
    
    // DB end

    // SW
    double ASW[4] = {0.0};
    vector<std::mutex> access(elements);
    
    vector<vector<double>> matrix(elements, vector<double>(clustering.size(), 0.0));
    
    
    // SW end
    
    // C INDEX
    unsigned int n=0;
    double S=0.0;

    //histo
    Hist h(metric, /*bins*/500000);

    vector<unsigned int> counts;
    vector<double> limits;
    h.prep_add(counts, limits);
    // C INDEX end

    // DUNN
    double max[4];
    double min[4];
    for(unsigned int i=0; i<4;i++){
        max[i]=metric_min;
        min[i]=metric_max;
    }
    // DUNN end
    
    if (!count){
    //inter-cluster distance
    #pragma omp parallel default(shared) num_threads(4)
    {
    int tid = omp_get_thread_num();
    
    for (unsigned int i=0; i < clusters-1; i++){
        for (unsigned int j=i+1; j < clusters; j++){
            #pragma omp for schedule(dynamic)
            for (auto a = clustering[i].begin(); a < clustering[i].end(); a++){
                for (auto b = clustering[j].begin(); b < clustering[j].end(); b++){
                    double dist = metric->distance(*a, *b);
                    // *a from cluster i
                    // *b from cluster j
                    unsigned int a_ = id[*a];
                    unsigned int b_ = id[*b];
                    // DUNN
                    if (dist < min[tid]) {
                        min[tid] = dist;
                    }
                    // DUNN end
                    
                    // C INDEX
                    h.add(dist, counts);
                    // C INDEX end
                    
                    // DB
                    
                    //DB end
                    
                    // SW
                    access[a_].lock();
                    matrix[a_][j]+=dist;// / clustering[j].size();
                    access[a_].unlock();
                    
                    access[b_].lock();
                    matrix[b_][i]+=dist;// / clustering[i].size();
                    access[b_].unlock();
                    // SW end
                }
            }
        }
    }
    
    //intra-cluster distance
    for (unsigned int i=0; i < clustering.size(); i++){
        #pragma omp for schedule(dynamic) reduction(+:S,n)
        for (auto a=clustering[i].begin(); a < clustering[i].end()-1; a++){
            for (auto b=a+1; b < clustering[i].end(); b++){

                double dist = metric->distance(*a, *b);
                // *a and *b from cluster i
                unsigned int a_ = id[*a];
                unsigned int b_ = id[*b];
                //cout << dist << endl;
                // C INDEX
                S+=dist;
                n++;
                // C INDEX end
                
                // DUNN
                if (dist != 0.0){
                    if (dist > max[tid]) {
                        max[tid] = dist;
                    }
                } else {
                    if (0.0000001 > max[tid]) {
                        max[tid] = 0.0000001;
                    }
                }
                // DUNN end
                
                // C INDEX
                h.add(dist, counts);
                // C INDEX end
                
                // SW
                // DB
                access[a_].lock();
                matrix[a_][i] += dist;
                db_sqr[a_] += dist*dist;
                access[a_].unlock();
                
                access[b_].lock();
                matrix[b_][i] += dist;
                db_sqr[b_] += dist*dist;
                access[b_].unlock();
                //DB end
                // SW end
            }
        }
    }
    }//end omp
    } // end if count
    else {
                unsigned long int cal_counter = 0;

        for (unsigned int i=0; i < clustering.size(); i++){
            cal_counter+=(clustering[i].size() * (clustering[i].size()-1)) / 2;
        }
        cout << "IntraCounter: " << cal_counter << endl;
        
        for (unsigned int i=0; i < clusters-1; i++){
            for (unsigned int j=i+1; j < clusters; j++){
                cal_counter+=clustering[i].size() * clustering[j].size();
            }
        }
        
        cout << "Counter: " << cal_counter << endl;
    }
    
    //
    //C INDEX
    //
    double S_min = 0.0;
    unsigned int n_min = n;
    for (unsigned int i=0; i<counts.size(); i++){
        double mean = (limits[i] + limits[i+1]) / 2.0;
        if (counts[i] < n_min) {
            S_min += counts[i]*mean;//limits[i];//
            n_min -= counts[i];
        } else {
            S_min += n_min*mean;//limits[i];//
            break;
        }
    }
    
    double S_max = 0.0;
    unsigned int n_max = n;
    for (unsigned int i=counts.size()-1; i>=0; i--){
        double mean = (limits[i] + limits[i-1]) / 2.0;
        if (counts[i] < n_max) {
            S_max += counts[i]*mean;//limits[i];//
            n_max -= counts[i];
        } else {
            S_max += n_max*mean;//limits[i];//
            break;
        }
    }
    cout << S << " " << S_min << " " << S_max << " " << n << endl;
    c = (S - S_min) / (S_max - S_min);
    // C INDEX end
    
    //
    // DUNN
    //
    double min_v, max_v;
    min_v = *min_element(min, min+3);
    max_v = *max_element(max, max+3);
    dunn = min_v / max_v;
    cout << "Dunn: " << min_v << "/" << max_v << " -> " << dunn << endl;
    // DUNN end

    //
    // SW
    //
    Hist sw_h(metric, /*bins*/5, -1, 1);

    vector<unsigned int> sw_counts;
    vector<double> sw_limits;
    sw_h.prep_add(sw_counts, sw_limits);
    
    //#pragma omp parallel default(shared) num_threads(4)
    //{
    
    int tid = omp_get_thread_num();
    
    for (unsigned int i=0; i<clusters; i++){
    //    #pragma omp for schedule(dynamic)
        for (auto a=clustering[i].begin(); a<clustering[i].end(); a++){
            
            unsigned int a_ = id[*a];
            
            double min = metric_max;
            for (unsigned int j=0; j<clusters; j++){
                if (i!=j){
                    if (matrix[a_][j]/clustering[j].size() < min) {
                        min = matrix[a_][j]/clustering[j].size();
                    }
                }
            }

            double a_i = matrix[a_][i] / (clustering[i].size()-1);
            double b_i = min;
            
            double sw_i;
            if (a_i != b_i) {
                sw_i = (b_i - a_i) / std::max(b_i, a_i);
            } else {
                sw_i = 0.0;
            }
            ASW[tid] += sw_i;
            sw_h.add(sw_i, sw_counts);
        }
    }
    
    for (unsigned int i : sw_counts){
        cout << i << " ";
    }
    cout << endl;
    //}// end omp

    // DB
    // get rhos
    
    vector<unsigned int> cluster_medoid(clusters);
    
    for (unsigned int i=0; i < clusters; i++){
        double min = metric_max * clustering[i].size();
        unsigned int min_id = elements+1;
        
        //for (auto a=clustering[i].begin(); a<clustering[i].end(); a++){
        for(unsigned int a : clustering[i]){
            unsigned int a_ = id[a];
            
            if (db_sqr[a_] < min) {
                min = db_sqr[a_];
                min_id = a;
            }
        }
        
        cluster_medoid[i] = min_id;
    }

    db = 0.0;
    for (unsigned int i=0; i < clusters; i++){
        double max=0;
        for (unsigned int j=0; j < clusters; j++){
            if (i!=j){
                unsigned int a = cluster_medoid[i];
                unsigned int b = cluster_medoid[j];
                
                unsigned int a_ = id[a];
                unsigned int b_ = id[b];
                    
                double p_i = matrix[a_][i] / (clustering[i].size()-1);
                double p_j = matrix[b_][j] / (clustering[j].size()-1);
                
                double dist = metric->distance(a, b);
                double db_i = (p_i + p_j);
                
                if (dist != 0){
                    db_i /= dist;
                }// else {
                 //   db_i *= 1000000;
                //}
                
                if (db_i > max){
                    max = db_i;
                }
            }
        }
        db += max;
        //cout << max << " ";
    }
    //cout << db << endl;
    db /= clusters;
    
    // DB end
    
    //cout << ASW[0] + ASW[1] + ASW[2] + ASW[3] << " " << clustered_nodes << endl;
    sw = (ASW[0] + ASW[1] + ASW[2] + ASW[3]) / clustered_nodes;
    // SW end
};
