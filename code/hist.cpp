#include <algorithm>
#include <mutex>
#include <omp.h>
#include <iostream>

#include "hist.h"
#include "metric.h"

using namespace std;

Hist::Hist(Met* metric, unsigned int num_bins, double min, double max) :
        metric(metric),
        min_val(min),
        max_val(max),
        num_bins(num_bins) {};

Hist::Hist(Met* metric, unsigned int num_bins) :
        metric(metric),
        min_val(metric->getMin()),
        max_val(metric->getMax()),
        num_bins(num_bins) {};

void Hist::hist(vector<unsigned int> &counts, vector<double> &limits, unsigned int elements){
    //std::mutex* access = new std::mutex[num_bins];
    
    counts.resize(num_bins, 0);
    limits.resize(num_bins+1, 0.0);
    
    vector<unsigned int> t_counts[4];
    for (unsigned int i=0; i<4; i++){
        t_counts[i].resize(num_bins, 0);
    }
    
    bin_width = (max_val - min_val) / num_bins;
    
    for (unsigned int i=0; i < limits.size(); i++){
        limits[i] = min_val + i * bin_width;
    }
    cout << "hist\n";
    
    #pragma omp parallel default(shared) num_threads(4)
    {
    int tid = omp_get_thread_num();
    #pragma omp for schedule(dynamic)
    for(unsigned int i=0; i<elements-1; i++){
        for(unsigned int j=i+1; j<elements; j++){
            double dist = metric->distance(i, j);
            unsigned int bin_id = (unsigned int)((dist - min_val) / bin_width);
            if ((0 <= bin_id) && (bin_id < num_bins)){
                t_counts[tid][bin_id]++;
            } else if (dist == max_val) {
                t_counts[tid][bin_id-1]++;
            } else {
                cout << "wrong bin\n";
                cout << bin_id << " " << num_bins << endl;
                cout << "dist: " << dist << " max: " << max_val << endl;
            }
        }
    }
    }
    
    for (unsigned int i=0; i<4; i++){
        for (unsigned int k=0; k<t_counts[i].size(); k++){
            counts[k] += t_counts[i][k];
        }
    };
    //delete[] access;
};

void Hist::prep_add(vector<unsigned int> &counts, vector<double> &limits){
    counts.resize(num_bins, 0);
    limits.resize(num_bins+1, 0.0);
    access = new std::mutex[num_bins];

    bin_width = (max_val - min_val) / num_bins;
    //cout << bin_width << " " << max_val << " " << min_val << " " << num_bins << endl;
    for (unsigned int i=0; i < limits.size(); i++){
        limits[i] = min_val + i * bin_width;
    }
};
unsigned int Hist::add(double dist, vector<unsigned int> &counts){
    unsigned int bin_id = (unsigned int)((dist - min_val) / bin_width);
    if ((0 <= bin_id) && (bin_id < num_bins)){
        access[bin_id].lock();
        counts[bin_id]++;
        access[bin_id].unlock();
    } else if (dist == max_val) {
        access[bin_id-1].lock();
        counts[bin_id-1]++;
        access[bin_id-1].unlock();
    } else {
        cout << "wrong bin\n";
        cout << bin_id << " " << num_bins << endl;
        cout << "dist: " << dist << " max: " << max_val << endl;
    }
    
    return bin_id;
};
