#ifndef HIST_H
#define HIST_H

//#include <algorithm>
//#include <mutex>
//#include <omp.h>
//#include <iostream>
#include <vector>
#include <mutex>

#include "metric.h"

using namespace std;

class Hist {
public:
    Hist(Met* metric, unsigned int num_bins, double min, double max);
    Hist(Met* metric, unsigned int num_bins);
    
    void hist(vector<unsigned int> &counts, vector<double> &limits, unsigned int elements);
    
    void prep_add(vector<unsigned int> &counts, vector<double> &limits);
    unsigned int add(double dist, vector<unsigned int> &counts);
    
private:
    Met* metric;
    std::mutex* access;// = new mutex[elements];

    double min_val;
    double max_val;
    unsigned int num_bins;
    double bin_width;
};

//http://stackoverflow.com/questions/4515874/searching-for-fast-efficient-histogram-algorithm

//bin_width = (max-min)/num_bins;
//bin_id = (int)((value - val_min) / bin_width);

#endif
