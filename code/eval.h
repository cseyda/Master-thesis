#ifndef EVAL_H
#define EVAL_H

//#include <sstream>
#include <string>
#include <vector>
//#include <omp.h>
//#include <algorithm>
//#include <iostream>

#include "clustering.h"
#include "metric.h"

using namespace std;

class Eval {
public:
    Eval(Met* metric, Clustering &c);
    
    ~Eval();

    void all(double &dunn, double &db, double &c, double &sw, bool count = false);//, unsigned int sample_make, unsigned int sample_skip);
    
    void setMetric(Met* new_metric);
    
private:
    Met* metric;
    vector<vector<unsigned int>> clustering;
};

#endif
