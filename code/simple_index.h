#ifndef SIMPLEINDEX_H
#define SIMPLEINDEX_H

#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>

#include <string>
#include <vector>
//#include <tuple>
#include <utility>
#include <mutex>

using namespace std;

template <typename Dist, typename Matrixtype, typename Distancetype>
class SimpleIndex {
public:

    SimpleIndex(){};
    ~SimpleIndex(){};
    
    unsigned int search(unsigned int i, double eps, vector<unsigned int>& indices, vector<double>& dists){
        double dist;
        for (auto it = index[i].begin(); it != index[i].end(); it++){
            dist = _factor * get<1>(*it);
            if (dist <= eps){
                indices.push_back(get<0>(*it));
                dists.push_back(dist);
            }
        }
        
        return indices.size();
    };
    
    bool loadIndex(const string index_name, unsigned int elements, double factor){
        index.resize(elements);
        _factor = 1.0 / factor;
        
        ifstream is;
        is.open(index_name);
        
        if (is.is_open()){
            string line;
            unsigned int j;
            //Distancetype dist;
            int dist;
            
            unsigned int i = 0;
            while(getline(is, line)){
                if (line != "") {
                    istringstream iss(line);
                    
                    while (iss >> j >> dist) {
                        //Distancetype dist2 = (Distancetype) ((int)dist - 33);
                        index[i].push_back(make_pair(j, (Distancetype) dist));
                        index[j].push_back(make_pair(i, (Distancetype) dist));
                    }
                }
                i++;
            }
            is.close();
            return true;
        } else {
            return false;
        }
    };
    
    bool saveIndex(const string index_name){
        ofstream is;
        is.open(index_name);
        
        if (is.is_open()){
            //write
            unsigned int i = 0;
            for (auto it_in = index.begin(); it_in != index.end(); it_in++){
                for (auto it = it_in->begin(); it != it_in->end(); it++){
                    if (get<0>(*it) > i){
                        is << get<0>(*it) << " " << ((int) get<1>(*it)) << " ";
                        //cout << (int) get<1>(*it) << endl;
                    }
                }
                is << endl;
                i++;
            }
            is.close();
            return true;
        } else {
            return false;
        }
    };
    
    void buildIndex(double eps, unsigned int length, unsigned int elements, Matrixtype &dataset, Dist _distance, double factor){
        index.resize(elements);
        _factor = 1.0 / factor;
        
        std::mutex* access = new std::mutex[elements];
        
        Distancetype eps_factor = (Distancetype) (eps * factor);
        
        //Dist _distance;
        #pragma omp parallel default(shared) num_threads(4)
        {
        #pragma omp for schedule(dynamic)
        for(unsigned int i=0; i<elements-1; i++){
            for(unsigned int j=i+1; j<elements; j++){
                Distancetype dist = (Distancetype) (factor * _distance(dataset[i], dataset[j], length));
                //cout << dist << " " << eps_factor << endl;
                if (dist < eps_factor){
                    access[i].lock();
                    index[i].push_back(make_pair(j, dist));
                    access[i].unlock();
                    
                    access[j].lock();
                    index[j].push_back(make_pair(i, dist));
                    access[j].unlock();
                }
            }
        }
        }
        
        delete[] access;
    };
    
private:
    vector<vector<pair<unsigned int, Distancetype> > > index;
    double _factor;
    //Dist distance_;
};

#endif
