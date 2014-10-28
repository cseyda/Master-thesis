#ifndef MONGO_H
#define MONGO_H

#include <string>
#include <vector>

#include "dbclient.h"
#include "node.h"

using namespace std;
using namespace mongo;

class DB {
public:
    DB();
    
    ~DB();
    
    void nodes(vector<Node*> &vecNodes, string collection);
    
    void bigrams(vector<string> &bigram_vec, string collection);

    void jaccard(vector<vector<unsigned int>> &jac_vec, string collection);

private:
    DBClientConnection conn;
    string addr;
    string port;
};


class DB2 {
public:
    DB2(string collection);
    ~DB2();
    
    unsigned int size();
    
    void get(unsigned int id, double& longitude, double& latitude, vector<string>& words);
    void get(unsigned int id, double& longitude, double& latitude, vector<unsigned int>& words);
    
    string getCollection();
    
    void jaccard(vector<vector<unsigned int>> &jac_vec, string kind);
    
    bool graph(vector<unsigned int> &ids, vector<unsigned int> &edges, unsigned int &id_count, unsigned int max_distance);
    
private:
    DBClientConnection conn;
    string addr;
    string port;
    string collection;
    
    unsigned int _size;
};
#endif
