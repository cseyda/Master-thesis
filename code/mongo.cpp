#include <string>
#include <vector>

#include "dbclient.h"
#include "node.h"
#include "mongo.h"

using namespace std;
using namespace mongo;


DB::DB() : addr("127.0.0.1:"), port("27017") {
    string errmsg;
    if ( ! conn.connect( addr + port , errmsg ) ) {
        cout << "couldn't connect : " << errmsg << endl;
        throw;
    }
    conn.setWriteConcern(W_NONE);
};

DB::~DB(){};

void DB::nodes(vector<Node*> &vecNodes, string collection){
    const string db = "tweets";
    unique_ptr<DBClientCursor> cursor = conn.query(db + "." + collection, BSONObj());
    unsigned int count = conn.count(db + "." + collection, BSONObj());
    vecNodes.resize(count);
    
    while( cursor->more() ) {
        BSONObj obj = cursor->next();
        
        unsigned int i = obj["_id"].Int();
        
        vector<BSONElement> vec = obj["loc"].Array();
        double x = vec[0].Double();
        double y = vec[1].Double();
        
        vec = obj["tag"].Array();
        string tag = " ";
        
        Node* pt = new Node(x, y, i);
        for (auto it = vec.begin() ; it != vec.end(); ++it){
            tag += it->String() + " ";
            pt->tags.push_back(it->String());
        }
        
        
        pt->tag = tag;
        
        //vecNodes.push_back(pt);
        vecNodes[i] = pt;
    }
};

void DB::bigrams(vector<string> &bigram_vec, string collection){
    const string db = "tweets";
    unique_ptr<DBClientCursor> cursor = conn.query(db + "." + "stats", QUERY("_id" << collection));

    while( cursor->more() ) {
        mongo::BSONObj obj = cursor->next();
        
        obj.getObjectField("bigrams").Vals(bigram_vec);
    }
};

void DB::jaccard(vector<vector<unsigned int>> &jac_vec, string collection){
    const string db = "tweets";
    unique_ptr<DBClientCursor> cursor = conn.query(db + "." + collection, BSONObj());
    while( cursor->more() ) {
//bigram
//words
        BSONObj obj = cursor->next();
        //vector<BSONElement> vec = obj["bigram"].Array();
        vector<int> tmp_mongo;
        obj.getObjectField("bigram").Vals(tmp_mongo);
        
        vector<unsigned int> tmp_vec;
        for (auto val : tmp_mongo){
            tmp_vec.push_back(val);
        }
        
        jac_vec.push_back(tmp_vec);
    }
};



//
// DB2
//

DB2::DB2(string collection): 
            addr("127.0.0.1:"), 
            port("27017"),
            collection(collection) {
    string errmsg;
    if ( ! conn.connect( addr + port , errmsg ) ) {
        cout << "couldn't connect : " << errmsg << endl;
        throw;
    }
    conn.setWriteConcern(W_NONE);
    
    const string db = "tweets";
    _size = conn.count(db + "." + collection, BSONObj());
};

DB2::~DB2(){};
    
unsigned int DB2::size(){
    return _size;
};

string DB2::getCollection(){
    return collection;
};
    
void DB2::get(unsigned int id, double& longitude, double& latitude, vector<string>& words){
    const string db = "tweets";
    unique_ptr<DBClientCursor> cursor = conn.query(db + "." + collection, QUERY("_id" << id));

    while( cursor->more() ) {
        mongo::BSONObj obj = cursor->next();
        
        obj.getObjectField("tag").Vals(words);
        //auto w_tmp = obj["tag"].Array();
        //for (unsigned int i=0; i<w_tmp.size(); i++){
        //    words[i] += w_tmp[i].Int();
        //}
        
        vector<BSONElement> vec = obj["loc"].Array();
        longitude = vec[0].Double();
        latitude  = vec[1].Double();
    }
};

void DB2::get(unsigned int id, double& longitude, double& latitude, vector<unsigned int>& words){
    const string db = "tweets";
    unique_ptr<DBClientCursor> cursor = conn.query(db + "." + collection, QUERY("_id" << id));

    while( cursor->more() ) {
        mongo::BSONObj obj = cursor->next();
        
        //obj.getObjectField("words").Vals(words);
        
        vector<int> tmp_mongo;
        obj["words"].Obj().Vals(tmp_mongo);
        
        //vector<unsigned int> tmp_vec;
        for (auto val : tmp_mongo){
            words.push_back(val);
        }
        
        
        vector<BSONElement> vec = obj["loc"].Array();
        longitude = vec[0].Double();
        latitude  = vec[1].Double();
    }
};

void DB2::jaccard(vector<vector<unsigned int>> &jac_vec, string kind){
    const string db = "tweets";
    unique_ptr<DBClientCursor> cursor = conn.query(db + "." + collection, BSONObj());
    
    while( cursor->more() ) {
//bigram
//words
        BSONObj obj = cursor->next();
        //vector<BSONElement> vec = obj["bigram"].Array();
        vector<int> tmp_mongo;
        
        //obj.getObjectField(kind).Vals(tmp_mongo);
        obj[kind].Obj().Vals(tmp_mongo);
        
        vector<unsigned int> tmp_vec;
        for (auto val : tmp_mongo){
            tmp_vec.push_back(val);
        }
        
        jac_vec.push_back(tmp_vec);
    }
};

bool DB2::graph(vector<unsigned int> &ids, vector<unsigned int> &edges, unsigned int &id_count, unsigned int max_distance){
    const string db = "graphs";
    unsigned int c = conn.count(db + ".graph", QUERY("_id" << collection).obj);
    
    if (c == 1){
        unique_ptr<DBClientCursor> cursor = conn.query(db + ".graph", QUERY("_id" << collection));
    
        while( cursor->more() ) {
            BSONObj obj = cursor->next();
            vector<int> tmp_mongo;
            obj["orig2triang"].Obj().Vals(tmp_mongo);
            
            //vector<unsigned int> tmp_vec;
            for (auto val : tmp_mongo){
                ids.push_back(val);
            }
            id_count = obj["combined_points"].Int();
        }
        
        cursor = conn.query(db + "." + collection, BSONObj());
        
        while( cursor->more() ) {
            BSONObj obj = cursor->next();
            
            if (static_cast<unsigned int>(obj["dist"].Int()) < max_distance){
                vector<BSONElement> vec = obj["edge"].Array();
                
                edges.push_back(vec[0].Int());
                edges.push_back(vec[1].Int());
            }
        }
        
        return true;
    
    } else {
        id_count = size();
        for(unsigned int i=0; i<id_count; i++){
            ids.push_back(i);
        }
        
        return false;
    }
};
