#include <vector>
#include <cmath>
#include <string>
#include <utility> // std::pair, std::make_pair
#include <set>
#include <unordered_set>
#include <unordered_map>
#include <map>
#include <fstream>
#include <queue>
#include <algorithm> //sort

#include "graph.h"
#include "dbscan.h"
//#include "node.h"
#include "metric.h"
#include "clustering.h"
#include "hist.h"

#include <lemon/list_graph.h>
//#include <lemon/concepts/graph.h>
#include <lemon/concepts/digraph.h>
#include <lemon/maps.h>
#include <lemon/path.h>


#include <Eigen/SparseCore>


#include <unsupported/Eigen/MatrixFunctions>

//
// CGAL http://www.cgal.org/Manual/latest/doc_html/cgal_manual/Triangulation_2/Chapter_main.html#Subsection_37.11.2
//
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned, K>    Vb;
typedef CGAL::Triangulation_data_structure_2<Vb>                    Tds;
typedef CGAL::Delaunay_triangulation_2<K, Tds>                      Delaunay;
typedef Delaunay::Point                                             Point;


template <typename T>
void FreeVector( T & t ) {
    T tmp;
    t.swap( tmp );
}



class SequenceGen {
public:
    SequenceGen (unsigned int start = 0) : current(start) { }
    unsigned int operator() () { return current++; }
private:
    unsigned int current;
};

class Comp{
public:
    Comp(vector<double>& v) : _v(&v) {}
    
    bool operator()(int i, int j){
        return (&(_v)[i] < &(_v[j]));
    }
private:
    vector<double>* _v;
};

Graph::Graph(CombinedMetric* metrics, GraphHandler &handler){
    //structure nodes
    vector<unsigned int> nodesFound;
    vector<double> distances;
    cout << "GraphMaker" << endl;
    bool many;
    unsigned int edges = 0;
    
    for(unsigned int i=0; i<metrics->getCount(); i++){
        many = metrics->Met::getNeighbors(i, nodesFound, distances); //??? try other
        
        vector<unsigned int> indices(nodesFound.size());
        std::generate(indices.begin(), indices.end(), SequenceGen(0));
        
        std::sort(indices.begin(), indices.end(), Comp(distances));
        

        if(many){
            for(unsigned int j=0; j<8; j++){
                handler.addEdge(i, nodesFound[indices[j]]);
                edges++;
            }
        } else {
            for(unsigned int j=0; j<nodesFound.size(); j++){
                handler.addEdge(i, nodesFound[j]);
                edges++;
            }
        }
        
        nodesFound.clear();
        distances.clear();
    }
    
    cout << edges << " Structure" << endl;
    
    //attribute nodes
    vector<Met*>* metrics_vec = metrics->metrics();
    for (unsigned int i=0; i<metrics_vec->size(); i++){
        
        DBSCAN scan(*(metrics_vec->at(i)));
        Clustering c;
        unsigned int clusters = scan.run(c);
        
        cout << metrics_vec->at(i)->getName() << " " << clusters  << endl;
        
        for(unsigned int cluster=0; cluster < clusters; cluster++){
            for(unsigned int k=0; k<c.clustering[cluster].size(); k++){
                unsigned int node = c.clustering[cluster][k];
                handler.addNodeAttribute(node, cluster);
            }
        }
        
        handler.addAttributeLimit(clusters);
        c.clustering.clear();
    }
    
    cout << "Attribute" << endl;
};

static const double DEG_TO_RAD = 0.017453292519943295769236907684886;
static const double EARTH_RADIUS_IN_METERS = 6372797.560856;

unsigned int diminishing_returns(unsigned int val) {
    double scale = 1.0;
    double mult = val / scale;
    double trinum = (sqrt(8.0 * mult + 1.0) - 1.0) / 2.0;
    return (unsigned int) trinum * scale;
};

inline double gd(double dist){
    return dist * EARTH_RADIUS_IN_METERS * DEG_TO_RAD;
    //return dist;
};
//geo base
//text augmented
Graph2::Graph2(DB2& db, vector<vector<pair<unsigned int, float>>> *scan, double w_1, double c, unsigned int dist_limit, unsigned int w_threshold, bool compressed){
    cout << "GraphMaker" << endl;
    double w_0 = 1.0;
    
    vector<unsigned int> nodesFound;
    vector<double> distances;
    
    //pygraph
    vector<unsigned int> id;
    vector<unsigned int> py_edges;
    unsigned int id_count;
    
    bool pygraph = db.graph(id, py_edges, id_count, dist_limit);
    
    if (compressed){
        if (!pygraph){
            throw "no graph saved";
        }
    };
    
    //CGAL
    std::vector< std::pair<Point,unsigned int> > cgal_points;
    
    // Graph test
    typedef Eigen::Triplet<float> T;
    std::vector<T> triplets;
    triplets.reserve(db.size()*10);
    
    //vector<unsigned int> edges;
    
    double min_val = 1.0;
    
    vector<vector<unsigned int>> words;
    words.resize(db.size());
    
    for(unsigned int i=0; i<db.size(); i++){
        //CGAL
        double lon, lat;
        db.get(i, lon, lat, words[i]);
        
        if (!pygraph){
            cgal_points.push_back(std::make_pair(Point(lon,lat), i));
        }
    }
    
    //
    // CGAL
    // start
    
    vector<vector<unsigned int>> cgal_edges(id_count);
    
    cout << "TransMap " << id_count << endl;
    
    vector<vector<unsigned int>> trans_map(id_count);
    for (unsigned int i=0; i<id.size(); i++){
        trans_map[id[i]].push_back(i);
    }
    
    //storage of edge_count for every real node for probability calc
    //[id][w_o / w_1]
    vector<vector<unsigned int>> edgecount(id_count, vector<unsigned int>(2));
    
    if (!compressed){
        cout << "Delaunay" << endl;
        
        LocF geo(db, 100/*eps*/, 8/*knn*/);
        Delaunay dt;
        
        dt.insert(cgal_points.begin(), cgal_points.end());
        
        for(Delaunay::Finite_edges_iterator it = dt.finite_edges_begin(); it != dt.finite_edges_end(); ++it){
            Delaunay::Edge e=*it;
    
            unsigned int i1 = e.first->vertex((e.second+1)%3)->info();
            unsigned int i2 = e.first->vertex((e.second+2)%3)->info();
            
            unsigned int dist = static_cast<unsigned int>(gd(geo.distance(trans_map[i1][0], trans_map[i2][0])));
        
            if (dist < dist_limit){  
                py_edges.push_back(i1);
                py_edges.push_back(i2);
            }
        }
    }
    
    
    
    cout << "counting" << endl;
    for (unsigned int i=0; i<py_edges.size(); i+=2){
        unsigned int m_i1 = py_edges[i];
        unsigned int m_i2 = py_edges[i+1];
        
        edgecount[m_i1][0]+=diminishing_returns(trans_map[m_i2].size());
        edgecount[m_i2][0]+=diminishing_returns(trans_map[m_i1].size());
        
        cgal_edges[m_i1].push_back(m_i2);
        cgal_edges[m_i2].push_back(m_i1);
    }
    
    //cout << "a" << endl;
    //
    // CGAL
    // end
    
    vector<set<unsigned int>> tags; // != words
    
    //
    // TEXT
    //start
    //cout << "b" << endl;
    
    if (w_1 != 0.0){
        for (unsigned int i=0; i<db.size(); i++){
            for(unsigned int tag_id : words[i]){
                if (tag_id >= tags.size()){
                    tags.resize(tag_id+1);
                }
                tags[tag_id].insert(i);
                
            }
        }
        //}
        //cout << "c" << endl;
        
        //
        //filter word count
        for (unsigned int i=0; i<tags.size(); i++){
            if (tags[i].size() <= w_threshold){
                for (unsigned int nid : tags[i]){
                    //deleting i from word vectors
        //http://stackoverflow.com/questions/3385229/
        //c-erase-vector-element-by-value-rather-than-by-position
                    auto it_end=std::remove(words[nid].begin(), words[nid].end(), i);
                    words[nid].erase(it_end, words[nid].end());
                }
                tags[i].clear();
            }
        }
        //cout << "d" << endl;
        
        
        // count aug edges
        //for (unsigned int i=0; i<db.size(); i++){
        //    edgecount[id[i]][1]+=words[i].size();
        //}
    } else {
        for (auto &it:words){
            it.clear();
        }
    }
    
    
    
    
    
    
    // [c_id][tag_id] = count
    vector<std::unordered_map<unsigned int, unsigned int>> tag_map(id_count);
    for (unsigned int w=0; w<words.size(); w++){
        for (auto tid: words[w]){
            tag_map[id[w]][tid]++;
        }
    }
    
    //DM on all tag counts
    for (auto &w : tag_map){
        for (auto &count : w){
            count.second = diminishing_returns(count.second);
        }
    }
    
    //edited edge count for augmented nodes
    vector<unsigned int> tags_size(tags.size(), 0);
    //edited edge count for every node
    for (unsigned int i=0; i<tag_map.size(); i++){
        for (auto count : tag_map[i]){
            edgecount[i][1] += count.second;
            tags_size[count.first] += count.second;
        }
    }
    
    //edited edge count for augmented nodes
    vector<unsigned int> tag_size(tags.size());
    
    
    
    std::unordered_map<unsigned int, float> to_text;
    std::unordered_map<unsigned int, float> to_node;
    //cout << "e" << endl;
    
    //calculating probabilities and collect results (many edges from one word to another)
    for (unsigned int c_id=0; c_id < trans_map.size(); c_id++){
        to_text.clear();
        to_node.clear();
        
        double w_sum = w_0*edgecount[c_id][0] + w_1*edgecount[c_id][1];
        double p_0 = w_0 / w_sum;
        double p_1 = w_1 / w_sum;
        //cout << w_sum << " " << p_0 << " " << p_1 << endl;
        //cout << edgecount[c_id][0] << " " << edgecount[c_id][1] << endl;
        
        
        //text
        for (auto count : tag_map[c_id]){
             to_text[id_count + count.first] += count.second * p_1;
             to_node[id_count + count.first] += count.second / tags_size[count.first];
        }
        
        
        
        //nodes
        for (unsigned int i2 : cgal_edges[c_id]){
            to_text[i2] += p_0 * diminishing_returns(trans_map[i2].size());
        }
    
        //submitting to matrix
        for (auto it=to_text.begin(); it!=to_text.end(); ++it){
            triplets.push_back(
                T(  it->first,
                    c_id,
                    it->second
                )
            );
            
            if (min_val > it->second){
                min_val = it->second;
            }
        }
        //cout << "i" << endl;
    
        for (auto it=to_node.begin(); it!=to_node.end(); ++it){
            triplets.push_back(
                T(  c_id,
                    it->first,
                    it->second
                )
            );
            
            if (min_val > it->second){
                min_val = it->second;
            }
        }
        //cout << "j" << endl;
    }
    //
    //
    
    
    //cout << "geo edges " << edge_count_sum << endl;
    //cout << "word edges " << triplets.size() - edge_count_sum << endl;
    
    //
    // TEXT
    //end
    
    
    
    cout << "Building DistMat" << endl;
    cout << "Min Val " << min_val << endl;
  
    typedef Eigen::SparseMatrix<float> Sparse;
    unsigned int matSize = id_count+tags.size()+1;
    Sparse P_A(matSize, matSize);
    Sparse P_1(matSize, matSize);
    Sparse R_A(matSize, matSize);
    
    
    P_1.setFromTriplets(triplets.begin(), triplets.end());
    
    FreeVector(triplets);
    
    P_A.setIdentity();
    /*
    for (unsigned int k=0; k<id_count; ++k){
        double val = 0.0;
        for (Sparse::InnerIterator it(P_1,k); it; ++it){
            val += it.value();
        }
        cout << k << ": " << val << " ";
    }
    cout << endl;
    */
    
    cout << "Multiplying DistMat" << endl;
    for(unsigned int i=1; i<20; i++){
        P_A = (P_A * P_1).pruned(min_val, 1);
        R_A += P_A * (c * pow(1.0-c, i));
        
        cout << "Iteration" << i << " done." << endl;
        cout << "NonZeros " << R_A.nonZeros() << endl;
    }
      
    
    //cout << "Max " << max << endl;
    //histo
    Hist h(0, /*bins*/200, 0.0, 1.0);
    
    vector<unsigned int> counts;
    vector<double> limits;
    h.prep_add(counts, limits);
    
    R_A += Sparse(R_A.transpose());
    R_A /= 2;
    
    for (unsigned int k=0; k<id_count; ++k){
        for (Sparse::InnerIterator it(R_A,k); it; ++it){
            if (it.row() < (int)id_count){
                double val = it.value();
                if (val > 1.0){
                    val = 1.0;
                }
                
                h.add(val, counts);
                
                if (scan != 0){
                    scan->at(k).push_back(make_pair(it.row(), 1.0-val));
                }
            }
        }   
    }
    
    
    for(unsigned int l=0; l<200; l++){
        if (counts[l] != 0){
            cout << limits[l] << "-" << limits[l+1] << " : " << counts[l] << endl;
        }
    }
    
    
    /**/
    Hist h_sort(0, 200, 0.0, 1.0);
    vector<unsigned int> counts_sort;
    vector<double> limits_sort;
    h_sort.prep_add(counts_sort, limits_sort);
    
    unsigned int nn = 4;
    cout << endl << nn << "-dist" << endl;
    
    vector<float> sort_me;
    for (unsigned int k=0; k<id_count; ++k){
        for (Sparse::InnerIterator it(R_A,k); it; ++it){
            if (it.row() < (int)id_count){
                if (it.value() > 1.0){
                    sort_me.push_back(1.0);
                } else {
                    sort_me.push_back(it.value());
                }
            }
        }

        std::sort(sort_me.begin(), sort_me.end());

        if (sort_me.size() >= nn){
            h_sort.add((double) sort_me[sort_me.size()-nn], counts_sort);
        }
        sort_me.clear();
    }

    for(unsigned int l=0; l<200; l++){
        if (counts_sort[l] != 0){
            cout << limits_sort[l] << "-" << limits_sort[l+1] << " : " << counts_sort[l] << endl;
        }
    }
};





// Prototype
// GraphHandler
GraphHandler::GraphHandler(unsigned int elements)
        : node_count(elements){};
GraphHandler::~GraphHandler(){};




//
// LemonGraph
//
LemonGraph::LemonGraph(unsigned int node_count)
        : GraphHandler(node_count),
        nodes(node_count),
        attributes(g),
        node_repr(g),
        arc_weights(g),
        node_degree(g),
        node_ids(g), visited(node_count, false) {
    
    g.reserveNode(node_count);
    for (unsigned int i=0; i<node_count; i++){
        nodes[i] = g.addNode();
    }
    cout << "LemonGraph" << endl;
};

LemonGraph::~LemonGraph(){};

void LemonGraph::addEdge(unsigned int i, unsigned int j){
    bool connected = false;
    for (OutArcIt a(g, nodes[i]); a!=lemon::INVALID; ++a){
        if (nodes[j] == g.target(a)) {;
            connected = true;
            break;
        }
    }
    if (!connected){
        g.addArc(nodes[i], nodes[j]);
        g.addArc(nodes[j], nodes[i]);
    }
};

void LemonGraph::addEdge(unsigned int i, unsigned int j, double weight){
    addEdge(i, j);
};

void LemonGraph::addNodeAttribute(unsigned int node, unsigned int attribute){
    attributes[nodes[node]].push_back(attribute);
};

void LemonGraph::addAttributeLimit(unsigned int limit){
    attrLimits.push_back(limit);
};


struct VectorReader {
    vector<unsigned int> operator()(const std::string& str) {
        //std::istringstream is(str);
        vector<unsigned int> value;
        
        unsigned int begin = 0;
        unsigned int end   = 0;
        
        while(true){
            end = str.find(",", begin);
            
            //if (end == string::npos){
            if(end == 4294967295){
                //cout << "end " <<  end << " " << str.substr(begin) << endl;
                value.push_back(stoul(str.substr(begin)));
                break;
            }
            //cout << "begin " << end << " " << str.substr(begin, end-begin) << endl;
            value.push_back(stoul(str.substr(begin, end-begin)));
            begin = end+1;
        }
        
        return value;
    }
};

struct VectorWriter {
    std::string operator()(const vector<unsigned int>& value) {
        string output = "";
        
        for(unsigned int i : value){
            output += to_string(i) + ",";
        }
        output.pop_back();
        
        return output;
    }
};

void LemonGraph::saveGraph(string location){
    cout << "Saving Graph" << endl;
    lemon::digraphWriter(g, location).
        // write cap into 'attributes'
        //NodeAttributeMap("attributes", attributes).
        nodeMap("attributes", attributes, VectorWriter()).
        attribute("attrLimits", attrLimits, VectorWriter()).
        run();
};

void LemonGraph::loadGraph(string location){
    cout << "Loading Graph" << endl;
    try {
    // read the directed graph into g
    lemon::digraphReader(g, location).
        // read the 'attributes' NodeAttributeMap into attributes
        //NodeAttributeMap("attributes", attributes).
        nodeMap("attributes", attributes, VectorReader()).
        attribute("attrLimits", attrLimits, VectorReader()).
        run();
        
        unsigned int i=0;
        for (NodeIt n(g); n!=lemon::INVALID; ++n){
            nodes[i] = n;
            i++;
        }
        
        cout << "Loading Graph done" << endl;

    } catch (lemon::Exception& error) { // check if there was any error
        std::cerr << "Error: " << error.what() << std::endl;
        cout << "Loading Graph failed" << endl;

        //return -1;
    }
};


void LemonGraph::initAugmentedGraph(){
    // create attribute node
    for (unsigned int i=0; i<attrLimits.size(); i++){
        g.reserveNode(attrLimits[i]);
        for(unsigned int j = 0; j<attrLimits[i]; j++){
            Node node = g.addNode();
            nodes.push_back(node);
            node_repr[node] = make_pair(i+1,j);
        }
    }
    
    // init nodes
    for (NodeIt n(g); n!=lemon::INVALID; ++n){
        node_repr[n] = make_pair(0,0);
        
        //create arcs
        unsigned int thres = 0;
        for(unsigned int i=0; i<attributes[n].size(); i++){
            
            unsigned int attr_node = thres + attributes[n][i];
            g.addArc(n, nodes[attr_node]);
            g.addArc(nodes[attr_node], n);
            thres += attrLimits[i];
        }
    }
    
    // init weights
    w.push_back(1.0);
    for(unsigned int i=0; i<attrLimits.size(); i++){
        w.push_back(1.0);
    }
    
    updateEdges();
};

//updateAll
void LemonGraph::updateAll(vector<vector<unsigned int>> &clustering, vector<unsigned int> &medians){
    updateWeights(clustering, medians);
    updateEdges();
};

double LemonGraph::distance(unsigned int source, unsigned int target, unsigned int length){
    if (length == 0){
        return 0.0;
    }
    
    set<Node> visited;
    unordered_map<unsigned int, double> hashed;
    
    //visited.insert(nodes[source]);
    
    double c = 0.2;
    vector<double> reset(length-1);
    for(unsigned int i=1; i<length; i++){
        reset[length-1-i] = c * pow((1-c), i);
    }
    
    double sum = 0.0;
    double atm = 0.0;
    
    for (OutArcIt a(g, nodes[source]); a!=lemon::INVALID; ++a){
        atm = arc_weights[a];// * reset[length-1];
        
        if (nodes[target] == g.target(a)){ 
            sum += atm * reset[length-1];
        } else {
            sum += distance_(Arc(a), nodes[target], length-1, atm, reset, visited, hashed);
        }
    }
    
    return sum;
};

//
// private
//
double LemonGraph::distance_(Arc a, Node target, unsigned int length, double atm, vector<double> &reset, set<Node> &visited, unordered_map<unsigned int, double> &hashed){
    //cout << atm << " " << length << " ";
    double sum = 0.0;
    Node arc_target = g.target(a);
    //unsigned int edges = 0;
    
    //string pair_now = to_string(node_ids[arc_target]) + " " + to_string(length);
    unsigned int pair_now = nodes.size() * length + node_ids[arc_target];
    
    auto got = hashed.find(pair_now);
    if (got == hashed.end()) {
        Node visit_node;
        
        for (OutArcIt a_new(g, arc_target); a_new!=lemon::INVALID; ++a_new){
            visit_node = g.target(a_new);
            
            if (target == visit_node){
                sum += atm * reset[length-1];
            } else { //recursion
                if (length != 1){//reached walking distance
                    if (visited.find(visit_node) == visited.end()){
                        visited.insert(visit_node);
                        double new_atm = atm * arc_weights[a_new];
                        sum += distance_(Arc(a_new), target, length-1, new_atm, reset, visited, hashed);
                    }
                }
            }
        }
        
        hashed[pair_now] = sum / atm;
    } else {
        sum = atm * got->second;
    }
    
    //cout << length << " " << edges << endl;
    //cout << sum << endl;
    visited.erase(arc_target);
    
    return sum;
};



//
// all distance
//
unsigned int LemonGraph::bfs(unsigned int source, unsigned int length){

    //set<Node> visited;
    //vector<bool> visited(nodes.size(), false);
    
    //vector<pair<Node, unsigned int>> distance;
    //queue<pair<Node, unsigned int>> todo;
    
    //visited.clear();
    for(auto iter=visited.begin(); iter!= visited.end(); iter++){
        *iter = false;
    }
    distance_vec.clear();
    //todo.clear();
    
    //cout << source << " " << length << endl;
    for (OutArcIt a(g, nodes[source]); a!=lemon::INVALID; ++a){
        //cout << "a " << in.second << endl;
        if (!visited[node_ids[g.target(a)]]){
            visited[node_ids[g.target(a)]] = true;
            todo.push(make_pair(g.target(a), 1));
            distance_vec.push_back(make_pair(g.target(a), 1));
        }
    }
    
    while(!todo.empty()){
        auto atm = todo.front();
        
        if (atm.second != length){
            for (OutArcIt a(g, atm.first); a!=lemon::INVALID; ++a){
                Node node = g.target(a);
                if (!visited[node_ids[g.target(a)]]){
                    visited[node_ids[g.target(a)]] = true;
                    todo.push(make_pair(node, atm.second+1));
                    if (node_repr[node].first == 0){
                        distance_vec.push_back(make_pair(node, atm.second+1));
                    }
                }
            }
        }
        todo.pop();
    }
    
    return distance_vec.size();
    
};


/*
k: 7, avg: 502.292 min,max: 0 4086
k: 8, avg: 502.787 min,max: 0 4086
k: 9, avg: 503.044 min,max: 0 4086
k: 10, avg: 503.157 min,max: 0 4086
k: 11, avg: 503.209 min,max: 0 4086
k: 12, avg: 503.222 min,max: 0 4086
k: 13, avg: 503.224 min,max: 0 4086
k: 14, avg: 503.225 min,max: 0 4086
k: 15, avg: 503.225 min,max: 0 4086
k: 16, avg: 503.225 min,max: 0 4086
k: 17, avg: 503.225 min,max: 0 4086
k: 18, avg: 503.225 min,max: 0 4086
k: 19, avg: 503.225 min,max: 0 4086 
 
*/


/*
augmented Graph1
all nodes reaching all other after 4 hops

k: 1, avg: 10.571 min,max: 2 45546
k: 2, avg: 45319.1 min,max: 27 57536
k: 3, avg: 57489.8 min,max: 36235 57554
k: 4, avg: 57554 min,max: 57554 57554  
*/



bool LemonGraph::getNeighbors(unsigned int pt, vector<unsigned int>& rgpNodesFound, vector<double>* distances, double eps, unsigned int minPts){
    return false;
};

//update Attribute Edges
void LemonGraph::updateEdges(){
    for (NodeIt source(g); source!=lemon::INVALID; ++source){
        for (OutArcIt a(g, source); a!=lemon::INVALID; ++a){
            Node target = g.target(a);
            
            arc_weights[a] = p(source, target);
        }
    }
};

// update w weights
void LemonGraph::updateWeights(vector<vector<unsigned int>> &clustering, vector<unsigned int> &medians){
    if (clustering.size() != medians.size()){
        cout << "error\n";
    }
    
    vector<unsigned int> votes(w.size(), 0);
    unsigned int votes_sum = 0;
    for (unsigned int i=0; i<w.size(); i++){
        for (unsigned int j=0; j<clustering.size(); j++){
            for (unsigned int k=0; k<clustering[j].size(); k++){
                votes[i] += vote(i, nodes[clustering[j][k]], nodes[medians[j]]);
            }
        }
        votes_sum += votes[i];
    }
    
    double w_sum = 0;
    for (unsigned int i=0; i<w.size(); i++){
        w[i] = 0.5 * (w[i] + (double (votes[i] * w.size())) / double (votes_sum));
        w_sum += w[i];
    }
    
    cout << "W: " << w_sum << " elements: " << w.size() << endl;
};

double LemonGraph::p(Node i, Node j){
    double w_j = w[node_repr[j].first];
    
    double w_s = 0.0;
    if (node_repr[i].first == 0){
        w_s = w.size()-1;
    }
    
    return w_j / (node_degree[i] + w_s);
};

unsigned int LemonGraph::vote(unsigned int attr, Node i, Node j){
    unsigned int vote = 0;
    if (attributes[i][attr] == attributes[j][attr]){
        vote = 1;
    }
    return vote;
};


//
// Simple Graph
//
SimpleGraph::SimpleGraph(unsigned int elements) :
        GraphHandler(elements), edges(elements) {
    //vector<vector<pair<unsigned int, double>>> edges
};

SimpleGraph::~SimpleGraph(){};
    
void SimpleGraph::addEdge(unsigned int i, unsigned int j){
    addEdge(i, j, 0);
};

void SimpleGraph::addEdge(unsigned int i, unsigned int j, double weight){
    edges[i].push_back(make_pair(j, weight));
    edges[j].push_back(make_pair(i, weight));
};

void SimpleGraph::addNodeAttribute(unsigned int node, unsigned int attribute){
    cout << "not implemented\n";
};

void SimpleGraph::addAttributeLimit(unsigned int limit){
    cout << "not implemented\n";
};

void SimpleGraph::saveGraph(string location){
    ofstream myfile;
    myfile.open(location);
    
    unsigned int node1 = 0;
    for(auto edge : edges ){
        for(auto node2 : edge){
            if (node1 < node2.first)
                myfile << node1 << " " << node2.first << " " << node2.second << endl;
        }
        node1++;
    }
    
    myfile.close();
};

void SimpleGraph::loadGraph(string location){
    cout << "not implemented\n";
};
