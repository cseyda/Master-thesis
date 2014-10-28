#include <vector>

#include "node.h"
#include "metric.h"
#include "dbscan.h"
#include "mongo.h"
#include "hist.h"
#include "eval.h"
#include "graph.h"
#include "clustering.h"

#include <omp.h>
#include <algorithm>

//srand rand
#include <cstdlib>
#include <time.h>
#include <utility>

//#include <limits>

using namespace std;

int printhelp(){
    cout << "scan metric threshold [collection minPts]\n";
    cout << "     metric: 0 -> Loc\n";
    cout << "             1 -> Euclid\n";
    return 0;
};


int main(int argc, char* argv[]){
    // C M D
    // 1. metric || alpha/beta
    // 2. threshold
    //optional
    // 3. collection
    // 4. minPts
    
    
    //const int min_int = std::numeric_limits<int>::min();
    /*const unsigned int max_uint = std::numeric_limits<unsigned int>::max();
    cout << "max_uint: " << max_uint << endl;
    const unsigned long int max_ulint = std::numeric_limits<unsigned long int>::max();
    cout << "max_ulint: " << max_ulint << endl;
    const unsigned long long int max_ullint = std::numeric_limits<unsigned long long int>::max();
    cout << "max_ullint: " << max_ullint << endl;*/
    
    //int metric;
    //double eps;
    //string collection;
    unsigned int minPts = 8;
    
    string choise;
    string collection;
    string clustering_file;
    
    //for(int p=0; p<argc; p++){
    //    cout << p << ": " << argv[p] << endl;
    //}
        
    if (argc < 3){
        choise = "scan";
        collection = "a56c48c3df72452177dce28efd790ddc";//big
        //collection = "eec06569c0a325988dc790db3877e8ae";//small
        //return printhelp();
    } else {
        choise = argv[1];
        collection = argv[2];
        
        //metric = stoi(argv[1]);
        //eps = stod(argv[2]);
        
        //if (argc == 5){
        //    collection = argv[3];
        //    minPts = stoi(argv[4]);
        //} else {
        //    collection = "a56c48c3df72452177dce28efd790ddc";//big
            //collection = "eec06569c0a325988dc790db3877e8ae";//small
        //    minPts = 8;
        //}
    }
    
    //DB db;
    
    //vector<string> bigram_vec;
    //db.bigrams(bigram_vec, collection);
    
    //cout << vecNodes[47064]->GetLongitude() << " " << vecNodes[47064]->GetLatitude() << " " << vecNodes[47064]->index << endl;
    //cout << vecNodes[8154]->GetLongitude() << " " << vecNodes[8154]->GetLatitude() << " " << vecNodes[8154]->index << endl;
    
    DB2 db2(collection);
    
    if ("scan" == choise){
        //1 LocF: threshold, minPts
        //2 J:    threshold, minPts, type (bigram, word)
        
        //3 Random Walk:
        //4 Comb:
        cout << "scan" << endl;
        double threshold;
        unsigned int metric = stoi(argv[4]);
        string output = argv[3];
        string kind;
        
        vector<Met*> metrics;
        
        double w_1 = 1.0;
        double c = 0.1;
        unsigned int dist_limit = 0;
        DBSCAN* scan = 0;

        switch(metric) {
            case 1:
                threshold = stod(argv[5]);
                minPts = stoi(argv[6]);
                
                metrics.push_back(new LocF(db2, threshold, minPts));
                scan = new DBSCAN(*(metrics[0]));
                break;
            case 2:
                threshold = stod(argv[5]);
                minPts = stoi(argv[6]);
                kind = argv[7];
                
                if ("w"== kind){
                    kind = "words";
                } else {
                    kind = "bigram";
                }
                cout << kind << endl;
                //bigram
                //words
                metrics.push_back(new JaccardIndex(db2, threshold, minPts, kind));
                scan = new DBSCAN(*(metrics[0]));
                
                //JaccardIndex j(db2, threshold, minPts, kind);
                //scan = new DBSCAN(j);
                break;
            case 3:{
                threshold = stod(argv[5]);
                minPts = stoi(argv[6]);
                
                w_1 = stod(argv[7]);
                c = stod(argv[8]);
                dist_limit = stoi(argv[9]);
                unsigned int w_threshold = stoi(argv[10]);
                //bool compressed = compressed;
                metrics.push_back(new GraphMetric(db2, threshold, minPts, w_1, c, dist_limit, w_threshold, true));
                scan = new DBSCAN(*(metrics[0]));
                
                break;}
            case 4:                
                //LocF
                threshold = stod(argv[8]);
                minPts = stoi(argv[9]);
                
                metrics.push_back(new LocF(db2, threshold, minPts));
                
                //Jaccard
                threshold = stod(argv[11]);
                minPts = stoi(argv[12]);
                kind = argv[13];
                if ("w"== kind){
                    kind = "words";
                } else {
                    kind = "bigram";
                }
                
                metrics.push_back(new JaccardIndex(db2, threshold, minPts, kind, false));
                
                threshold = stod(argv[5]);
                minPts = stoi(argv[6]);
                
                scan = new DBSCAN(*(new CombinedMetric(db2, &metrics, threshold, minPts)));
            
                break;
            
            default:
                break;
        }
        
        //GraphMetric comb(vecNodes);
        
        cout << "run scan" << endl;
        Clustering clustering;
        /*
        DBSCAN scan(comb);
        cout << "run\n";
        cout << scan.run(clustering) << endl;
        */
        unsigned int clusters = scan->run(clustering);
        
        cout << "clusters " << clusters-1 << endl;
        cout << "clustered " << db2.size()-clustering.clustering[0].size() << endl;
        cout << "unclustered " << clustering.clustering[0].size() << endl;
        
        clustering.write_file(output, true);
        
    } else if ("eval" == choise) {
        clustering_file = argv[3];
        unsigned int sample_make = stoi(argv[4]);
        unsigned int sample_skip = stoi(argv[5]);
        
        bool count = false;
        if ("1" == string(argv[6])) {
            count = true;
        }
        LocF LF(db2, 10000.0, minPts);
        
        
        //CombinedMetric comb(&metrics);
        Clustering cluster;
        cluster.read_file(clustering_file, sample_make, sample_skip);
        
        sampling(cluster, sample_make, sample_make);
        Eval e(&LF, cluster);
        double dunn, db, c, sw;
        
        cout << "eval start" << endl;
        e.all(dunn, db, c, sw, count);//, db2.size());
        cout << "Eval Geo D " << dunn << endl;
        cout << "Eval Geo DB " << db << endl;
        cout << "Eval Geo C " << c << endl;
        cout << "Eval Geo SW " << sw << endl;
        
        JaccardIndex JI(db2, 0.5, minPts, "bigram", false);
        e.setMetric(&JI);
        e.all(dunn, db, c, sw, count);//, db2.size());
        cout << "Eval JBi D " << dunn << endl;
        cout << "Eval JBi DB " << db << endl;
        cout << "Eval JBi C " << c << endl;
        cout << "Eval JBi SW " << sw << endl;
        
    } else if ("stat" == choise) {
        cout << "Stat\n";
        Met *j = new JaccardIndex(db2, 0.5, minPts, "bigram");
        //Met *j = new LocF(vecNodes, 15000.0, minPts, collection);
        Hist h(j, /*bins*/4, 0, 1);
        
        vector<unsigned int> counts;
        vector<double> limits;
        h.hist(counts, limits, db2.size());
        
        for (auto it = counts.begin(); it != counts.end(); it++){
            cout << *it << endl;
        }
    } else if ("graph" == choise){
        vector<Met*> metrics;
        metrics.push_back(new LocF(db2, 10000.0, minPts));
        metrics.push_back(new JaccardIndex(db2, 0.5, minPts, "bigram"));
        
        CombinedMetric comb(db2, &metrics, 1.0, minPts);
        
        LemonGraph handler(db2.size());
        //SimpleGraph handler(vecNodes.size());
        //Graph graph(&comb, handler);
        //Graph2 graph(vecNodes, handler, 0);
        //handler.loadGraph("graph.lemon");
        //handler.saveGraph("graph.lemon");
        //handler.initAugmentedGraph();
        
        /*
        for (unsigned int k=1; k<20; k++){
            unsigned int sum = 0;
            unsigned int max=0;
            unsigned int min=999999;
            
            unsigned int res;
            for (unsigned i=0; i<vecNodes.size(); i++){
                res = handler.bfs(i, k);
                sum+=res;
                if (res > max){max = res;}
                if (res < min){min = res;}
            }
            cout << "k: " << k << ", avg: " << (double) sum/vecNodes.size() << " min,max: " << min << " " << max << endl;
        }
        */
        
        
        /**/
        cout << "Init AugGraph" << endl;
        //handler.initAugmentedGraph();
        cout << "Dist Test 4" << endl;
        cout << "Graph dist " << handler.distance(0, 3, 4) << endl;
        cout << "Dist Test 8" << endl;
        cout << "Graph dist " << handler.distance(0, 3, 8) << endl;
        cout << "Dist Test 12" << endl;
        cout << "Graph dist " << handler.distance(0, 3, 12) << endl;
        
        /*
        double max = 0.0;
        double sum = 0.0;
        
        #pragma omp parallel default(shared) num_threads(4)
        {
        #pragma omp for schedule(dynamic)
        for(unsigned int i=0; i<vecNodes.size()-1; i++){
            for(unsigned int j=i+1; j<vecNodes.size(); j++){
                double dist = handler.distance(i, j, 5);
                if (max < dist) {max = dist;}
                sum += dist;
            }
        }
        }
        cout << "max: " << max << endl;
        cout << "avg: " << (sum*2) / (vecNodes.size() * vecNodes.size()) << endl;
        */
    } else {
        cout << "wrong option\n";
    }
};
