#include <string>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <map>
#include <set>
#include <bitset>
#include <pthread.h>
#include <fstream>
#include <climits>
#include <random>
#include <algorithm>
#include <unordered_set>
#include <unordered_map>
#include <boost/threadpool.hpp>
#include <cstring>
#include <unistd.h>
#include <stdio.h>
#include <sys/types.h>

#define NUM_NODES 1000000
#define BUFFERSIZE 512
#define THREADS 32

//#define NUMBER_WALKS 500
//#define LENGTH_WALKS 10


using namespace std;
using namespace boost::threadpool;


struct Edge {
  unsigned int edge ;
  unsigned int node ;
} ;

unordered_map<unsigned int, vector<Edge>> graph ;
unordered_map<unsigned int, string> graph_dict;
vector<string> phenos_vec;
vector<string> nonphenos_vec;

random_device rd;
mt19937 rng(rd());
uniform_int_distribution<int> uni(0,INT_MAX);


ofstream fout;
boost::mutex mtx;


void build_graph(string fname) {
  char buffer[BUFFERSIZE];
  graph.reserve(NUM_NODES) ;
  ifstream in(fname);
  while(in) {
    in.getline(buffer, BUFFERSIZE);
    if(in) {
      Edge e ;
      unsigned int source = atoi(strtok(buffer, " "));
      e.node = atoi(strtok(NULL, " ")) ;
      e.edge = atoi(strtok(NULL, " ")) ;
      graph[source].push_back(e) ;
    }
  }
}


void get_graphdict(string fname) {
  char buffer[BUFFERSIZE];
  graph.reserve(NUM_NODES) ;
  ifstream in(fname);
  while(in) {
    in.getline(buffer, BUFFERSIZE);
    if(in) {
      string node_str = strtok(buffer, " ");
      int id = atoi(strtok(NULL, " "));
      graph_dict[id] = node_str;
    }
  }
}


void get_phenos_vec(string fname) {
  char buffer[BUFFERSIZE];
  graph.reserve(NUM_NODES) ;
  ifstream in(fname);
  while(in) {
    in.getline(buffer, BUFFERSIZE);
    if(in) {
      string node_str = strtok(buffer, " ");
      phenos_vec.push_back(node_str);
    }
  }
}


void get_nonphenos_vec(string fname) {
  char buffer[BUFFERSIZE];
  graph.reserve(NUM_NODES) ;
  ifstream in(fname);
  while(in) {
    in.getline(buffer, BUFFERSIZE);
    if(in) {
      string node_str = strtok(buffer, " ");
      nonphenos_vec.push_back(node_str);
    }
  }
}


void walk(unsigned int source, unsigned int num_walks, unsigned int length_walks) {
  vector<vector<unsigned int>> walks(num_walks);
  vector<Edge> adj_phenos;
  vector<Edge> adj_phenos_genes;

  if (graph[source].size()>0) { // if there are outgoing edges at all
    for (unsigned int i = 0 ; i < num_walks ; i++) {
      int count = length_walks;
      int current = source ;

      walks[i].push_back(source);
      while (count > 0) {
	if (graph[current].size() > 0 ){ // if there are outgoing edges
	  vector<Edge> neigbors = graph[current];
	  
	  //get phenotypes neigbor nodes
          for (unsigned int i = 0; i < neigbors.size(); i++) {
	      int neigbor = neigbors[i].node;
	      if (graph_dict[neigbor].find("MP_") == 0){
                  adj_phenos.push_back(neigbors[i]);
	      }
	  }

	  Edge next;
          unsigned int r;
          int target;
          int edge;

	  //A: if node is gene with phenotypes or pheno node, sample only from phenos neigbor
          if (adj_phenos.size() > 0 && find(phenos_vec.begin(), phenos_vec.end(),  
          graph_dict[current]) != phenos_vec.end() && graph_dict[current].find("MP_") == 0){

              r = uni(rng) % adj_phenos.size();
	      next = adj_phenos[r];
	  } 

	/*
	 //B: if node is gene with no phenotypes, one-hop sample from adj genes with phenos
	 else if (find(nonphenos_vec.begin(), nonphenos_vec.end(),  
          graph_dict[current]) != nonphenos_vec.end()){
	     for (unsigned int i=0; i < neigbors.size(); i++){
		 int neigbor = neigbors[i].node; 
	         if (find(phenos_vec.begin(), phenos_vec.end(),
		     graph_dict[neigbor]) != phenos_vec.end()){
		     adj_phenos_genes.push_back(neigbors[i]);	
	         }     
	     }
	     if (adj_phenos_genes.size() > 0){
	         r = uni(rng) % adj_phenos_genes.size();
		 next = adj_phenos_genes[i];    	 
	     }
	}*/


	//C: sample from any neigboring node
	else{ 
	     r = uni(rng) % graph[current].size();	
	     next = graph[current][r];
 	}

	
	mtx.lock(); //prevent threads from doing weird things
	target = next.node;
	edge = next.edge;
	walks[i].push_back(edge);
	walks[i].push_back(target);
	current = target;
 	mtx.unlock();
	}

        else {
          int edge = INT_MAX ; // null edge
          current = source ;
          walks[i].push_back(edge) ;
          walks[i].push_back(current) ;
	}
	count--;
      }
    }
  }

  stringstream ss;
  for(vector<vector<unsigned int>>::iterator it = walks.begin(); it != walks.end(); ++it) {
    for(size_t i = 0; i < (*it).size(); ++i) {
      if(i != 0) {
	ss << " ";
      }
      ss << (*it)[i];
    }
    ss << "\n" ;
  }
  mtx.lock() ;
  fout << ss.str() ;
  fout.flush() ;
  mtx.unlock() ;
}

void generate_corpus(unsigned int num_walks, unsigned int length_walks) {
  pool tp(THREADS);
  for ( auto it = graph.begin(); it != graph.end(); ++it ) {
    unsigned int source = it -> first ;
    tp.schedule(boost::bind(&walk, source, num_walks, length_walks)) ;
  }
  cout << tp.pending() << " tasks pending." << "\n" ;
  tp.wait() ;
}

int main (int argc, char *argv[]) {
  cout << "Building graph from " << argv[1] << "....\n";
  cout << "Building graph dict from " << argv[2] << "....\n";
  cout << "Reading genes with phenotypes..." << argv[3] << "\n";
  cout << "Reading genes with no phenotypes..." << argv[4] << "\n";
  get_graphdict(argv[2]); 
  get_phenos_vec(argv[3]);
  get_nonphenos_vec(argv[4]);
 
  unsigned int num_walks = atoi(argv[5]);
  unsigned int length_walks = atoi(argv[6]);
  build_graph(argv[1]);
  cout << "Number of nodes in graph: " << graph.size() << "\n" ;
  cout << "Writing walks to " << argv[7] << "\n" ;
  
  fout.open(argv[7]);
  generate_corpus(num_walks, length_walks);
  fout.close();
}

