#include <Rcpp.h>
#include "interaction.h"
#include "ht_matrix.h"
#include "dataset.h"
#include <vector>
#include <string>
#include <chrono>
#include <random>
#include <functional>
#include <stack>
#include <unordered_map>
#include <unordered_set>
#include <sys/types.h>
#include <unistd.h>


// [[Rcpp::plugins(cpp11)]]

using namespace std;
typedef _Bind<uniform_int_distribution<int>(mersenne_twister_engine<long unsigned int, 32ul, 624ul, 397ul, 31ul, 2567483615ul, 11ul, 4294967295ul, 7ul, 2636928640ul, 15ul, 4022730752ul, 18ul, 1812433253ul>)> Generator;

Generator create_PRNG(int nrows){
  auto seed = chrono::high_resolution_clock::now().time_since_epoch().count()*::getpid();
  Generator int_rand = bind(uniform_int_distribution<int>(0,nrows-1),mt19937(seed));
  return int_rand;
}

void insert_interaction(unordered_map<string,Interaction>& map,Interaction inter,
vector<Ht_matrix> const& matrices,vector<double> const& theta, vector<Dataset> const& all_datasets){
  string repres = inter.as_string();
  
  if(map.count(repres) > 0) // element already in map, do nothing
    return;
    
  if(inter.check_for_map(matrices,theta,all_datasets))
    map[repres] = inter;
}

unordered_set<int> random_instance(Dataset const& x, Generator& prng){
  // Select random row number
  int row = prng();
  
  return x[row];
}

// DFS search on tree. Places results in leaves parameter.
void tree(unordered_map<string,Interaction>& leaves,vector<Dataset> const& all_datasets, int cls,
vector<Ht_matrix> const& all_matrices,vector<double> const& theta,int depth,int split_nb,
int min_inter_sz,int branch,bool es,Generator& prng){
  stack<Interaction, vector<Interaction> > frontier;
  
  // Init root
  int nb_class = all_matrices.size();
  unordered_set<int> root_instance = random_instance(all_datasets[cls],prng);
  Interaction root(root_instance,0,nb_class);
  frontier.push(root);

  while(!frontier.empty()){
    Interaction parent = frontier.top();
    frontier.pop();
    int c_depth = parent.get_depth();
    int eff_branch = (c_depth <= 0 || c_depth%split_nb != 0) ? 1 : branch;
    
    int nb_valid = 0; // number of valid children for this parent node
    for(int i=0;i<eff_branch;++i){
      Interaction child(parent);
      child.intersect(random_instance(all_datasets[cls],prng));
      child.set_depth(c_depth + 1);
      
      if(child.size() >= min_inter_sz && child.check_prev(all_matrices,theta,all_datasets,es)){ // child is valid
        if(child.size() == min_inter_sz || child.get_depth() == depth){ // child is leaf
          insert_interaction(leaves,child,all_matrices,theta,all_datasets);
        }
        else{ // child is valid but not leaf
          frontier.push(child); 
        }
        nb_valid++;
      }
    }
    
    if(nb_valid == 0) // If no child is valid, add parent as leaf
        insert_interaction(leaves,parent,all_matrices,theta,all_datasets);
  }
}

// [[Rcpp::export]]
Rcpp::List cpp_naive_RIT(Rcpp::List const& datas,Rcpp::NumericVector const& theta,
int n_trees,int depth,int split_nb,int branch,int min_inter_sz,int L,bool es){
  
  unordered_map<string,Interaction> res;
  
  // Transform types
  vector<double> c_theta = Rcpp::as<vector<double> >(theta);
  
  // Generate Ht matrices and Dataset objects
  vector<Dataset> all_datasets(datas.size());
  vector<Ht_matrix> all_matrices(datas.size());
  int p = Rcpp::as<Rcpp::DataFrame>(datas[0]).size();
  for(int c=0;c<datas.size();++c){
    Dataset class_data(Rcpp::as<Rcpp::DataFrame>(datas[c]));
    all_datasets[c] = class_data;
    
    Ht_matrix mat(class_data,L,p);
    all_matrices[c] = mat;
  }
 
  for(int c=0;c<datas.size();++c){// Iterate over classes
    //Rcpp::Rcout << "Start trees in class " << c << endl;
    Generator gen = create_PRNG(all_datasets[c].size());
    
    for(int t=0;t<n_trees;++t){// Iterate over n_trees
      tree(res,all_datasets,c,all_matrices,c_theta,depth*split_nb,split_nb,min_inter_sz,branch,es,gen);
      Rcpp::checkUserInterrupt();
    }
  }
  
  // Transform to R type
  Rcpp::List ret(res.size());
  Rcpp::CharacterVector keys(res.size());
  int i = 0;
  for(auto& kv : res){
    keys[i] = kv.first;
    ret[i] = kv.second.as_List(all_matrices,all_datasets);
    i++;
  }
  ret.names() = keys;
  
  return ret;
}