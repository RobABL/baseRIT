#include <Rcpp.h>
#include "interaction.h"
#include "ht_matrix.h"
#include "dataset.h"
#include <vector>
#include <string>
#include <cmath>
#include <sstream>
#include <algorithm>
#include <unordered_set>

// [[Rcpp::plugins(cpp11)]]

using namespace std;

int Interaction::size() const{
    return items.size();
}

string Interaction::as_string() const{
  stringstream ss;
  
  int i = 0;
  for(int x : items){
    ss << x + 1; // Add 1 for R indexes
    if(i < size()-1) // Not last
      ss << " ";
      
    i++;
  }
  
  return ss.str();
}

Rcpp::List Interaction::as_List(vector<Ht_matrix> const& matrices,vector<Dataset> const& all_datasets){
  Rcpp::IntegerVector inter(size());
  
  // Interaction
  int i = 0;
  for(int x : items){
    inter[i] = x + 1; // Add 1 for R indexes
    i++;
  }
  
  // If prev not valid, compute them
  if(!prev_valid)
    compute_prev(matrices,all_datasets);
    
  // Add prevalences
  Rcpp::NumericVector p(prev.size());
  for(int c=0;c<prev.size();++c)
    p[c] = prev[c];
    
  // Result
  vector<string>info({"interaction","prevalence","depth"});
  Rcpp::List res(3);
  res[0] = inter;
  res[1] = p;
  res[2] = depth;
  res.names() = info;
  
  return res;
}

bool Interaction::check_for_map(vector<Ht_matrix> const& matrices,vector<double> const& theta,
vector<Dataset> const& all_datasets){
  bool low = check_prev(matrices,theta,all_datasets,true); // at least one low prevalence
  if(!low) // Optimization: all high
    return false;
    
  bool high = false; // at least one high prevalence
  
  for(int c=0;c<prev.size();++c){
    if(prev[c] > theta[c])
      high = true;
  }
  return low && high;
}

bool Interaction::check_prev(vector<Ht_matrix> const& matrices,vector<double> const& theta,
vector<Dataset> const& all_datasets,bool es){
  if(!es) // No early stopping
    return true;
    
  if(!prev_valid)
    compute_prev(matrices,all_datasets);
    
  for(int c=0;c<prev.size();++c){
    if(prev[c] < theta[c]) // Low prevalence in at least 1 class
      return true;
  }
  return false;
}

void Interaction::compute_prev(vector<Ht_matrix> const& matrices,vector<Dataset> const& all_datasets){
  int idx = 0; // class idx
  
  for(Ht_matrix const& ht : matrices){ // Iterate over 'classes'
    int L = ht.get_nb_perm();
    
    if(L > 0){ // L positive, use estimator
      double p1 = 0;
      double p2 = 0;
      int T;
      int m; // remember min_k H_{lk}
      for(int l=0;l<L;++l){
        T = 1;
        m = ht[*(items.begin())][l];
        for(int x : items){
          if(m != ht[x][l]){
            T = 0; // there is a pair k1,k2 for which ht[k1][l] != ht[k2][l]
            m = min(m,ht[x][l]); // min may have changed, check
          }
          if(m == 1 && T == 0)
            break;
        }
        
        p1 += T;
        p2 += m;
      }
      
      int n = ht.get_nb_inst();
      double pref = (n+1)/n;
      double recip = 1/(n+1);
      
      prev[idx] = pref*(double(L)/p2 - recip)*p1/double(L);
    }
    else{ // L negative or zero, compute prevs from dataset
    
      int n = all_datasets[idx].size();
      int count = 0;
      for(int i=0;i<n;++i){ //Iterate over instances in class idx
        bool matching = true;
        for(int x : items){
          if(all_datasets[idx][i].count(x) == 0){ // mismatch found
            matching = false;
            break;
          }
        }
        
        if(matching) // No mismatch found
          ++count; 
      }
      
      prev[idx] = double(count)/double(n);
    }
    
    // Next prevalence
    ++idx;
  }
  
  prev_valid = true;
}

void Interaction::intersect(unordered_set<int> const& instance){
    unordered_set<int> cpy(items);
    
    for(int x : cpy){
      if(instance.count(x) == 0){
        items.erase(x);
      }
    }
    
    // Check prev
    if(items.size() != cpy.size())
      prev_valid = false;
}

int Interaction::get_depth() const{
  return depth;
}

void Interaction::set_depth(int i){
  depth = i;
}

Interaction::Interaction(){
  depth = -1;
  prev_valid = false;
}

Interaction::Interaction(unordered_set<int> const& instance,int p_depth,int nb_class){
    items = instance;
    depth = p_depth;
    prev_valid = false;
    vector<double> new_prev(nb_class);
    prev = new_prev;
}
