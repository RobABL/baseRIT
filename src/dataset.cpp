#include <Rcpp.h>
#include <vector>
#include "dataset.h"
#include<unordered_set>

// [[Rcpp::plugins(cpp11)]]

using namespace std;

size_t Dataset::size() const{
  return data.size();
}

const unordered_set<int>& Dataset::operator[](size_t idx) const{
  return data[idx];
}

Dataset::Dataset(Rcpp::DataFrame const& class_data){
  int ncol = class_data.size();
  int nrow = Rcpp::as<vector<int> >(class_data[0]).size();
  vector<unordered_set<int> > active_data(nrow);
  
  for(int r=0;r<nrow;++r){
    unordered_set<int> active;
    for(int c=0;c<ncol;++c){
      bool elem = Rcpp::as<vector<int> >(class_data[c])[r];
      if(elem != 0)
        active.insert(c);
    }
    active_data[r] = active;
  }
  
  data = active_data;
}

Dataset::Dataset(){
  // Do nothing
}
