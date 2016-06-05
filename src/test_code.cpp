#include <Rcpp.h>
#include "dataset.h"
#include "interaction.h"
#include <unordered_set>
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

// [[Rcpp::export]]
void print_test(DataFrame const& data){
  Dataset active(data);
  for(int r=0;r<active.size();++r){
    std::unordered_set<int> instance = active[r];
    for(int x : instance){
      Rcout << x+1 << " ";
    }
    Rcout << std::endl;
  }
}

// [[Rcpp::export]]
void print_test2(DataFrame const& data,int idx0,int idx1,int nb_class){
  Dataset active(data);
  std::unordered_set<int> instance0 = active[idx0];
  std::unordered_set<int> instance1 = active[idx1];
  
  Interaction inter(instance0,0,nb_class);
  inter.intersect(instance1);
  Rcout << inter.as_string();
}