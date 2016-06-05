#include "ht_matrix.h"
#include <random>
#include "dataset.h"

// [[Rcpp::plugins(cpp11)]]

using namespace std;

int Ht_matrix::get_nb_perm() const{
  return nb_perm;
}

int Ht_matrix::get_nb_attr() const{
  return nb_attr;
}

int Ht_matrix::get_nb_inst() const{
  return nb_inst;
}

const vector<int>& Ht_matrix::operator[](size_t idx) const{
  return mat[idx];
}

Ht_matrix::Ht_matrix(){
  // Do nothing
}

// L is hyper-parameter, p is the total number of attributes after discretization
// n is the number of instances in current class
// class_data is the class data (active attribute indexes)
Ht_matrix::Ht_matrix(Dataset const& class_data,int L, int p) : nb_perm(L), nb_attr(p), nb_inst(class_data.size()) {
  if(L <= 0){ // No permutation, we will compute the prevalences without using the estimator (i.e. compute from datasets)
    return;
  }
  
  // Init empty
  vector<int> sub(L);
  vector<vector<int> > m(p,sub);
  
  // Dimensions
  int n = class_data.size();
  
  // Random generator
  random_device rd;
  mt19937 g(rd());
  
  // Indexes of instances
  vector<int> perm(n);
  for(int k=0;k<n;++k)
    perm[k] = k;
    
  for(int l=0;l<L;++l){
    shuffle(perm.begin(),perm.end(),g); // permutation
    for(int k=0;k<p;++k){
      int count = 0;
      int i = 0;
      while(count == 0 && i<n){
        count = class_data[perm[i]].count(k); // nb of times k appears in instance perm[i]
        i++;
      }
      
      if(i == n && count == 0)
        m[k][l] = n+1;
      else
        m[k][l] = i;
    }
  }
  
  mat = m;
}
