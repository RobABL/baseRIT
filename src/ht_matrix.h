#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include "dataset.h"

// [[Rcpp::plugins(cpp11)]]

using namespace std;

class Ht_matrix{
  private:
    vector<vector<int> > mat;
    int nb_perm;
    int nb_attr;
    int nb_inst;
    
  public:
    int get_nb_perm() const;
    int get_nb_attr() const;
    int get_nb_inst() const;
    const vector<int>& operator[](size_t) const;
    Ht_matrix();
    Ht_matrix(Dataset const&, int, int);
};

#endif
