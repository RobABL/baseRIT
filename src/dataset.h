#ifndef DATASET_H
#define DATASET_H

#include <Rcpp.h>
#include <vector>
#include <unordered_set>

// [[Rcpp::plugins(cpp11)]]

using namespace std;

class Dataset{
  private:
    vector<unordered_set<int> > data;
  
  public:
    size_t size() const;
    const unordered_set<int>& operator[](size_t idx) const;
    Dataset(Rcpp::DataFrame const& class_data);
    Dataset();
};

#endif
