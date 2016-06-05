#ifndef INTERACT_H
#define INTERACT_H

#include <vector>
#include <string>
#include <Rcpp.h>
#include "dataset.h"
#include "ht_matrix.h"
#include <unordered_set>

// [[Rcpp::plugins(cpp11)]]

using namespace std;

class Interaction{
    private:
        unordered_set<int> items;
        int depth;
        vector<double> prev;
        bool prev_valid;
        void compute_prev(vector<Ht_matrix> const&,vector<Dataset> const&);

    public:
        int size() const;
        void intersect(unordered_set<int> const&);
        bool check_prev(vector<Ht_matrix> const&,vector<double> const&,vector<Dataset> const&,bool);
        bool check_for_map(vector<Ht_matrix> const&,vector<double> const&,vector<Dataset> const&);
        string as_string() const;
        Rcpp::List as_List(vector<Ht_matrix> const&,vector<Dataset> const&);
        int get_depth() const;
        void set_depth(int);
        Interaction();
        Interaction(unordered_set<int> const&, int, int);
};

#endif
