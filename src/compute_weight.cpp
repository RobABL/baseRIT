#include <Rcpp.h>
#include <cmath>
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

// [[Rcpp::export]]
double compute_weight(DataFrame const& instance,IntegerVector const& interaction, List const& map){
   
   for(int i=0;i<interaction.size();++i){
     
     int disc_attr_idx = interaction[i];
     --disc_attr_idx;
     List corresp = as<List>(map[disc_attr_idx]);
     int orig_attr_idx = corresp["attr"];
     --orig_attr_idx;
     NumericVector value = as<NumericVector>(corresp["value"]);

     if(value.size() > 1){ // Discretized continuous attribute
       double variance = corresp["var"];
       double instance_value = instance[orig_attr_idx];
       double lb = value[0];
       double ub = value[1];
       
       if(instance_value < lb || instance_value > ub){ // value lower than lower bound or higher than higher bound
          return 0;
       }
       // If value is in interval, continue
     }
     else{ // Categorical attribute
       int instance_value = instance[orig_attr_idx];
       int v = value[0];
       if(instance_value != v)
        return 0;
     }
   }
   
   return 1;
}
