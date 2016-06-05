#include <Rcpp.h>
#include <cmath>
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

// [[Rcpp::export]]
double compute_weight(DataFrame const& instance,IntegerVector const& interaction, List const& map,double radius){
   double dist = 0;
   
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
       
       if(instance_value < lb){ // value lower than lower bound
          if(radius != 0)
            dist = std::max(((instance_value - lb)*(instance_value - lb))/variance,dist);
          else
            return 0;
       }
       if(instance_value > ub){ // value higher than higher bound
          if(radius != 0)
            dist = std::max(((instance_value - ub)*(instance_value - ub))/variance,dist);
          else
            return 0;
       }
       // If value is in interval, distance is zero for the current attribute
     }
     else{ // Categorical attribute
       int instance_value = instance[orig_attr_idx];
       int v = value[0];
       if(instance_value != v)
        return 0;
     }
   }
   
   if(radius != 0){
     double weight = std::exp(-dist/radius);
     return weight;
   }
   else{
     return 1;
   }
}
