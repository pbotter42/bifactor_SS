#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]

NumericMatrix calculateInformation(List traceLines,
                                   NumericVector slope,
                                   NumericVector qPoints) {
  NumericMatrix out(traceLines.size(), qPoints.size());
  for(int i = 0; i < traceLines.size(); ++i) {
    NumericVector temp_q = as<List>(as<List>(as<List>(traceLines)[i])[0])[0];
    NumericVector temp_p = as<List>(as<List>(as<List>(traceLines)[i])[1])[0];   
    //out.row(i) =  (slope[i]*slope[i]) * qPoints * temp_p * temp_q;
    out.row(i) =  (slope[i]*slope[i]) * temp_p * temp_q;
  }
  return out;
}
/*
1. by doing sommething like a cumulative product, i need to be able to calculate information
 for likert items
2. also there may be a float point precision issue bc things are alittle for SOME of the qpoints when i do by hand

 
 */