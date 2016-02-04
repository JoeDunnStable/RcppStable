#include "Stable.h"

void warning(string msg) {
  Rcpp::warning(msg);
}

void stop(string msg) {
  Rcpp::stop(msg);
}

// [[Rcpp::export]]
double test_g(double th, double x, double alpha, double beta) {
  g_class param(alpha, beta);
  param.set_x(x);
  Rcpp::Rcout<<param;
  return param.g(th);
}

// [[Rcpp::export]]
Rcpp::NumericVector test_sin(Rcpp::NumericVector x) {
  Rcpp::NumericVector ret(x.size());
  for (int i=0; i<x.size(); i++) ret[i]=sin(x[i]);
  return wrap(ret);
}
