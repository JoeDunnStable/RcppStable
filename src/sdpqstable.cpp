#include <RcppArmadillo.h>
#include "Stable.hpp"

// [[Rcpp::export]]
NumericVector sdstable(NumericVector x, double alpha, double beta, int log_flag,
                       double tol, double zeta_tol, int subdivisions, int verbose)
{
  BEGIN_RCPP
  arma::vec ax(x);
  arma::vec ret(ax.n_elem);
  arma::uword i;
  for (i=0; i<ax.n_elem; i++)
    ret(i)=sdstable1(x(i), alpha, beta, log_flag, tol, zeta_tol, subdivisions, verbose);
  return wrap(ret);
  END_RCPP;
}

// [[Rcpp::export]]
NumericVector spstable(NumericVector z, double alpha, double beta, int lower_tail, int log_p,
                       double dbltol, int subdivisions, int verbose) {
  BEGIN_RCPP
  arma::vec az(z);
  arma::vec ret(az.n_elem);
  arma::uword i;
  for (i=0; i<az.n_elem; i++)
    ret(i)=spstable1(az(i), alpha, beta, lower_tail, log_p,
                        dbltol, subdivisions, verbose);
  return wrap(ret);
  END_RCPP
}

// [[Rcpp::export]]
NumericVector sqstable(NumericVector p, double alpha, double beta, int lower_tail, int log_p,
                       double dbltol, double integ_tol, int subdivisions, int verbose) {
  BEGIN_RCPP
  arma::vec ap(p);
  arma::vec ret(ap.n_elem);
  arma::uword i;
  for (i=0; i<ap.n_elem; i++)
    ret(i)=sqstable1(ap(i), alpha, beta, lower_tail, log_p,
                     dbltol, integ_tol, subdivisions, verbose);
  return wrap(ret);
  END_RCPP
}
