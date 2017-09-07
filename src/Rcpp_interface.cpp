#include "stable_distribution_Vec.h"
#include "stable_distribution_fit.h"

using namespace stable_distribution;

using Rcpp::NumericVector;
using Rcpp::IntegerVector;
using Rcpp::DoubleVector;
using Rcpp::CharacterVector;
using Rcpp::DataFrame;
using Rcpp::Named;
using Rcpp::wrap;

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
NumericVector dstable_cpp(Eigen::VectorXd x, double alpha, double beta,
                       Eigen::VectorXd gamma, Eigen::VectorXd delta,
                       int pm, int log_flag,
                       double tol, int subdivisions, int verbose)
{
  BEGIN_RCPP
  bool noext = true;  // No extrapolation.  It's too unreliable.
  int N = 10;  // Us a 10 point Gausss and 21 point Gauss Kronrod rule
  double epsabs=0.;
  double epsrel=tol;
  IntegrationController ctl(noext, N, epsabs, epsrel, subdivisions, verbose);
  return wrap(pdf(x, alpha, beta, gamma, delta, pm, log_flag, ctl, verbose));
  END_RCPP
}

// [[Rcpp::export]]
NumericVector dstable_quick(Eigen::VectorXd x, double alpha, double beta,
                            Eigen::VectorXd gamma, Eigen::VectorXd delta,
                            int pm, int log_flag,
                       double tol, int subdivisions, int verbose)
{
  BEGIN_RCPP
  bool noext = true;  // No extrapolation.  It's too unreliable.
  int N = 10;  // Us a 10 point Gausss and 21 point Gauss Kronrod rule
  double epsabs=0.;
  double epsrel=tol;
  IntegrationController ctl(noext, N, epsabs, epsrel, subdivisions, verbose);
  return wrap(pdf_quick(x, alpha, beta, gamma, delta, pm, log_flag, ctl, verbose));
  END_RCPP
}

// [[Rcpp::export]]
NumericVector pstable_cpp(Eigen::VectorXd z, double alpha, double beta,
                          Eigen::VectorXd gamma, Eigen::VectorXd delta,
                          int pm, int lower_tail, int log_p,
                          double dbltol, int subdivisions, int verbose) {
  BEGIN_RCPP
  bool noext = true;  // No extrapolation.  It's too unreliable.
  int N = 10;  // Us a 10 point Gausss and 21 point Gauss Kronrod rule
  double epsabs=0.;
  double epsrel=dbltol;
  IntegrationController ctl(noext, N, epsabs, epsrel, subdivisions, verbose);
  return wrap(cdf(z, alpha, beta, gamma, delta, pm, lower_tail, log_p, ctl, verbose));
  END_RCPP
}

// [[Rcpp::export]]
NumericVector qstable_cpp(Eigen::VectorXd p, double alpha, double beta,
                          Eigen::VectorXd gamma, Eigen::VectorXd delta,
                          int pm, int lower_tail, int log_p,
                          double dbltol,
                          double integ_tol, int subdivisions, int verbose) {
  BEGIN_RCPP
  bool noext = true;  // No extrapolation.  It's too unreliable.
  int N = 10;  // Use a 10 point Gausss and 21 point Gauss Kronrod rule
  double epsabs=0.;
  double epsrel=integ_tol;
  IntegrationController ctl(noext, N, epsabs, epsrel, subdivisions, verbose);
  return wrap(quantile(p, alpha, beta, gamma, delta, pm, lower_tail, log_p,
                     dbltol, ctl, verbose));
  END_RCPP
}

// [[Rcpp::export]]
NumericVector ddx_sdstable(Eigen::VectorXd x, double alpha, double beta,
                       double tol, int subdivisions, int verbose)
{
  BEGIN_RCPP
  bool noext = true;  // No extrapolation.  It's too unreliable.
  int N = 10;  // Us a 10 point Gausss and 21 point Gauss Kronrod rule
  double epsabs=0.;
  double epsrel=tol;
  IntegrationController ctl(noext, N, epsabs, epsrel, subdivisions, verbose);
  return wrap(std_ddx_pdf(x, alpha, beta, ctl, verbose));
  END_RCPP
}

// [[Rcpp::export]]
NumericVector rstable_cpp(double alpha, double beta,
                          Eigen::VectorXd gamma, Eigen::VectorXd delta,
                          int pm,
                          Eigen::VectorXd u1, Eigen::VectorXd u2)
{
  BEGIN_RCPP
  if (u1.size() != u2.size()) throw std::range_error("srstable: u1 & u2 sizes differ");
  size_t n = u1.size();
  // We need an integration controller to calculate the mode
  bool noext = true;  // No extrapolation.  It's too unreliable.
  int N = 10;  // Us a 10 point Gausss and 21 point Gauss Kronrod rule
  double epsabs=0.;
  double epsrel=64 * std::numeric_limits<double>::epsilon();
  int subdivisions = 1000;
  int verbose = 0;
  IntegrationController ctl(noext, N, epsabs, epsrel, subdivisions, verbose);
  Vec ret(n);
  ret=random_stable(alpha, beta, gamma, delta, pm, ctl, verbose, u1, u2);
  return wrap(ret);
  END_RCPP
}

// [[Rcpp::export]]
double sdstableMode(double alpha, double beta,
                    double dbltol,
                    double tol, int subdivisions, int verbose)
{
  bool noext = true;  // No extrapolation.  It's too unreliable.
  int N = 10;  // Us a 10 point Gausss and 21 point Gauss Kronrod rule
  double epsabs=0.;
  double epsrel=tol;
  IntegrationController ctl(noext, N, epsabs, epsrel, subdivisions, 0);
  StandardStableDistribution std_stable_dist(alpha, beta, ctl, 0);
  return std_stable_dist.mode(dbltol, verbose).first;
}

// [[Rcpp::export]]
DataFrame stable_fit_cpp(Eigen::VectorXd y, std::string type, bool quick) {
  BEGIN_RCPP
  bool noext = true;  // No extrapolation.  It's too unreliable.
  int N = 10;  // Us a 10 point Gausss and 21 point Gauss Kronrod rule
  double epsabs=0.;
  double epsrel=64*std::numeric_limits<double>::epsilon();
  int subdivisions = 1000;
  int verbose = 0;
  IntegrationController ctl(noext, N, epsabs, epsrel, subdivisions, verbose);
  double dbltol=1e-12;

  std::vector<FitResult> results = stable_fit(y, ctl, dbltol, type, quick, verbose);
  int m = results.size();
  DoubleVector alpha(m);
  DoubleVector beta(m);
  DoubleVector gamma(m);
  DoubleVector delta(m);
  DoubleVector two_ll_n(m);
  IntegerVector pm(m);
  IntegerVector n(m);
  CharacterVector method(m);
  DoubleVector q_kurt(m);
  DoubleVector q_skew(m);
  DoubleVector q_scale(m);
  DoubleVector q_location(m);
  CharacterVector convergence(m);
  IntegerVector iterations(m);
  DoubleVector cpu_time(m);
  for (int i = 0; i<m; ++i){
    alpha[i]=results.at(i).alpha;
    beta[i]=results.at(i).beta;
    gamma[i]=results.at(i).gamma;
    delta[i]=results.at(i).delta;
    two_ll_n[i]=results.at(i).two_ll_n;
    pm[i]=0;
    n[i]=results.at(i).n;
    method[i]=results.at(i).method;
    q_kurt[i]=results.at(i).q_kurt;
    q_skew[i]=results.at(i).q_skew;
    q_scale[i]=results.at(i).q_scale;
    q_location[i]=results.at(i).q_location;
    convergence[i]=results.at(i).convergence;
    iterations[i]=results.at(i).iterations;
    cpu_time[i]=results.at(i).cpu_time;
  }
  return DataFrame::create(Named("method", method),
                          Named("alpha",alpha),
                          Named("beta",beta),
                          Named("gamma",gamma),
                          Named("delta",delta),
                          Named("pm",pm),
                          Named("two_ll_n",two_ll_n),
                          Named("n", n),
                          Named("q_kurt",q_kurt),
                          Named("q_skew",q_skew),
                          Named("q_scale",q_scale),
                          Named("q_location",q_location),
                          Named("convergence",convergence),
                          Named("iterations",iterations),
                          Named("cpu_time",cpu_time));
  END_RCPP
}

// [[Rcpp::export]]
DataFrame g_map_dataframe(double x, double alpha, double beta) {
  BEGIN_RCPP
  // We need an integration controller and verbose indicator
  // to form a StandardStableDistribution but we won't actually use them in this routine
  int subdivisions=1000;
  int verbose = 0;
  bool noext = true;  // No extrapolation.  It's too unreliable.
  int N = 10;  // Us a 10 point Gausss and 21 point Gauss Kronrod rule
  double epsabs=0.;
  double epsrel=64*std::numeric_limits<double>::epsilon();
  IntegrationController ctl(noext, N, epsabs, epsrel, subdivisions, verbose);
  StandardStableDistribution std_stable_dist(alpha, beta, ctl, verbose);
  std_stable_dist.pdf(x,false);  // Just to force the mapping process
  int n = std_stable_dist.points.size();
  int m = std_stable_dist.last;
  DoubleVector th(n+m+1);
  DoubleVector g_th(n+m+1);
  CharacterVector type(n+m+1);
  for (int i=0; i<n; ++i) {
    th(i) = std_stable_dist.points.at(i);
    g_th(i) = std_stable_dist.g(std_stable_dist.points.at(i));
    type(i) = "g_map";
  }
  double b_max=-std::numeric_limits<double>::infinity();
  for (int i=0; i<m; ++i) {
    th(i+n) = std_stable_dist.controller->subs.at(i).a;
    g_th(i+n) = std_stable_dist.g(std_stable_dist.controller->subs.at(i).a);
    b_max = std::max(b_max,std_stable_dist.controller->subs.at(i).b);
    type(i+n) = "subinterval";
  }
  th(n+m) = b_max;
  g_th(n+m) = std_stable_dist.g(b_max);
  type(n+m) = "subinterval";
  return DataFrame::create(Named("theta", th),
                   Named("g",g_th),
                   Named("type",type));
  END_RCPP
}

