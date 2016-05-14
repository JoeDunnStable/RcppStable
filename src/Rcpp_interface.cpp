#include "stable_double.h"
#include "sdpqstable_double.h"

using Rcpp::NumericVector;
using Rcpp::IntegerVector;
using Rcpp::DoubleVector;
using Rcpp::CharacterVector;
using Rcpp::DataFrame;
using Rcpp::Named;
using Rcpp::wrap;

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
NumericVector sdstable(Eigen::VectorXd x, double alpha, double beta, int log_flag,
                       double tol, double zeta_tol, int subdivisions, int verbose)
{
  BEGIN_RCPP
  bool noext = true;  // No extrapolation.  It's too unreliable.
  double wtg_level = 0.;  // We'll go after the subinterval with the largest error.
  int N = 10;  // Us a 10 point Gausss and 21 point Gauss Kronrod rule
  double epsabs=0.;
  double epsrel=tol;
  integr_ctl_double ctl(noext, wtg_level, N, epsabs, epsrel, subdivisions, verbose);
  return wrap(sdstable(x,alpha,beta, log_flag, ctl, zeta_tol, verbose));
  END_RCPP
}

// [[Rcpp::export]]
NumericVector sdstable_quick(Eigen::VectorXd x, double alpha, double beta, int log_flag,
                       double tol, double zeta_tol, int subdivisions, int verbose)
{
  BEGIN_RCPP
  bool noext = true;  // No extrapolation.  It's too unreliable.
  double wtg_level = 0.;  // We'll go after the subinterval with the largest error.
  int N = 10;  // Us a 10 point Gausss and 21 point Gauss Kronrod rule
  double epsabs=0.;
  double epsrel=tol;
  integr_ctl_double ctl(noext, wtg_level, N, epsabs, epsrel, subdivisions, verbose);
  return wrap(sdstable_quick(x, alpha, beta, log_flag, ctl, zeta_tol, verbose));
  END_RCPP
}

// [[Rcpp::export]]
NumericVector spstable(Eigen::VectorXd z, double alpha, double beta, int lower_tail, int log_p,
                       double dbltol, int subdivisions, int verbose) {
  BEGIN_RCPP
  bool noext = true;  // No extrapolation.  It's too unreliable.
  double wtg_level = 0.;  // We'll go after the subinterval with the largest error.
  int N = 10;  // Us a 10 point Gausss and 21 point Gauss Kronrod rule
  double epsabs=0.;
  double epsrel=dbltol;
  double zeta_tol=0.;  // Not used by pstable
  integr_ctl_double ctl(noext, wtg_level, N, epsabs, epsrel, subdivisions, verbose);
  return wrap(spstable(z, alpha, beta, lower_tail, log_p, ctl, zeta_tol, verbose));
  END_RCPP
}

// [[Rcpp::export]]
NumericVector sqstable(Eigen::VectorXd p, double alpha, double beta, int lower_tail, int log_p,
                       double dbltol, double integ_tol, int subdivisions, int verbose) {
  BEGIN_RCPP
  bool noext = true;  // No extrapolation.  It's too unreliable.
  double wtg_level = 0.;  // Go after the subinterval with the largest error first.
  int N = 10;  // Use a 10 point Gausss and 21 point Gauss Kronrod rule
  double epsabs=0.;
  double epsrel=integ_tol;
  double zeta_tol = 0.;
  integr_ctl_double ctl(noext, wtg_level, N, epsabs, epsrel, subdivisions, verbose);
  return wrap(sqstable(p, alpha, beta, lower_tail, log_p,
                     dbltol, ctl, zeta_tol, verbose));
  END_RCPP
}

// [[Rcpp::export]]
NumericVector ddx_sdstable(Eigen::VectorXd x, double alpha, double beta,
                       double tol, double zeta_tol, int subdivisions, int verbose)
{
  BEGIN_RCPP
  bool noext = true;  // No extrapolation.  It's too unreliable.
  double wtg_level = 0.;  // We'll go after the subinterval with the largest error.
  int N = 10;  // Us a 10 point Gausss and 21 point Gauss Kronrod rule
  double epsabs=0.;
  double epsrel=tol;
  integr_ctl_double ctl(noext, wtg_level, N, epsabs, epsrel, subdivisions, verbose);
  return wrap(ddx_sdstable(x, alpha, beta, ctl, zeta_tol, verbose));
  END_RCPP
}

// [[Rcpp::export]]
NumericVector srstable(double alpha, double beta, Eigen::VectorXd u1, Eigen::VectorXd u2)
{
  BEGIN_RCPP
  if (u1.size() != u2.size()) throw std::range_error("srstable: u1 & u2 sizes differ");
  size_t n = u1.size();
  Vec ret(n);
  for (size_t i=0; i<n; ++i) ret(i)=srstable1(alpha, beta, u1(i), u2(i));
  return wrap(ret);
  END_RCPP
}

// [[Rcpp::export]]
double sdstableMode(double alpha, double beta,
                    double dbltol,
                    double tol, double zeta_tol, int subdivisions, int verbose)
{
  bool noext = true;  // No extrapolation.  It's too unreliable.
  double wtg_level = 0.;  // We'll go after the subinterval with the largest error.
  int N = 10;  // Us a 10 point Gausss and 21 point Gauss Kronrod rule
  double epsabs=0.;
  double epsrel=tol;
  integr_ctl_double ctl(noext, wtg_level, N, epsabs, epsrel, subdivisions, verbose);
  return stable_mode(alpha, beta, dbltol, ctl, zeta_tol, verbose);
}

// [[Rcpp::export]]
DataFrame stable_fit_cpp(Eigen::VectorXd y, std::string type, bool quick) {
  BEGIN_RCPP
  bool noext = true;  // No extrapolation.  It's too unreliable.
  double wtg_level = 0.;  // We'll go after the subinterval with the largest error.
  int N = 10;  // Us a 10 point Gausss and 21 point Gauss Kronrod rule
  double epsabs=0.;
  double epsrel=64*std::numeric_limits<double>::epsilon();
  int subdivisions = 1000;
  int verbose = 0;
  integr_ctl_double ctl(noext, wtg_level, N, epsabs, epsrel, subdivisions, verbose);
  double dbltol=1e-12;

  std::vector<fit_result> results = stable_fit(y, ctl, dbltol, type, quick, verbose);
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


