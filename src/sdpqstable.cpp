#include "stable_double.h"
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

vec sdstable(const vec x, const double alpha, const double beta, const int log_flag,
             integr_ctl_double& ctl, double zeta_tol, int verbose){
  vec ret(x.n_elem);
  uword i;
  g_double_class param(alpha, beta, ctl, zeta_tol, verbose);
  for (i=0; i<x.n_elem; i++)
    ret(i)=param.sdstable1(x(i), log_flag);
  return ret;
}

// [[Rcpp::export]]
NumericVector sdstable(NumericVector x, double alpha, double beta, int log_flag,
                       double tol, double zeta_tol, int subdivisions, int verbose)
{
  BEGIN_RCPP
  vec ax(x);
  bool noext = true;  // No extrapolation.  It's too unreliable.
  double wtg_level = 0.;  // We'll go after the subinterval with the largest error.
  int N = 10;  // Us a 10 point Gausss and 21 point Gauss Kronrod rule
  double epsabs=0.;
  double epsrel=tol;
  integr_ctl_double ctl(noext, wtg_level, N, epsabs, epsrel, subdivisions, verbose);
  return wrap(sdstable(ax,alpha,beta, log_flag, ctl, zeta_tol, verbose));
  END_RCPP
}

vec spstable(vec z, double alpha, double beta, int lower_tail, int log_p,
                       integr_ctl_double& ctl, double zeta_tol, int verbose) {
  vec ret(z.n_elem);
  uword i;
  g_double_class param(alpha, beta, ctl, zeta_tol, verbose);
  for (i=0; i<z.n_elem; i++)
    ret(i)=param.spstable1(z(i), lower_tail, log_p);
  return ret;
}

// [[Rcpp::export]]
NumericVector spstable(NumericVector z, double alpha, double beta, int lower_tail, int log_p,
                       double dbltol, int subdivisions, int verbose) {
  BEGIN_RCPP
  vec az(z);
  bool noext = true;  // No extrapolation.  It's too unreliable.
  double wtg_level = 0.;  // We'll go after the subinterval with the largest error.
  int N = 10;  // Us a 10 point Gausss and 21 point Gauss Kronrod rule
  double epsabs=0.;
  double epsrel=dbltol;
  double zeta_tol=0.;  // Not used by pstable
  integr_ctl_double ctl(noext, wtg_level, N, epsabs, epsrel, subdivisions, verbose);
  return wrap(spstable(az, alpha, beta, lower_tail, log_p, ctl, zeta_tol, verbose));
  END_RCPP
}

vec sqstable(const vec p, const double alpha, const double beta, const int lower_tail, int log_p,
                       const double dbltol, integr_ctl_double& ctl,
                       const double zeta_tol, const int verbose) {
  vec ret(p.n_elem);
  uword i;
  g_double_class param(alpha, beta, ctl, zeta_tol, verbose);
  for (i=0; i<p.n_elem; i++)
    ret(i)=param.sqstable1(p(i), lower_tail, log_p, dbltol);
  return ret;
}

// [[Rcpp::export]]
NumericVector sqstable(NumericVector p, double alpha, double beta, int lower_tail, int log_p,
                       double dbltol, double integ_tol, int subdivisions, int verbose) {
  BEGIN_RCPP
  vec ap(p);
  bool noext = true;  // No extrapolation.  It's too unreliable.
  double wtg_level = 0.;  // Go after the subinterval with the largest error first.
  int N = 10;  // Use a 10 point Gausss and 21 point Gauss Kronrod rule
  double epsabs=0.;
  double epsrel=integ_tol;
  double zeta_tol = 0.;
  integr_ctl_double ctl(noext, wtg_level, N, epsabs, epsrel, subdivisions, verbose);
  return wrap(sqstable(ap, alpha, beta, lower_tail, log_p,
                     dbltol, ctl, zeta_tol, verbose));
  END_RCPP
}

vec ddx_sdstable(const vec x, const double alpha, const double beta,
             integr_ctl_double& ctl, double zeta_tol, int verbose){
  vec ret(x.n_elem);
  uword i;
  g_double_class param(alpha, beta, ctl, zeta_tol, verbose);
  for (i=0; i<x.n_elem; i++)
    ret(i)=param.ddx_sdstable1(x(i));
  return ret;
}

// [[Rcpp::export]]
NumericVector ddx_sdstable(NumericVector x, double alpha, double beta,
                       double tol, double zeta_tol, int subdivisions, int verbose)
{
  BEGIN_RCPP
  vec ax(x);
  bool noext = true;  // No extrapolation.  It's too unreliable.
  double wtg_level = 0.;  // We'll go after the subinterval with the largest error.
  int N = 10;  // Us a 10 point Gausss and 21 point Gauss Kronrod rule
  double epsabs=0.;
  double epsrel=tol;
  integr_ctl_double ctl(noext, wtg_level, N, epsabs, epsrel, subdivisions, verbose);
  return wrap(ddx_sdstable(ax, alpha, beta, ctl, zeta_tol, verbose));
  END_RCPP
}

vec srstable(const double alpha, const double beta, vec u1, vec u2){
  vec ret(u1.n_elem);
  uword i;
  for (i=0; i<u1.n_elem; i++)
    ret(i)=srstable1(alpha, beta, u1[i], u2[i]);
  return ret;
}

// [[Rcpp::export]]
NumericVector srstable(double alpha, double beta, NumericVector u1, NumericVector u2)
{
  BEGIN_RCPP
  vec au1(u1), au2(u2);
  return wrap(srstable(alpha, beta, au1, au2));
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

vec dt(vec x, double df, int log_p=0){
  vec ret(x.n_elem);
  uword i;
  for (i=0; i<x.n_elem; i++)
    ret(i)=Rf_dt(x(i),df,log_p);
  return ret;
}

vec pt(vec x, double df, int lower_tail=1, int log_p=0){
  vec ret(x.n_elem);
  uword i;
  for (i=0; i<x.n_elem; i++)
    ret(i)=Rf_pt(x(i),df,lower_tail,log_p);
  return ret;
}

vec qt(vec p, double df, int lower_tail=1, int log_p=0){
  vec ret(p.n_elem);
  uword i;
  for (i=0; i<p.n_elem; i++)
    ret(i)=Rf_qt(p(i),df,lower_tail,log_p);
  return ret;
}

#include "cubicspline.h"

// [[Rcpp::export]]
NumericVector cubic_spline(NumericVector x, NumericVector x_knots, NumericVector y_knots,
                           bool isclamped, NumericVector deriv) {
  cubicspline tst(x_knots,y_knots,isclamped,deriv);
  tst.get_knots();
  tst.get_coefs();
  return wrap(tst(x));
}

// [[Rcpp::export]]
NumericVector sdstable_quick(NumericVector x, double alpha, double beta,
                       double tol, double zeta_tol, int subdivisions, int verbose)
{
  BEGIN_RCPP
  vec ax(x);
  bool noext = true;  // No extrapolation.  It's too unreliable.
  double wtg_level = 0.;  // We'll go after the subinterval with the largest error.
  int N = 10;  // Us a 10 point Gausss and 21 point Gauss Kronrod rule
  double epsabs=0.;
  double epsrel=tol;
  integr_ctl_double ctl(noext, wtg_level, N, epsabs, epsrel, subdivisions, verbose);
  uword n = ax.n_elem;
  vec ret(n);
  if (n>200 && alpha !=2){
    ret.fill(R_NegInf);
    // the splines don't work well near the mode so we'll use dstable for everything near the mode
    double x_mode=stable_mode(alpha,beta,1-1e-11,ctl, zeta_tol,verbose);
    uvec ordx = sort_index(ax);
    uvec sel_exact= (ax>x_mode-.1) % (ax<x_mode+.1);
    sel_exact = sel_exact + (ordx<10) + (ordx>=n-10);
    uvec sel_zero(n,fill::zeros);
    if (alpha<1) {
      if (beta==1)
        sel_zero=ax<=-tan(M_PI_2*alpha);
      else if (beta==-1)
        sel_zero=ax>=tan(M_PI_2*alpha);
    }
    uvec sel_high = (sel_exact == 0) % (ax > x_mode) % (sel_zero == 0);
    uvec sel_low = (sel_exact == 0) % (ax < x_mode) % (sel_zero == 0);
    if (sum(sel_high)<85) {
        sel_exact = sel_exact + sel_high;
    } else {
        double x_high_break1=x_mode+1;
        double x_high_break2=x_mode+10;
        vec x_high_inner = join_cols(linspace(x_mode+.1,x_high_break1,30),
                        linspace(x_high_break1,x_high_break2,50).tail_rows(49));
        double pt_high_break2=Rf_pt(x_high_break2,alpha,1,0);
        vec pt_high_outer = linspace(pt_high_break2,1-.000001,
                                     fmax(10,500*(1-pt_high_break2))).tail_rows(9);
        vec pt_high_inner(pt(x_high_inner,alpha));
        vec pt_knot_high=join_cols(pt_high_inner,pt_high_outer);
        vec x_high_outer(qt(pt_high_outer,alpha));
        vec x_knot_high = join_cols(x_high_inner,x_high_outer);
        vec y_knot_high = sdstable(x_knot_high,alpha,beta,TRUE,ctl ,zeta_tol,verbose)
                                - dt(x_knot_high,alpha,TRUE);
        uvec sel_high_finite = find_finite(y_knot_high);
        pt_knot_high = pt_knot_high(sel_high_finite);
        x_knot_high = x_knot_high(sel_high_finite);
        y_knot_high = y_knot_high(sel_high_finite);
        cubicspline spline_high(pt_knot_high,y_knot_high,FALSE,zeros(2));
        vec pt_high(pt(ax(find(sel_high)),alpha));
        ret(find(sel_high)) = spline_high(pt_high)+dt(ax(find(sel_high)),alpha,TRUE);
      }
      if (sum(sel_low)<85) {
        sel_exact = sel_exact + sel_low;
      } else {
        double x_low_break2 = x_mode-10;
        double x_low_break1 = x_mode-1;
        vec x_low_inner=join_cols(linspace(x_low_break2,x_low_break1,50),
                      linspace(x_low_break1, x_mode-.1,30).tail_rows(29));
        double pt_low_break_2=Rf_pt(x_low_break2,alpha,1,0);
        vec pt_low_outer=linspace(0.000001,pt_low_break_2,
                                  fmax(10,500*pt_low_break_2)).head_rows(9);
        vec pt_low_inner(pt(x_low_inner,alpha));
        vec pt_knot_low = join_cols(pt_low_outer,pt_low_inner);
        vec x_low_outer(qt(pt_low_outer,alpha));
        vec x_knot_low = join_cols(x_low_outer,x_low_inner);
        vec y_knot_low = sdstable(x_knot_low, alpha, beta, TRUE, ctl, zeta_tol, verbose)
                         - dt(x_knot_low, alpha, TRUE);
        uvec sel_low_finite = find_finite(y_knot_low);
        pt_knot_low = pt_knot_low(sel_low_finite);
        x_knot_low = x_knot_low(sel_low_finite);
        y_knot_low = y_knot_low(sel_low_finite);
        cubicspline spline_low(pt_knot_low,y_knot_low,FALSE,zeros(2));
        vec pt_low(pt(ax(find(sel_low)),alpha));
        ret(find(sel_low))=spline_low(pt_low)+dt(ax(find(sel_low)),alpha,TRUE);
      }

      ret(find(sel_exact))=sdstable(ax(find(sel_exact)),alpha,beta,TRUE, ctl, zeta_tol, verbose);
  } else {
    ret=sdstable(ax, alpha, beta, TRUE, ctl, zeta_tol, verbose);
  }
  return wrap(ret);
  END_RCPP;
}

