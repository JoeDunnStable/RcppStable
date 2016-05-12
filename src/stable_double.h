#ifndef STABLE_double_H
#define STABLE_double_H

#include "dqagp_double.h"
#include <iostream>
#include <vector>

#include <RcppArmadillo.h>
#define cout Rcpp::Rcout
#define cerr Rcpp::Rcerr
// using std::cout;
// using std::cerr;
using std::endl;
using std::ostream;
using std::vector;
using std::string;

class g_double_class {
public:
  static double Machine_eps;
  static const double alpha_small_dstable;
  static double large_exp_arg;
  static double pi;
  static double pi2;
  static double PosInf;
  static double NegInf;
  double alpha;
  double beta_input;      //As input.  Sign might be flipped in the computation
  double x_input;
  // variables independent of x used when alpha != 1
  double zeta;
  double theta0_x_gt_zeta;
  double cat0;
  // variable dependent on x used when alpha != 1
  double beta;
  double theta0;
  double x_m_zet;
  double add_l;
  double add_r;
  double th_max;         //Both th_l and th_r range from 0 to th_max
  double c2;
  double c_ddx;
  enum dist {Cauchy=1, normal=2, fin_support=3, other=4};
  dist dist_type;
  enum fun {fun_g_l=1, fun_g_r=2, fun_ga1_r=3};
  fun fun_type;

  // variable used when alpha = 1
  double abs_x;
  double i2b;
  double p2b;
  double ea;
  double u0;

  // variables controling the integrator
  integr_ctl_double* ctl_final;
  double zeta_tol;
  int verbose;

  // the results of the map_g process. Ussed to establish the initial subintervals.
  bool good_theta2;
  double theta2, g_theta2, g_dd_theta2;
  vector<double> points;

  // work area for integrator
  static vector<subinterval_double> subs;

  // additional output from integrator
  double abserr;
  int neval;
  int ier;

  g_double_class(double alpha, double beta, integr_ctl_double& ctl, double zeta_tol, int verbose) : alpha(alpha), beta_input(beta),
      ctl_final(&ctl), zeta_tol(zeta_tol), verbose(verbose), x_input(NAN) {
      subs.resize(ctl_final->limit);
    if (alpha!=1){
      zeta = -beta_input*tan(alpha*pi2);
      theta0_x_gt_zeta = atan(beta_input*tan(alpha*pi2))/alpha;
      theta0_x_gt_zeta = fmin(pi2,fmax(-pi2,theta0_x_gt_zeta));
      cat0=cos(alpha*theta0_x_gt_zeta);
    } else {
      c2=pi2*fabs(1/(2*beta_input));
      c_ddx=-c2*pi2/beta_input;
    }

  };

  void set_x(double x){
      if (x!=x_input) {
          x_input = x;
          if (!isfinite(x)) {
              if (verbose)
                  cout << "set_x: x is not finite." << endl;
              return;
          }
          if (alpha == 1 && beta_input == 0) {
              if (verbose) {
                  cout << "set_x: Cauchy distribution" <<endl;
              }
              dist_type=Cauchy;
              return;
          }
          if (alpha == 2) {
              if (verbose) {
                  cout << "set_x: normal distribution with std dev of sqrt(2)" <<endl;
              }
              dist_type=normal;
              return;
          }
          if (alpha < 1 && ((beta_input==1 && x<=zeta) || (beta_input==-1 && x>=zeta))) {
              if (verbose) {
                  cout << "set_x: outside of distribution support" << endl;
              }
              dist_type=fin_support;
              return;
          }
          dist_type=other;
          if (alpha!=1){
              x_m_zet=fabs(x-zeta);
              if (x>=zeta){
                beta=beta_input;
                theta0 =theta0_x_gt_zeta;
              } else {
                beta=-beta_input;
                theta0 = -theta0_x_gt_zeta;
              }
              th_max=(pi2)+theta0;
              add_l=theta0-(pi2);
              add_r=fmax(0.,pi-alpha*(theta0+(pi2)));
              if (alpha<1) {
                  if (beta==1)
                      add_l=0;
                  else if (beta==-1)
                      add_l=pi;
              } else if (alpha>1) {
                  if (beta==-1)
                      add_r=0;
              }
              c2 = (alpha/(pi * fabs(alpha-1) * x_m_zet));
              c_ddx = c2/((alpha-1)*(x-zeta));
              if (x_m_zet<fmax(1,fabs(zeta)))
                  fun_type=fun_g_l;
              else
                  fun_type=fun_g_r;
          } else { // alpha = 1
              abs_x=fabs(x);
              if (x >= 0) {
                  beta = beta_input;
              } else {
                  beta = -beta_input;
              }
              i2b=1/(2*beta);
              p2b=pi*i2b;
              u0=-((beta>=0)-(beta<=0));
              ea = -p2b*abs_x;
              th_max=2;
              fun_type=fun_ga1_r;
          }
          map_g();
      }
  }

  double g(double th) const{
    switch(fun_type){
    case fun_g_l:
      return g_l(th);
    case fun_g_r:
      return g_r(th);
    case fun_ga1_r:
      return ga1_r(th);
    }
  }

  double g_l(double th_l) const;   //th_l is th+theta0
  double g_r(double th_r) const;   //th_r is pi/2-th
  double ga1_r(double u_r) const;  //u_r is 1-u

  double dlng_dth(double th) {
    switch(fun_type){
      case fun_g_l:
        return dlng_dth_l(th);
      case fun_g_r:
        return dlng_dth_r(th);
      case fun_ga1_r:
        return dlnga1_du_r(th);
    }
  }

  double dlng_dth_l(double th_l);
  double dlng_dth_r(double th_r);
  double dlnga1_du_r(double u_r);

   void map_g();
  void th_guess(const double &value, double &lower, double &g_lower, double &upper, double &g_upper);
  double integrate_dstable(bool log_flag);
  double integrate_ddx_dstable();
  double integrate_pstable();
  double sdstable1(double x, int log_flag);
  double spstable1(double z, int lower_tail, int log_p);
  double sqstable1(double p, int lower_tail, int log_p, double dbltol);
  double ddx_sdstable1(double x);

  friend ostream& operator<<(ostream& os, const g_double_class& gc);

}; // g_double_class

double srstable1(double alpha, double beta, double u1, double u2);

// This class contains everything needed to integrate f(g) except the endpoints
// and starting subintervals of the integration range.
class Int_f_of_g_double {
private:
    g_double_class *param;
    static string msgs[];
    double (*f)(double, g_double_class*);

public:
    double f_of_g(double th) {return (*f)(param->g(th), param);}
    Int_f_of_g_double(double (*f)(double, g_double_class*), g_double_class* param, integr_ctl_double* ctl) : f(f),
    param(param) {
        subinterval_double::initialize(ctl);
    }
    string msg() {return msgs[ier];};
    g_double_class* get_param() {return param;}
    double result;
    double abserr;
    int neval;
    int ier;
    int last;
    double operator() ();
};

// Functor passed to toms768 solve to find x such that g(x) == value
class g_double_solve {
private:
  g_double_class* param;
  double value;
  bool log_flag;
  double g_min;
  double g_max;

public:
  g_double_solve(double value, g_double_class* param, bool log_flag=false) :
    param(param), value(value), log_flag(log_flag), g_min(std::numeric_limits<double>::min()),
    g_max(std::numeric_limits<double>::max()) {}
  double operator()(const double th) {
    double g_ = param->g(th);
    g_=fmax(g_min,fmin(g_,g_max));
    if (param->verbose >=3)
      cout << "theta = " << th
            << ", g(theta) = " << g_ << endl;
    if (log_flag)
      return log(g_)-value;
    else
      return g_ - value;
  }
  void set_value(double value_in) {value=value_in;}
}; //g_double_solve

double dPareto(double x, double alpha, double beta, bool log_flag);
double pPareto(double x, double alpha, double beta, bool lower_tail, bool log_p);

double stable_mode(double alpha, double beta, double dbltol, integr_ctl_double& ctl, double zeta_tol, int verbose);

// Functor passed to toms748_solve to find x such that ddx_dstable(x) == value
class ddx_sdstable_double_solve {
private:
  double value;
  g_double_class param;
  double tol;
public:
  ddx_sdstable_double_solve(double value, double alpha, double beta, double tol, integr_ctl_double& ctl, double zeta_tol, int verbose):
  value(value),
  param(alpha, beta, ctl, zeta_tol, verbose),
  tol(tol){}
  double operator()(const double x) {
    double ret=param.ddx_sdstable1(x);
      if (param.verbose >=3)
      cout << "x = " << x
            << ", ddx_dstable(x) = " << ret << endl;
      return ret - value;
  }
};

// These classes are used in several place when running toms748_solve
class rel_eps_tolerance_double
{
public:
  rel_eps_tolerance_double(double eps) : eps(eps) {};
  inline bool operator()(const double& a, const double& b)
  {
    return fabs(a - b) <= (eps * (std::min)(fabs(a), fabs(b)));
  }
private:
  double eps;
};

class abs_eps_tolerance_double
{
public:
  abs_eps_tolerance_double(double eps_in) : eps(eps_in){};
  inline bool operator()(const double& a, const double& b)
  {
    bool r=fabs(a - b) <= eps;
    return r;
  }
private:
  double eps;
};



#endif

