#ifndef STABLE_H
#define STABLE_H

#include <RcppArmadillo.h>
#include <cmath> //For trig functions
#include <iostream>
#include <string>

using Rcpp::Rcout;
using std::endl;
using std::string;

//#include <Rmath.h>

//using namespace Rcpp;

class g_class {
public:
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
  enum fun {fun_g_l=1, fun_g_r=2, fun_ga1_r=3};
  fun fun_type;

  // variable used when alpha = 1
  double abs_x;
  double i2b;
  double p2b;
  double ea;
  double u0;

  g_class(double alpha, double beta) : alpha(alpha), beta_input(beta) {
    if (alpha!=1){
      zeta = -beta_input*tan(alpha*M_PI_2);
      theta0_x_gt_zeta = atan(beta_input*tan(alpha*M_PI_2))/alpha;
      theta0_x_gt_zeta = fmin(M_PI_2,fmax(-M_PI_2,theta0_x_gt_zeta));
      cat0=cos(alpha*theta0_x_gt_zeta);
    } else {
      c2=M_PI_2*fabs(1/(2*beta_input));
      c_ddx=-c2*M_PI_2/beta_input;
    }

  };

  void set_x(double x){
    x_input=x;
    if (alpha!=1){
      x_m_zet=fabs(x-zeta);
      if (x>=zeta){
        beta=beta_input;
        theta0 =theta0_x_gt_zeta;
      } else {
        beta=-beta_input;
        theta0 = -theta0_x_gt_zeta;
      }
      th_max=M_PI_2+theta0;
      add_l=theta0-M_PI_2;
      add_r=fmax(0.,M_PI-alpha*(theta0+M_PI_2));
      if (alpha<1) {
        if (beta==1)
          add_l=0;
        else if (beta==-1)
          add_l=M_PI;
      } else if (alpha>1) {
        if (beta==-1)
          add_r=0;
      }
      c2 = (alpha/(M_PI * fabs(alpha-1) * x_m_zet));
      c_ddx = c2/((alpha-1)*(x_input-zeta));
      if (x_m_zet<1e-2*fmax(1,fabs(zeta)))
        fun_type=fun_g_l;
      else
        fun_type=fun_g_r;
    } else { // alpha = 1
      abs_x=fabs(x_input);
      if (x >= 0) {
        beta = beta_input;
      } else {
        beta = -beta_input;
      }
      i2b=1/(2*beta);
      p2b=M_PI*i2b;
      u0=-((beta>=0)-(beta<=0));
      ea = -p2b*abs_x;
      th_max=2;
      fun_type=fun_ga1_r;
    }
  }

  double g(double th) {
    switch(fun_type){
    case fun_g_l:
      return g_l(th);
    case fun_g_r:
      return g_r(th);
    case fun_ga1_r:
      return ga1_r(th);
    }
  }

  double g_l(double th_l);   //th_l is th+theta0
  double g_r(double th_r);   //th_r is pi/2-th
  double ga1_r(double u_r);  //u_r is 1-u

  double integrate_g_exp_m_g(bool log_flag, double dbltol, int subdivisions,
                                       double zeta_tol, int verbose);
  double integrate_g_1_m_alpha_g_exp_m_g(double dbltol, int subdivisions,
                             double zeta_tol, int verbose);
  double integrate_exp_m_g(bool giveI, double dbltol, int subdivisions, int verbose = 0);
  double sdstable1(double x, int log_flag,
                   double tol, double zeta_tol, int subdivisions, int verbose);
  double spstable1(double z, int lower_tail, int log_p,
                   double dbltol, int subdivisions, int verbose);
  double sqstable1(double p, int lower_tail, int log_p,
                            double dbltol, double integ_tol, int subdivisions, int verbose);
  double ddx_sdstable1(double x, double tol, double zeta_tol,
                       int subdivisions, int verbose);

  friend std::ostream& operator<<(std::ostream& os, const g_class& gc);

}; // g_class

// Functor passed to toms768 solve to find x such that g(x) == value
class g_solve {
private:
  g_class* param;
  double value;
  bool log_flag;
  int verbose;

public:
  g_solve(double value, g_class* param, int verbose, bool log_flag=false) :
    param(param), value(value), log_flag(log_flag), verbose(verbose) {}
  double operator()(const double th) {
    double g_ = param->g(th);
    g_=fmax(1e-300,fmin(g_,1e100));
    if (verbose >=3)
      Rcout << "theta = " << th
            << ", g(theta) = " << g_ << std::endl;
    if (log_flag)
      return log(g_)-value;
    else
      return g_ - value;
  }
  void set_value(double value_in) {value=value_in;}
}; //g_solve

void g_exp_m_g(double* th, int n, void* ext);
void exp_m_g(double* th, int n, void* ext);
void one_m_exp_m_g(double* th, int n, void* ext);
void g_1_m_alpha_g_exp_m_g(double* th, int n, void* ext);

// Functor passed to toms748_solve to find x such that g1(x) == value
class g_exp_m_g_solve {
private:
  double value;
  g_class *param;
  int verbose;
public:
  g_exp_m_g_solve(double value,g_class *param, int verbose) :
    value(value),
    param(param),
    verbose(verbose){}
  double operator()(const double th) {
    double g1_=th;
    g_exp_m_g(&g1_,1,param);
    g1_=fmin(g1_,1e100);
    if (verbose >=3)
      Rcout << "theta = " << th
            << ", g1(theta) = " << g1_ << std::endl;
    return g1_ - value;
  }
};

// Functor passed to toms748_solve to find x such that g1(x) == value
class exp_m_g_solve {
private:
  double value;
  g_class *param;
  int verbose;
public:
  exp_m_g_solve(double value,g_class *param, int verbose) :
  value(value),
  param(param),
  verbose(verbose){}
  double operator()(const double th) {
    double g1_=th;
    exp_m_g(&g1_,1,param);
    if (verbose >=3)
      Rcout << "theta = " << th
            << ", exp(-g(theta)) = " << g1_ << std::endl;
      return g1_ - value;
  }
};


// Functor passed to toms748_solve to find x such that ddx_dstable(x) == value
class ddx_sdstable_solve {
private:
  double value;
  g_class param;
  double tol;
  double zet_tol;
  int subdivisions;
  int verbose;
public:
  ddx_sdstable_solve(double value, double alpha, double beta,
                              double tol, double zet_tol, int subdivisions,
                              int verbose) :
  value(value),
  param(alpha, beta),
  tol(tol),
  zet_tol(zet_tol),
  subdivisions(subdivisions),
  verbose(verbose){}
  double operator()(const double x) {
    double ret=param.ddx_sdstable1(x, tol, zet_tol, subdivisions, verbose);
    if (verbose >=3)
      Rcout << "x = " << x
            << ", g_1_m_alpha_g_exp_m_g(x) = " << ret << std::endl;
      return ret - value;
  }
};

// These classes are used in several place when running toms748_solve
class rel_eps_tolerance
{
public:
  rel_eps_tolerance(double eps) : eps(eps) {};
  inline bool operator()(const double& a, const double& b)
  {
    return fabs(a - b) <= (eps * (std::min)(fabs(a), fabs(b)));
  }
private:
  double eps;
};

class abs_eps_tolerance
{
public:
  abs_eps_tolerance(double eps_in) : eps(eps_in){};
  inline bool operator()(const double& a, const double& b)
  {
    bool r=fabs(a - b) <= eps;
    return r;
  }
private:
  double eps;
};

void warning(string msg);
void stop(string msg);

#endif

