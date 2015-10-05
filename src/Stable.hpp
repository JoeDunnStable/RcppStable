#ifndef STABLE_HPP
#define STABLE_HPP

#include <cmath> //For trig functions
#include <Rmath.h>
using namespace Rcpp;

class g_param {
public:
  double alpha;
  double at0;
  double cat0;
  double x_m_zet;
};

double g(double th, g_param* param);

/* These are the ancillary functions used by the R function dstable */
double fct1(double x, double zeta,
                   double alpha, double beta, double theta0,
                   bool log_flag, double tol, long subdivisions,
                   double zeta_tol, int verbose = 0);

double fct2(double x, double beta, bool log_flag,
                   double dbltol, int subdivisions, int verbose);

double sdstable1(double x, double alpha, double beta, int log_flag,
                       double tol, double zeta_tol, int subdivisions, int verbose);

//' @title tan(pi/2*x), for x in [-1,1] with correct limits
//'   i.e. tanpi2(-/+ 1) == -/+ Inf
//' @param x numeric vector
//' @return numeric vector of values tan(pi/2*x)
//' @author Martin Maechler
inline double tanpi2(double x) {
  double r;
  double vec[] = {0,R_PosInf,0,R_NegInf};
  if(x == round(x))
    r = vec[(int)x%4];
  else
    r = tan(M_PI_2* x);
  return r;
}

//' @title cos(pi/2*x), for x in [-1,1] with correct limits
//'   i.e. cospi2(+- 1) == 0
//' @param x numeric vector
//' @return numeric vector of values cos(pi/2*x)
//' @author Martin Maechler
inline double cospi2(double x) {
  double vec[]={1,0,-1,0};
  double r;
  if (x == round(x))
    r = vec[(int)x]; // 1 or 0 - iff x \in [-1,1] !
  else
    r = cos(M_PI_2* x);
  return r;
}

/* These are the ancillary function used by the R function pstable */
double spstable1(double z, double alpha, double beta, int lower_tail, int log_p,
                 double dbltol, int subdivisions, int verbose);

double FCT1(double x, double zeta, double alpha, double theta0,
                   bool giveI, double dbltol, int subdivisions, int verbose = 0);

double FCT2(double x, double beta, double dbltol, int subdivisions,
                   bool giveI = 0, int verbose = 0);

double sqstable1(double p, double alpha, double beta, int lower_tail, int log_p,
                 double dbltol, double integ_tol, int subdivisions, int verbose);

// These classes are used in several place when running toms748_solve
class rel_eps_tolerance
{
public:
  rel_eps_tolerance(double eps_in) : eps(eps_in) {};
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

#endif

