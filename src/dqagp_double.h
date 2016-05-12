//
//  dqagp.h
//
//
//  Created by Joseph Dunn on 11/19/15.
//
//

#ifndef dqagp_double_h
#define dqagp_double_h
#include <vector>
#include <cstddef>
#include <cmath>

using std::vector;
using std::ostream;

typedef void integr_fn_double(double *x, int n, void *ex);
/* vectorizing function   f(x[1:n], ...) -> x[]  {overwriting x[]}. */

class integr_ctl_double {
public:
    bool  noext;
    double wtg_level;
    int N;
    vector<double> w_gauss;
    vector<double> x_kronrod;
    vector<double> w_kronrod;
    double epsabs;
    double epsrel;
    int limit;
    int verbose;
    integr_ctl_double(bool noext, double wtg_level, int N, double epsabs, double epsrel, int limit, int verbose);
    friend ostream& operator<< (ostream& os, integr_ctl_double& ctl);
};


class subinterval_double {
public:
  double a;
  double b;
  double r;
  double e;
  double rabs;
  double defabs;
  int level;
  int ndin;
  static double epmach;
  static double uflow;
  static integr_ctl_double* ctl;
  static bool initialized;
  static void initialize(integr_ctl_double* ctl);
  subinterval_double() : a(0), b(0), r(0), e(0), ndin(0){};
  subinterval_double (double a, double b, double r, double e, int ndin):
    a(a), b(b), r(r), e(e), ndin(ndin) {};
  bool operator<(const subinterval_double rhs) const {
    if ((!is_divisible()) && rhs.is_divisible()) return true;
    else if (is_divisible() && !rhs.is_divisible()) return false;
    else if (ndin < rhs.ndin) return true;
    else if (ndin > rhs.ndin) return false;
    else if (log(e)-ctl->wtg_level*level
               < log(rhs.e) -ctl->wtg_level*rhs.level )return true;
    else return false;
  }
  bool is_divisible() const{
    return e!=0 && (fmax(fabs(a), fabs(b)) > (1 +
        100 * epmach) * (fabs((a+b)/2) + 1000 * uflow));
  }
  void dqk(integr_fn_double f, void* ex);
};

void print_subs(const vector<subinterval_double> subs, const int last, const vector<double> points);

void
dqagpe(
       integr_fn_double f,
       void* ex,
       const vector<double>& points,
       double& result,
       double& abserr,
       int& neval,
       int& ier,
       vector<subinterval_double>& subs,
       int& last);





#endif /* dqagp_h */
