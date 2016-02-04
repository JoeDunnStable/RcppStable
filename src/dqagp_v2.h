//
//  dqagp.h
//
//
//  Created by Joseph Dunn on 11/19/15.
//
//

#ifndef dqagp_h
#define dqagp_h
#include <array>
#define LIMIT 1000

typedef void integr_fn(double *x, int n, void *ex);
/* vectorizing function   f(x[1:n], ...) -> x[]  {overwriting x[]}. */

class subinterval {
public:
  double a;
  double b;
  double r;
  double e;
  double rabs;
  double defabs;
  int level;
  int ndin;
  static const double epmach;
  static const double uflow;
  static const double wtg_level;
  subinterval() : a(0), b(0), r(0), e(0), ndin(0) {};
  subinterval (double a, double b, double r, double e, int ndin ):
    a(a), b(b), r(r), e(e), ndin(ndin) {};
  bool operator<(const subinterval rhs) const {
    if ((!is_divisible()) && rhs.is_divisible()) return true;
    else if (is_divisible() && !rhs.is_divisible()) return false;
    else if (ndin < rhs.ndin) return true;
    else if (ndin > rhs.ndin) return false;
    else if (log(e)-wtg_level*level
               < log(rhs.e) -wtg_level*rhs.level )return true;
    else return false;
  }
  bool is_divisible() const{
    return (fmax(fabs(a), fabs(b)) > (0.1e+01 +
        0.1e+03 * epmach) * (fabs((a+b)/2) + 0.1e+04 * uflow));
  }
  void
    dqk21_v2(
      integr_fn f, void* ex);
};



void
dqagpe_v2(
       integr_fn f,
       void* ex,
       double const& a,
       double const& b,
       int const& npts2,
       const std::array<double, 100>& points,
       double const& epsabs,
       double const& epsrel,
       int const& limit,
       double& result,
       double& abserr,
       int& neval,
       int& ier,
       std::array<subinterval,LIMIT>& subs,
       int& last);



#endif /* dqagp_h */
