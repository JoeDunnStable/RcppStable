#include "Stable.h"
#include <iostream>
#include <boost/math/tools/toms748_solve.hpp>
#include <boost/math/distributions/students_t.hpp>

using boost::math::tools::toms748_solve;
using boost::math::students_t;
using boost::math::cdf;
using boost::math::quantile;
static const double Machine_eps = std::numeric_limits<double>::epsilon();
static const double PosInf = std::numeric_limits<double>::infinity();
static const double NegInf = -PosInf;

// Functor passed to toms748 solve to find q for unit stable distribution
class p_solve {
private:
  double p;
  g_class* param;
  int lower_tail;
  int log_p;
  int maxiter;
  int verbose;
  double integ_tol;
  int subdivisions;

public:
  p_solve(double p, g_class* param,
          int lower_tail, int log_p,
          double integ_tol, int subdivisions, int verbose) :
  p(p), param(param), lower_tail(lower_tail) , log_p(log_p),
  verbose(verbose), integ_tol(integ_tol), subdivisions(subdivisions){}
  double operator()(const double q) {
    if (verbose)
      Rcout << "Calling spstable with parmerters" << std::endl
           << "q = " << q << ", alpha = " << param->alpha
           << ", beta = " << param->beta_input << std::endl
           << "targeting p = " << p;

      double r = param->spstable1(q, lower_tail, log_p,
                                  integ_tol, subdivisions, false)-p;
      if (verbose)
        Rcout << ", Resulting delta = " << r << std::endl;
      return r;
  }
};

namespace boost { namespace math { namespace tools {

template <class F, class T, class Tol, class Policy>
std::pair<T, T> bracket_and_solve_root2(F f, const T& guess, T factor, bool rising, Tol tol, boost::uintmax_t& max_iter, const Policy& pol)
{
  BOOST_MATH_STD_USING
  static const char* function = "bracket_and_solve_root2<%1%>";
  //
  // Set up inital brackets:
  //
  T a = guess;
  T b = a;
  T fa = f(a);
  T fb = fa;
  //
  // Set up invocation count:
  //
  boost::uintmax_t count = max_iter - 1;

  int step = 32;

  if((fa < 0) == rising)
  {
    //
    // Zero is to the right of b, so walk upwards
    // until we find it:
    //
    while((boost::math::sign)(fb) == (boost::math::sign)(fa))
    {
      if(count == 0)
        return boost::math::detail::pair_from_single(policies::raise_evaluation_error(function, "Unable to bracket root, last nearest value was %1%", b, pol));
      //
      // Heuristic: normally it's best not to increase the step sizes as we'll just end up
      // with a really wide range to search for the root.  However, if the initial guess was *really*
      // bad then we need to speed up the search otherwise we'll take forever if we're orders of
      // magnitude out.  This happens most often if the guess is a small value (say 1) and the result
      // we're looking for is close to std::numeric_limits<T>::min().
      //
      if((max_iter - count) % step == 0)
      {
        factor *= 2;
        if(step > 1) step /= 2;
      }
      //
      // Now go ahead and move our guess by "factor":
      //
      a = b;
      fa = fb;
      b += factor;
      fb = f(b);
      --count;
      BOOST_MATH_INSTRUMENT_CODE("a = " << a << " b = " << b << " fa = " << fa << " fb = " << fb << " count = " << count);
    }
  }
  else
  {
    //
    // Zero is to the left of a, so walk downwards
    // until we find it:
    //
    while((boost::math::sign)(fb) == (boost::math::sign)(fa))
    {
      if(count == 0)
        return boost::math::detail::pair_from_single(policies::raise_evaluation_error(function, "Unable to bracket root, last nearest value was %1%", a, pol));
      //
      // Heuristic: normally it's best not to increase the step sizes as we'll just end up
      // with a really wide range to search for the root.  However, if the initial guess was *really*
      // bad then we need to speed up the search otherwise we'll take forever if we're orders of
      // magnitude out.  This happens most often if the guess is a small value (say 1) and the result
      // we're looking for is close to std::numeric_limits<T>::min().
      //
      if((max_iter - count) % step == 0)
      {
        factor *= 2;
        if(step > 1) step /= 2;
      }
      //
      // Now go ahead and move are guess by "factor":
      //
      b = a;
      fb = fa;
      a -= factor;
      fa = f(a);
      --count;
      BOOST_MATH_INSTRUMENT_CODE("a = " << a << " b = " << b << " fa = " << fa << " fb = " << fb << " count = " << count);
    }
  }
  max_iter -= count;
  max_iter += 1;
  std::pair<T, T> r = toms748_solve(
    f, a, b, fa, fb,
    tol, count, pol);
  max_iter += count;
  BOOST_MATH_INSTRUMENT_CODE("max_iter = " << max_iter << " count = " << count);
  BOOST_MATH_LOG_COUNT(max_iter)
    return r;
}

template <class F, class T, class Tol>
inline std::pair<T, T> bracket_and_solve_root2(F f, const T& guess, const T& factor, bool rising, Tol tol, boost::uintmax_t& max_iter)
{
  return bracket_and_solve_root2(f, guess, factor, rising, tol, max_iter, policies::policy<>());
}

} // namespace tools
} // namespace math
} // namespace boost

double q_guess(double p,double alpha,double beta,int lower_tail,int log_p);

// Returns a quantile for the unit stable distribution.
double g_class::sqstable1(double p, int lower_tail, int log_p,
                          double dbltol, double integ_tol, int subdivisions, int verbose) {

  p_solve p_s(p, this, lower_tail, log_p,
              integ_tol, subdivisions,verbose);
  std::pair<double,double> r;
  rel_eps_tolerance tol(dbltol);
  double guess = q_guess(p,alpha,beta,lower_tail,log_p);
  if (verbose)
    Rcout << "Guess for q " << guess << std::endl;
  double factor = 2.;
  bool rising = lower_tail;
  boost::uintmax_t maxiter = 1000;
  r=boost::math::tools::bracket_and_solve_root2(p_s,guess,factor,rising,tol,maxiter);
  if (verbose)
    Rcout << "r.first = " << r.first << ", r.second - " << r.second
         << " in " << maxiter << " iterations" << std::endl;
    return (r.first+r.second)/2;
}

//Functor which contains an approximation for stable p as a function of p for
//student t with df=alpha
class pt_solve {
private:
  double rplus;
  double rminus;
  double knot1;
  double knot2;
  double a0;
  double a1;
  double a2;
  double a3;
  double p;
public:
  pt_solve(double pp,double alpha,double beta, int lower_tail, int log_p){
    double c_stable_plus = sin(M_PI_2*alpha )*tgamma(alpha)/M_PI*alpha*(1+beta);
    double c_stable_minus = sin(M_PI_2*alpha )*tgamma(alpha)/M_PI*alpha*(1-beta);
    double c_t=tgamma((alpha+1)/2)/(sqrt(alpha*M_PI)*tgamma(alpha/2))*pow(alpha,((alpha+1)/2));
    rplus=c_stable_plus/c_t;
    rminus=c_stable_minus/c_t;
    // construct a cubic spline for the mapping of pt to pstable
    students_t tdist(alpha);
    knot1=(alpha<1 && beta==1) ? cdf(tdist,-tan(M_PI_2*alpha)) : .01;
    knot2=(alpha<1 && beta==-1) ? cdf(tdist,tan(M_PI_2*alpha)) : .99;
    double dk=knot2-knot1;
    a0 = rminus*knot1;
    a1 = rminus;
    double b0 = 1-rplus*(1-knot2)-rminus*knot2;
    double b1 = rplus-rminus;
    a2 = -(dk*b1-3*b0)/(dk*dk);
    a3 = (dk*b1-2*b0)/(dk*dk*dk);
    /*    Rcout << "rminus = " << rminus << ",knot1 = " << knot1 << std::endl
               << "rplus = " << rplus << ", knot2 = " << knot2 << std::endl
               << "a0 = " << a0 << ", a1 = " << a1 << ", a2 = " << a2 << ", a3 = " << a3 << std::endl;
    */
    p = (log_p) ? exp(pp) : pp;
    p = (lower_tail) ? p : 1-p;
  }
  double operator () (double pt) {
    if (pt<=knot1)
      return pt*rminus-p;
    else if (pt>=knot2)
      return 1-rplus*(1-pt)-p;
    else {
      double ptmk1 = pt-knot1;
      double r = (((a3*ptmk1)+a2)*ptmk1+a1)*ptmk1+a0;
      return r-p;
    }
  }

};

// [[Rcpp::export]]
double guess_test(double p, double alpha, double beta, int lower_tail,int log_p){
  pt_solve pt_s(0,alpha,beta,lower_tail,log_p);
  return pt_s(p);
}

// [[Rcpp::export]]
double q_guess(double p,double alpha,double beta,int lower_tail,int log_p){
  pt_solve pt_s(p,alpha,beta,lower_tail,log_p);
  double lower = 0;
  double upper = 1;
  rel_eps_tolerance tol(1e-6);
  uintmax_t maxiter = 200;
  std::pair<double,double> r;
  r = toms748_solve(pt_s,lower,upper,tol,maxiter);
  double pt_=(r.first+r.second)/2;
  students_t tdist(alpha);
  return quantile(tdist, pt_);
}

