/* If using sourceCpp run Sys.setenv("PKG_CXXFLAGS"="-I /opt/local/include") */

#include <Rcpp.h>
#include "Stable.h"
#include <iostream>
#include <boost/math/tools/toms748_solve.hpp>
#include <R_ext/Applic.h> // For integr_fn and Rdqags

static const double Machine_eps = std::numeric_limits<double>::epsilon();

// Function to integrate: pstable(..)= f(..) = c2 * \int_{-\theta_0}^{\pi/2} g1_p(u) du
void exp_m_g(double * th, int n, void * ext) {
  g_class * param = (g_class *) ext;
  int i;

  for (i=0; i<n; i++){
    double g_th = param->g(th[i]);
    // g1_p :=  exp(-g(.))
    th[i] = exp(-g_th);
  }
}

class Int_exp_m_g {
private:
  g_class *param;
  double abs_tol;
  double rel_tol;
  int subdivisions;
  int verbose;
  static std::string msgs[];

public:
  Int_exp_m_g(g_class* param, double tol, int subdivisions, int verbose) :
    param(param), abs_tol(tol), rel_tol(tol), subdivisions(subdivisions), verbose(verbose) {};

  double result;
  double abserr;
  int neval;
  int ier;
  double operator() (double lower, double upper) {
    int lenw = 4*subdivisions;
    int last;
    int iwork[subdivisions];
    double work[lenw];

    if (verbose>=4)
      Rcout << std::endl
                << "Rdqags(g1_p" << "with param.alpha = " << param->alpha
                << ", param.beta = " << param->beta
                << ", param.cat0 = " << param->cat0
                << ", param.x_m_zet = " << param->x_m_zet << ", " << std::endl
                << "lower = " << lower
                << ", upper = " << upper << ", " <<std::endl
                << "abs_tol " << abs_tol
                << ", rel_tol " << rel_tol << std::endl;

      Rdqags(exp_m_g, (void *) param, &lower, &upper, &abs_tol, &rel_tol,
             &result, &abserr, &neval, &ier,
             &subdivisions, &lenw, &last, iwork, work);
             if (ier>0) {
               if (ier<6){
                 warning(msgs[ier]);
               } else {
                 stop(msgs[ier]);
               }
             }
             if (verbose>=3)
               Rcout << "Integral of exp_m_g from theta = " << lower
                         << " to theta = " << upper
                         << " = " << result
                         << " with absolute error = " << abserr
                         << ", subintervals = " << last << std::endl;

             return result;
  };
};

std::string Int_exp_m_g::msgs[]={"OK","Maximum subdivisions reached","Roundoff error detected",
                            "Bad integrand behavior","Roundoff error in the extrapolation table",
                            "Integral probably divergent","Input is invalid"};


// ------------------------------------------------------------------------------
//' Auxiliary for pstable()  (for alpha != 1)
double g_class::integrate_exp_m_g(bool giveI, double dbltol, int subdivisions, int verbose)
{
  boost::uintmax_t max_iter;
  double lower;
  double upper;

  if(verbose)
    Rcout << "integrate_exp_m_g(" << *this << std::endl << "giveI = " << giveI << std::endl;

// as g() is montone, the integrand  exp(-g(.)) is too ==> maximum is at the boundary
// however, integration can be inaccuracte when g(.) quickly jumps from Inf to 0
  double eps=1e-10;
  double th_eps, th_1_m_eps, th1, th2;
  double exp_m_g_lo=0;
  exp_m_g(&exp_m_g_lo, 1, (void*)this);
  double exp_m_g_hi=th_max;
  exp_m_g(&exp_m_g_hi, 1, (void*)this);
  if (verbose)
    Rcout << "integrate_exp_m_g: exp(-g(0) = " << exp_m_g_lo << std::endl
          << ", exp(-g(" << th_max << ") = " << exp_m_g_hi << std::endl;
  bool have_eps = FALSE, have_1_m_eps = FALSE;
  bool do1=TRUE, do2=TRUE, do3=TRUE;
  if ((exp_m_g_hi-eps)*(exp_m_g_lo-eps)<0){
    g_solve g_s(-log(eps), this, verbose);
    lower = 0;
    upper = th_max;
    max_iter = 200;
    rel_eps_tolerance tol(1e-8);
    std::pair<double,double> ur;
    ur = boost::math::tools::toms748_solve(g_s, lower, upper, tol, max_iter);
    th_eps = (ur.first+ur.second)/2;
    have_eps=TRUE;
    if(verbose)
      Rcout << " exp(-g)(" << th_eps <<") = " << eps << ", iterations = " << max_iter << std::endl;
  }
  if ((exp_m_g_hi-(1-eps))*(exp_m_g_lo-(1-eps))<0){
    g_solve g_s(-log(1-eps), this, verbose);
    lower = (exp_m_g_lo==0) ? th_eps : 0;
    upper = (exp_m_g_hi==0) ? th_eps : th_max;
    max_iter = 200;
    rel_eps_tolerance tol(1e-8);
    std::pair<double,double> ur;
    ur = boost::math::tools::toms748_solve(g_s, lower, upper, tol, max_iter);
    th_1_m_eps = (ur.first+ur.second)/2;
    have_1_m_eps=TRUE;
    if(verbose)
      Rcout << "exp(-g)(" << th_1_m_eps <<") = " << 1-eps << ", iterations = " << max_iter << std::endl;
  }
  if (have_1_m_eps){
    if (exp_m_g_hi==0) {
      th1=th_1_m_eps;
      th2=th_eps;
    } else {
      th1=th_eps;
      th2=th_1_m_eps;
    }
  } else if (have_eps){
    if (exp_m_g_hi==0){
      do1=FALSE;
      th1=0;
      th2=th_eps;
    } else {
      do3=FALSE;
      th1=th_eps;
      th2=th_max;
    }
  } else {
    do1=FALSE;
    do2=FALSE;
    th2=0;
  }
  Int_exp_m_g int_exp_m_g(this, dbltol, subdivisions, verbose);
  double r1, r2, r3;
  double rerr = 0;
  if (do1) {
    r1=int_exp_m_g(0,th1);
    rerr+=int_exp_m_g.abserr;
  }
  else
    r1=0;
  if (do2) {
    r2=int_exp_m_g(th1,th2);
    rerr+=int_exp_m_g.abserr;
  } else
    r2=0;
  if (do3){
    r3=int_exp_m_g(th2,th_max);
    rerr+=int_exp_m_g.abserr;
  } else
    r3=0;
  if(verbose)
    Rcout << "--> Int r1= " << r1 << std::endl
          << "    Int r2= " << r2 << std::endl
          << "    Int r3= " << r3 << std::endl;
  double c0, c1;
  if (alpha < 1){
    c0 = .5 - theta0/M_PI;
    c1 = 1/M_PI;
  } else if (alpha > 1){
    c0=1;
    c1=-1/M_PI;
  } else{ // alpha == 1
    c0 = 1;
    c1 = -.5;
  }
  if(giveI) {
    // { ==> alpha > 1 ==> c0 = 1; c1 = -1/pi}
    // return (1 - F) = 1 - (1 -1/pi * r) = r/pi :
    return -(r1+r2+r3)*c1;
  } else {
    return c0 + c1* (r1+r2+r3);
  }
} // {.FCT1}

// ------------------------------------------------------------------------------


inline double retValue(double F, int useLower, int log_p) {
  return (useLower) ? ((log_p) ? log(F)    : F)
                    : ((log_p) ? log1p(-F) : 1 - F);

}

double g_class::spstable1(double x,
                 int lower_tail, int log_p,
                 double tol, int subdivisions, int verbose) {
  if(x==R_PosInf)
    return retValue(1, lower_tail, log_p);
  else if (x==R_NegInf)
    return retValue(0, lower_tail, log_p);
  set_x(x);
  if (alpha !=1) {
    //-------->>>  identically as in .fct1() for dstable() above: <<=----------
    // FIXME: also provide "very small alpha" case, as in .fct1()
    if (!isfinite(zeta)) stop("Zeta is not finite");
    bool finSupp = (fabs(beta) == 1 && alpha < 1);
    if(finSupp) {
      // has *finite* support    [zeta, Inf)  if beta ==  1
      //                         (-Inf, zeta] if beta == -1
      if(beta_input == 1 && x <= zeta)
        return retValue(0., lower_tail, log_p);
      else if(beta_input == -1 && x >= zeta)
        return retValue(1., lower_tail, log_p);
    // else .. one of the cases below
    }

    if(fabs(x - zeta) < 2 * Machine_eps) {
      // FIXME? same problem as dstable
      double r = (lower_tail) ? (.5 - theta0_x_gt_zeta/M_PI) : (.5 + theta0_x_gt_zeta/M_PI);
      if (verbose)
        Rcout << "x ~ zeta: " << ".5 "
                  << ((lower_tail) ? " - " : " + ") << "th0/pi = " << r << std::endl;
      return (log_p) ? log(r) : r;
    }
    int useLower = ((x > zeta && lower_tail) ||
                 (x < zeta && !lower_tail));
    // FIXME: for alpha > 1 -- the following computes F1 = 1 -c3*r(x)
    // and suffers from cancellation when 1-F1 is used below:
    bool giveI = !useLower && alpha > 1; // if TRUE, int_exp_m_g() return 1-F
    double _F1 = integrate_exp_m_g(giveI, tol, subdivisions, verbose);
    if(giveI)
      return (log_p) ? log(_F1) : _F1;
    else
      return retValue(_F1, useLower, log_p);
  } else { // alpha = 1
    bool useL;
    if(beta_input>=0)
      useL = !lower_tail;
    else {
      useL=lower_tail;
    }

    bool giveI = !useL && !log_p;
    if(giveI)
      useL = TRUE;
    double _F2 = integrate_exp_m_g(giveI, tol, subdivisions, verbose);
    return retValue(_F2,useL,log_p);
  }
}

// Functor passed to toms748 solve to find q for unit stable distribution
class p_solve {
private:
  double p;
  g_class param;
  int lower_tail;
  int log_p;
  int maxiter;
  int verbose;
  double integ_tol;
  int subdivisions;

public:
  p_solve(double p, double alpha, double beta,
          int lower_tail, int log_p,
          double integ_tol, int subdivisions, int verbose) :
    p(p), param(alpha, beta), lower_tail(lower_tail) , log_p(log_p),
    verbose(verbose), integ_tol(integ_tol), subdivisions(subdivisions){}
  double operator()(const double q) {
    if (verbose)
      Rcout << "Calling spstable with parmerters" << std::endl
            << "q = " << q << ", alpha = " << param.alpha << ", beta = " << param.beta << std::endl
               << "targeting p = " << p;

    double r = param.spstable1(q, lower_tail, log_p,
                     integ_tol, subdivisions, FALSE)-p;
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
double sqstable1(double p, double alpha, double beta, int lower_tail, int log_p,
                 double dbltol, double integ_tol, int subdivisions, int verbose) {

  p_solve p_s(p, alpha, beta, lower_tail, log_p,
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
    knot1=(alpha<1 && beta==1) ? Rf_pt(-tan(M_PI_2*alpha),alpha,1,0) : .01;
    knot2=(alpha<1 && beta==-1) ? Rf_pt(tan(M_PI_2*alpha),alpha,1,0) : .99;
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
  r = boost::math::tools::toms748_solve(pt_s,lower,upper,tol,maxiter);
  double pt_=(r.first+r.second)/2;
  return Rf_qt(pt_,alpha,TRUE,FALSE);
}
