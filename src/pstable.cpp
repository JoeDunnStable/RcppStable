/* If using sourceCpp run Sys.setenv("PKG_CXXFLAGS"="-I /opt/local/include") */

#include <Rcpp.h>
#include "Stable.hpp"
#include <iostream>
#include <boost/math/tools/toms748_solve.hpp>
#include <R_ext/Applic.h> // For integr_fn and Rdqags

static const double Machine_eps = std::numeric_limits<double>::epsilon();

// Functor passed to toms748 solve to find x where g(x) becomes finite
class g_solve_inf {
private:
  double value;
  g_param* param;
  int verbose;

public:
  g_solve_inf(double value_in, g_param* param_in, int verbose_in) {
    value=value_in;
    param=param_in;
    verbose=verbose_in;
  }
  double operator()(const double th) {
    double g_ = g(th,param);
    if (isfinite(g_))
      return -1;
    else
      return 1;
  }
};

// Function to integrate: pstable(..)= f(..) = c2 * \int_{-\theta_0}^{\pi/2} g1_p(u) du
void g1_p(double * th, int n, void * ext) {
  g_param * param = (g_param *) ext;
  int i;

  for (i=0; i<n; i++){
    double g_th = g(th[i], param);
    // g1_p :=  exp(-g(.))
    th[i] = exp(-g_th);
  }
}

class Int_g1_p {
private:
  g_param *param;
  double abs_tol;
  double rel_tol;
  int subdivisions;
  int verbose;
  static std::string msgs[];

public:
  Int_g1_p(g_param* param_in, double tol_in, int subdivisions_in, int verbose_in){
    param=param_in;
    abs_tol=tol_in;
    rel_tol=tol_in;
    subdivisions=subdivisions_in;
    verbose=verbose_in;
  };

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
                << ", param.at0 = " << param->at0
                << ", param.cat0 = " << param->cat0
                << ", param.x_m_zet = " << param->x_m_zet << ", " << std::endl
                << "lower = " << lower
                << ", upper = " << upper << ", " <<std::endl
                << "abs_tol " << abs_tol
                << ", rel_tol " << rel_tol << std::endl;

      Rdqags(g1_p, (void *) param, &lower, &upper, &abs_tol, &rel_tol,
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
               Rcout << "Integral of g1_p from theta = " << lower
                         << " to theta = " << upper
                         << " = " << result
                         << " with absolute error = " << abserr
                         << ", subintervals = " << last << std::endl;

               return result;
  };
};

std::string Int_g1_p::msgs[]={"OK","Maximum subdivisions reached","Roundoff error detected",
                            "Bad integrand behavior","Roundoff error in the extrapolation table",
                            "Integral probably divergent","Input is invalid"};


// ------------------------------------------------------------------------------
//' Auxiliary for pstable()  (for alpha != 1)
double FCT1(double x, double zeta, double alpha, double theta0,
                   bool giveI, double dbltol, int subdivisions, int verbose)
{
  boost::uintmax_t max_iter;
  double lower;
  double upper;

  if(!isfinite(x)){
    return (giveI) ? 0 : 1;
  }
  if (!isfinite(zeta)) stop("Zeta is not finite");
  double x_m_zet = fabs(x - zeta);
//-------->>>  identically as in .fct1() for dstable() above: <<=----------
// FIXME: also provide "very small alpha" case, as in .fct1()
  if(x < zeta) theta0 = -theta0;
  double at0 = alpha*theta0;
  double cat0 = cos(at0);
  g_param param;
  param.alpha=alpha;
  param.at0=at0;
  param.cat0=cat0;
  param.x_m_zet = x_m_zet;
  if(verbose)
    Rcout << "FCT1(x =" << x <<", zeta=" << zeta << ", th0 = " << theta0
              << ", giveI = " << giveI << std::endl;

// as g() is montone, the integrand  exp(-g(.)) is too ==> maximum is at the boundary
// however, integration can be inaccuracte when g(.) quickly jumps from Inf to 0
// _BUT_  empirically I find that good values l.th / u.th below are *INDEPENDENT* of x,
//  double l_th = -theta0 + 1e-6 * fabs(-theta0);
  double l_th = -theta0;
  if(alpha > 1 && g(l_th, &param) == R_PosInf) {
    g_solve_inf g_s(0,&param, verbose);
    lower = l_th;
    upper = M_PI_2;
    max_iter = 200;
    rel_eps_tolerance tol(1e-8);
    std::pair<double,double> ur;
    ur = boost::math::tools::toms748_solve(g_s, lower, upper, tol, max_iter);
    l_th = (ur.first+ur.second)/2;
    if(verbose)
      Rcout << " g(-th0)=Inf: unirt(" << max_iter
                << " it) -> l.th=" << l_th << std::endl;
  }
//  double u_th = M_PI_2*(1-1e-6);
  double u_th = M_PI_2;
  if(alpha < 1 && g(u_th, &param) == R_PosInf) {
    g_solve_inf g_s(0,&param, verbose);
    lower = l_th;
    upper = u_th;
    max_iter = 200;
    rel_eps_tolerance tol(1e-8);
    std::pair<double,double> ur;
    ur = boost::math::tools::toms748_solve(g_s, lower, upper, tol, max_iter);
    u_th = (ur.first+ur.second)/2;
    if(verbose)
      Rcout << " g(pi/2)=Inf: unirt(" << max_iter <<" it) -> u.th= " << u_th << std::endl;

  }
  double eps = 1e-8;
  g_solve g_s(-log(eps),&param,FALSE,verbose);
  std::pair<double,double> th1_pair;
  max_iter=200;
  rel_eps_tolerance tol(1e-8);
  lower =l_th;
  upper = u_th;
  th1_pair = boost::math::tools::toms748_solve(g_s, lower, upper, tol, max_iter);
  double th1=(th1_pair.first+th1_pair.second)/2;
  if (verbose){
    double g_th1;
    g_th1=g(th1,&param);
    Rcout << "theta1 = " << th1
          << ", g(theta1) = " << g_th1 << std::endl;
    Rcout << "Region 1 is " << (th1-l_th) << " long. " << std::endl
          << "Region 2 is " << (u_th-th1) << " long. " << std::endl;
  }
  Int_g1_p int_g1_p(&param, dbltol, subdivisions, verbose);
  double r1 = int_g1_p(l_th,th1);
  double r2 = int_g1_p(th1,u_th);
  if(verbose)
    Rcout << "--> Int r1= " << r1 << std::endl
          << "    Int r2= " << r2 << std::endl;
  if(giveI) {
    // { ==> alpha > 1 ==> c1 = 1; c3 = -1/pi}
    // return (1 - F) = 1 - (1 -1/pi * r) = r/pi :
    return (r1+r2)/M_PI;
  } else {
    double c1 = (alpha < 1) ? (.5 - theta0/M_PI) : 1.;
    double c3 = (alpha < 1) ? (1/M_PI) : (-1/M_PI);
    // FIXME: for alpha > 1, F = 1 - |.|*r(x)
    //    <==> cancellation iff we eventually want 1 - F() [-> 'lower.tail']
    return c1 + c3* (r1+r2);
  }
} // {.FCT1}

// ------------------------------------------------------------------------------

typedef struct {
  double u0;
  double p2b;
  double ea;
  bool giveI;
} ga1_p_param;

//' g(), g is strictly monotone;
//' ga1_p(u) := original_g(u*pi/2) for alpha = 1
//'	 for beta > 0: increasing from g(-1) = 0   to  g(+1) = Inf
//'	 for beta < 0: decreasing from g(-1) = Inf to  g(+1) = 0
// original_g :
// g = function(th) {
//     h = p2b+ th # == g'/beta where g' := pi/2 + beta*th
//     (h/p2b) * exp(ea + h*tan(th)) / cos(th)
// }
double ga1_p(double u, ga1_p_param* param) {
  double r = u;
  if (fabs(u-param->u0) < 1e-10){
    r= 0;
    return r;
  } else {
    double th = u*M_PI_2;
    double h = param->p2b + th; // == g'/beta where g' := pi/2 + beta*th = pi/2* (1 + beta*u)
    r= (h/param->p2b) * exp(param->ea + h*tanpi2(u)) / cospi2(u);
    return r;
  }
}

// Used to find the point where ga1_p equals a finite value
class ga1_p_solve {
private:
  double value;
  ga1_p_param* param;
  bool giveI;
  int verbose;

public:
  ga1_p_solve(double value_in, ga1_p_param* param_in, int verbose_in) {
    value=value_in;
    param=param_in;
    verbose=verbose_in;
  }
  double operator()(const double th) {
    double g_ = ga1_p(th,param);
    if (verbose >=3)
      Rcout << "ga1_p_solve: p2b= " << param->p2b << ", ea = "<< param->ea
            << ", u0 = " << param->u0 << std::endl
            << "th= " << th << ", g(th) = "<< g_ << " vs target " << value << std::endl;
    return g_-value;
  }
};

// Used to find the point where ga1_p becomes infinite
class ga1_solve_inf {
private:
  double value;
  ga1_p_param* param;
  bool giveI;
  int verbose;

public:
  ga1_solve_inf(double value_in, ga1_p_param* param_in, int verbose_in) {
    value=value_in;
    param=param_in;
    verbose=verbose_in;
  }
  double operator()(const double th) {
    double g_ = ga1_p(th,param);
    if (verbose >=3)
    Rcout << "ga1_solve_inf: p2b= " << param->p2b << ", ea = "<< param->ea
              << ", u0 = " << param->u0 << std::endl
              << "th= " << th << ", g(th) = "<< g_
              << ", isfinite(g_) " << isfinite(g_) << std::endl;
    if (isfinite(g_))
      return -1;
    else
      return 1;
  }
};

// Function to Integrate; u is a non-sorted vector!
void g2_p(double*u, int n, void* ext) {
  ga1_p_param* param = (ga1_p_param*) ext;
  // g2_p = exp(-g(.))
  for (int i=0; i<n; i++)
    if (param->giveI)
      u[i]=expm1(- ga1_p(u[i],param));
    else
      u[i]=exp(- ga1_p(u[i],param));
}

//Public version of ga1_p.  Used for debugging
// [[Rcpp::export]]
NumericVector ga1_p_public(double u,double x,double beta) {
  double i2b = 1/(2*beta);
  double p2b = M_PI*i2b; // = M_PI/(2 beta)

  double ea = -p2b* (x);
  if(!isfinite(ea)){
    double r = (ea < 0) ? R_PosInf : 0;
    return wrap(r);
  }

  //t0 = -pi2# g(t0) == 0  mathematically, but not always numerically
  double u0 = -1; // g(u0) == 0  mathematically, but not always numerically
  ga1_p_param param;
  param.p2b=p2b;
  param.ea=ea;
  param.giveI=FALSE;
  param.u0=u0;
  double r = ga1_p(u,&param);
  return wrap(r);
}

  class Int_g2_p {
private:
  ga1_p_param *param;
  double abs_tol;
  double rel_tol;
  int subdivisions;
  int verbose;
  static std::string msgs[];

public:
  Int_g2_p(ga1_p_param* param_in, double tol_in, int subdivisions_in, int verbose_in){
    param=param_in;
    abs_tol=tol_in;
    rel_tol=tol_in;
    subdivisions=subdivisions_in;
    verbose=verbose_in;
  };

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
                << "Rdqags(g2_p" << "with param.u0 = " << param->u0
                << ", param.p2b = " << param->p2b
                << ", param.ea = " << param->ea
                << "lower = " << lower
                << ", upper = " << upper << ", " <<std::endl
                << "abs_tol " << abs_tol
                << ", rel_tol " << rel_tol << std::endl;

      Rdqags(g2_p, (void *) param, &lower, &upper, &abs_tol, &rel_tol,
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
               Rcout << "Integral of g2_p from theta = " << lower
                         << " to theta = " << upper
                         << " = " << result
                         << " with absolute error - " << abserr << std::endl;

               return result;
  };
};

std::string Int_g2_p::msgs[]={"OK","Maximum subdivisions reached","Roundoff error detected",
                            "Bad integrand behavior","Roundoff error in the extrapolation table",
                            "Integral probably divergent","Input is invalid"};

//' Auxiliary for spstable()  only used when alpha == 1 :
//' @param x numeric *scalar*
//' @param beta  >= 0 here
//' @param tol
//' @param subdivisions
double FCT2(double x, double beta, double dbltol, int subdivisions,
                   bool giveI, int verbose)
{
  Rcout.precision(12);
  boost::uintmax_t max_iter;
  double lower;
  double upper;
  double i2b = 1/(2*beta);
  double p2b = M_PI*i2b; // = M_PI/(2 beta)

  double ea = -p2b* (x);
  if(!isfinite(ea)){
    double p = (ea < 0) ? 1 : 0;
    return (!giveI) ? p : 1-p;
  }

  //t0 = -pi2# g(t0) == 0  mathematically, but not always numerically
  double u0 = -1; // g(u0) == 0  mathematically, but not always numerically
  ga1_p_param param;
  param.p2b=p2b;
  param.ea=ea;
  param.giveI=giveI;
  param.u0=u0;
  if(verbose)
    Rcout << ".FCT2(x= " << x <<", beta = " << beta
              << " , giveI =" << giveI <<"): ";


// g(-u0) == +Inf {at other end}, mathematically ==> exp(-g(.)) == 0
// in the outer tails, the numerical integration can be inaccurate,
// because g(.) jumps from 0 to Inf,  but is 0 almost always
//   <==> g1_p(.) = exp(-g(.)) jumps from 1 to 0 and is 1 almost everywhere
//  ---> the integration "does not see the 0" and returns too large..
  double u_ = 1;
//  double uu = u_* (1-1e-6);
  double uu = u_;
  if(ga1_p(uu,&param)== R_PosInf) {
    ga1_solve_inf ga1_s_inf(0,&param, verbose);
    lower = -1;
    upper = uu;
    max_iter = 200;
    abs_eps_tolerance tol(1e-8);
    std::pair<double,double> ur;
    // Determine the point where ga1_p becomes infinite
    ur = boost::math::tools::toms748_solve(ga1_s_inf, lower, upper, tol, max_iter);
    u_ = ur.first;  // Take the finite side
    if(verbose)
      Rcout << " g(" << uu << ")=Inf:" << std::endl <<
                " unirt(" << max_iter << " it) -> u.first = " << ur.first
                          << ", ur.second = " << ur.second << std::endl;
  }
  double eps=1e-4;
  bool do1;
  double th1;
  if((do1=(ga1_p(-1,&param)<eps && ga1_p(u_,&param)>eps))) {
    ga1_p_solve ga1_p_s(eps,&param, verbose);
    lower = -1;
    upper = u_;
    max_iter = 200;
    abs_eps_tolerance tol(1e-8);
    std::pair<double,double> th1_pair;
    // Determine the point where ga1_p passes through eps
    th1_pair = boost::math::tools::toms748_solve(ga1_p_s, lower, upper, tol, max_iter);
    th1 = (th1_pair.first + th1_pair.second)/2;
    if(verbose)
      Rcout << " g(" << th1 << ") = "<< th1 << std::endl <<
        " unirt(" << max_iter << " it) -> th1_pair.first = " << th1_pair.first
                  << ", th1_pair.second = " << th1_pair.second << std::endl;
  }


//' g2_p(.) = exp(-g(.)) is strictly monotone .. no need for 'theta2' !
//'
  Int_g2_p int_g2_p(&param, dbltol, subdivisions, verbose);
  double r1, r2, r3;
  if (do1) {
    r1 = int_g2_p(-1,th1);
    r2 = int_g2_p(th1,u_);
  } else {
    r1=0;
    r2 = int_g2_p(-1, u_);
  }
  // JLD: Added regions r1 and r2 are to force integrate to recognize the important region
  // region 3 contributes to the total when giveI=TRUE & integrand = 1-exp(-g())
  r3 = int_g2_p(u_,1);
  double r = (r1+r2+r3)/2;
  if(verbose)
    Rcout << "--> Int r= " << r << std::endl;
  return (giveI) ? -r : r;
}// {.FCT2}

inline double retValue(double F, int useLower, int log_p) {
  return (useLower) ? ((log_p) ? log(F)    : F)
                    : ((log_p) ? log1p(-F) : 1 - F);

}



double spstable1(double z, double alpha, double beta, int lower_tail, int log_p,
                       double tol, int subdivisions, int verbose) {
  if (alpha !=1) {
    double tanpa2 = tan(M_PI_2*alpha);
    double zeta = -beta * tanpa2;
    double theta0 = fmin(fmax(-M_PI_2, atan(-zeta) / alpha), M_PI_2);
    bool finSupp = (fabs(beta) == 1 && alpha < 1);
    if(finSupp) {
      // has *finite* support    [zeta, Inf)  if beta ==  1
      //                         (-Inf, zeta] if beta == -1
      if(beta == 1 && z <= zeta)
        return retValue(0., lower_tail, log_p);
      else if(beta == -1 && z >= zeta)
        return retValue(1., lower_tail, log_p);
    // else .. one of the cases below
    }

    if(fabs(z - zeta) < 2 * Machine_eps) {
      // FIXME? same problem as dstable
      double r = (lower_tail) ? (.5 - theta0/M_PI) : (.5 + theta0/M_PI);
      if (verbose)
        Rcout << "z ~ zeta: " << ".5 "
                  << ((lower_tail) ? " - " : " + ") << "th0/pi = " << r << std::endl;
      return (log_p) ? log(r) : r;
    }
    int useLower = ((z > zeta && lower_tail) ||
                 (z < zeta && !lower_tail));
    // FIXME: for alpha > 1 -- the following computes F1 = 1 -c3*r(x)
    // and suffers from cancellation when 1-F1 is used below:
    bool giveI = !useLower && alpha > 1; // if TRUE, .FCT1() return 1-F
    double _F1 = FCT1(z, zeta, alpha, theta0,
                      giveI, tol, subdivisions, verbose);
    if(giveI)
      return (log_p) ? log(_F1) : _F1;
    else
      return retValue(_F1, useLower, log_p);
  } else { // alpha = 1
    bool useL;
    if(beta >= 0)
      useL = lower_tail;
    else {
      beta = -beta;
      z = -z;
      useL=!lower_tail;
    }
    bool giveI = !useL && !log_p;
    if(giveI)
      useL = TRUE;
    double _F2 = FCT2(z,beta, tol, subdivisions, giveI, verbose);
    return retValue(_F2,useL,log_p);
  }
}

// Functor passed to toms748 solve to find q for unit stable distribution
class p_solve {
private:
  double p;
  double alpha;
  double beta;
  int lower_tail;
  int log_p;
  int maxiter;
  int verbose;
  double integ_tol;
  int subdivisions;

public:
  p_solve(double p_in, double alpha_in, double beta_in,
          int lower_tail_in, int log_p_in,
          double integ_tol_in, int subdivisions_in, int verbose_in) {
    p=p_in;
    alpha=alpha_in;
    beta=beta_in;
    lower_tail=lower_tail_in;
    log_p=log_p_in;
    verbose=verbose_in;
    integ_tol=integ_tol_in;
    subdivisions=subdivisions_in;
  }
  double operator()(const double q) {
    if (verbose)
      Rcout << "Calling spstable with parmerters" << std::endl
               << "q = " << q << ", alpha = " << alpha << ", beta = " << beta << std::endl
               << "targeting p = " << p;

    double r = spstable1(q, alpha, beta, lower_tail, log_p,
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
