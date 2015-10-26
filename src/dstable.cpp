/* If using sourceCpp run Sys.setenv("PKG_CXXFLAGS"="-I /opt/local/include") */

#include <Rcpp.h>
#include "Stable.h"
#include <Rmath.h>
#include <limits>
#include <iostream>
#include <boost/math/tools/toms748_solve.hpp>
#include <R_ext/Applic.h> // For integr_fn and Rdqags

/* This is an ancillary functions use by the R function fct1 */
double dstable_smallA(double x, double alpha, double beta, bool log_flag=FALSE);

static const double Machine_eps = std::numeric_limits<double>::epsilon();
static const double alpha_small_dstable = 1e-17;
static const double large_exp_arg = 708.39641853226407875e0;

/* The input parameters for g other than theta */
integr_fn g1;

double g(double th, g_param* param) {
// g(-pi/2) or g(pi/2) could become  NaN --> work around
  if (fabs(M_PI_2 - ((param->alpha > 1)? th :-th)) < 64 * Machine_eps) return 0.;
  double att = param->at0 + param->alpha * th; // = alpha*(theta0 + theta)
  double r = pow(param->cat0 * cos(th) * pow(param->x_m_zet/sin(att),param->alpha),1/(param->alpha-1)) * cos(att - th);
  if (isnan(r) || r<0){
    if ((param->alpha>1 && fabs(th-M_PI_2)<1e-5) ||(param->alpha<1 && fabs(att)<1e-5))
      r=0;
    else if ((param->alpha>1 && fabs(att)<1e-5) ||(param->alpha <1 && fabs(th-M_PI_2)<1e-5))
      r=R_PosInf;
    else
      r=NAN;
  }
  return r;
}

//' @title  x*exp(-x)  numerically stable, with correct limit 0 for x --> Inf
//' @param x  numeric
//' @return x*exp(x)
//' @author Martin Maechler
inline double x_exp_m_x(double x) {
  double r;
  if(isnan(x))
    r = NAN;
  else if(x > large_exp_arg) // e.g. x == Inf
    r = 0;
  else
    r = x*exp(-x);
  return r;
}

// Function to integrate: dstable(..)= f(..) = c2 * \int_{-\theta_0}^{\pi/2} g1(u) du
void g1(double * th, int n, void * ext) {
  g_param * param = (g_param *) ext;
  int i;

  for (i=0; i<n; i++){
    double g_th = g(th[i], param);
    // g1 :=  g(.) exp(-g(.))
    th[i] = x_exp_m_x(g_th);
  }
}

// This class contains everything need to integrate g1 except the endpoints
// of the interval.  The operator() double(a,b) actually performs the integration
// using the R version of quadpack, Rdqags, which is exported from the R interface
class Int_g1 {
private:
  g_param *param;
  double abs_tol;
  double rel_tol;
  int subdivisions;
  int verbose;
  static std::string msgs[];

public:
  Int_g1(g_param* param_in, double tol_in, int subdivisions_in, int verbose_in){
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
                << "Rdqags(g1" << "with param.alpha = " << param->alpha
                << ", param.at0 = " << param->at0
                << ", param.cat0 = " << param->cat0
                << ", param.x_m_zet = " << param->x_m_zet << ", " << std::endl
                << "lower = " << lower
                << ", upper = " << upper << ", " <<std::endl
                << "abs_tol " << abs_tol
                << ", rel_tol " << rel_tol << std::endl;

    Rdqags(g1, (void *) param, &lower, &upper, &abs_tol, &rel_tol,
           &result, &abserr, &neval, &ier,
           &subdivisions, &lenw, &last, iwork, work);
    if (ier>0) {
      if (ier<6){
         double endpt[] = {lower,upper};
         g1(endpt,2,(void *)param);
         // The min and max of the result given that the integrand is monotonic
         double r_min=std::min(endpt[0],endpt[1])*(upper-lower);
         double r_max=std::max(endpt[0],endpt[1])*(upper-lower);
         if (result < r_min) {
           result = r_min;
           abserr = fabs(endpt[1]-endpt[0])*(upper-lower);
           warning("result was less than indicated by endpoints\n");
         } else if (result > r_max) {
           result = r_max;
           abserr = fabs(endpt[1]-endpt[0])*(upper-lower);
           warning("result was greater than indicated by endpoints\n");
         }
         warning(msgs[ier]);
      } else {
        stop(msgs[ier]);
      }
    }
    if (verbose>=3)
      Rcout << "Integral of g1 from theta = " << lower
                << " to theta = " << upper
                << " = " << result
                << ", with absolute error = " << abserr
                << ", subintervals = " << last << std::endl;

    return result;
  };
};

std::string Int_g1::msgs[]={"OK","Maximum subdivisions reached","Roundoff error detected",
        "Bad integrand behavior","Roundoff error in the extrapolation table",
        "Integral probably divergent","Input is invalid"};

double fct1(double x, double zeta,
                    double alpha, double beta, double theta0,
                    bool log_flag, double dbltol, int subdivisions,
                    double zeta_tol, int verbose) {
  // --- dstable(x, alpha, beta, ..)  for alpha < 2 ---
  // For  x = zeta, have special case formula [Nolan(1997)];
  // need to use it also for x ~= zeta, i.e., x.m.zet := |x - zeta| < delta :
  if (!isfinite(zeta)) stop("!is_finite(zeta)");
  boost::uintmax_t max_iter;
  rel_eps_tolerance tol(dbltol);

  double x_m_zet = fabs(x - zeta);
  // For  x = zeta, have special case formula [Nolan(1997)];
  // need to use it also for x ~= zeta, i.e., x.m.zet := |x - zeta| < delta :
  double f_zeta = (log_flag?
      lgamma(1 + 1/alpha) + log(cos(theta0)) - (log(M_PI) + log1p(pow(zeta,2))/(2 * alpha)):
      tgamma(1 + 1/alpha) * cos(theta0)/(M_PI * pow(1 + pow(zeta,2),(1/(2 * alpha)))));
// Modified: originally was	 if (z == zeta),
// then (D.W.)   if (x.m.zet < 2 * .Machine$double.eps)
// then (M.M.)   if (x.m.zet <= 1e-5 * abs(x))
  if (isfinite(x) && x_m_zet <= zeta_tol * (zeta_tol + std::max(fabs(x),
                                                 fabs(zeta)))) {
    if (verbose)
      Rcout << ".fct1(" << x
                << "~=" << zeta
                << " Using f.zeta()" << std::endl;
    return f_zeta;
  }
  // the real check should be about the feasibility of g() below, or its integration

  bool smallAlpha = alpha < alpha_small_dstable;
  if (x < zeta) {
    theta0 = -theta0; // see Nolan(1997), Thm.1 (c)
    if (smallAlpha) {
      beta = -beta;
      x = -x;
    }
  }
  if (smallAlpha) {
    // here, *MUST* have  __ x > zeta __
    if (verbose)
      Rcout << ".fct1(" << x << " , "
                << zeta
                << ",..): small alpha=" << alpha << "\n";
    return dstable_smallA(x, alpha, beta, log_flag);
  }
// constants ( independent of integrand g1(th) = g*exp(-g) ):
// zeta = -beta * tan(pi*alpha/2)
// theta0 = (1/alpha) * atan( beta * tan(pi*alpha/2))
// x.m.zet = abs(x - zeta)
//-------->>>  identically as in .FCT1() for pstable() below: <<=----------
  double a_1 = alpha - 1;
  double at0 = alpha * theta0;
  double cat0 = cos(at0);
  g_param param;
  param.alpha=alpha;
  param.at0=at0;
  param.cat0=cat0;
  param.x_m_zet=x_m_zet;
  //' g() is strictly monotone -- Nolan(1997) ["3. Numerical Considerations"]
  //'     alpha >= 1  <==>  g() is falling, ie. from Inf --> 0;  otherwise growing from 0 to +Inf
  double c2 = (alpha/(M_PI * fabs(a_1) * x_m_zet));
  double g_pi2 = g(M_PI_2,&param);
  double g_mth0 = g(-theta0,&param);
  if ((alpha >= 1 && (g_pi2 > large_exp_arg || g_mth0 == 0)) ||
      (alpha < 1 && (g_mth0 > large_exp_arg || g_pi2  == 0))) {
    if (verbose)
      Rcout << ".fct1(" << x
                << " , "<< zeta
                <<  ",..): g() is 'Inf' (or 0) ==> result 0" << std::endl;
    if (log_flag)
      return R_NegInf;
    else
      return 0;
  }
  double g_;
  if (alpha >= 1)
    g_=g(-theta0+1e-06*fabs(theta0),&param);
  else
    g_=g(M_PI_2*(1-1e-06),&param);
  if (isnan(g_)) {
    Rcout << "g(" << ((alpha>=1)?-theta0+1e-06*fabs(theta0):M_PI_2*(1-1e-6)) << ") is NaN" << std::endl;
    if (fmax(x_m_zet, x_m_zet/fabs(x)) < 0.01){
      if (verbose)
        Rcout << "x_m_zet is close to zero.  Using f.zeta" << std::endl;
      return f_zeta;
    }
  }
  g_=g(M_PI_2,&param);
  if (verbose >=2)
    Rcout << std::endl
          << "-theta0 = "<< -theta0 << ", g(-theta0) = " << g_mth0 << std::endl
          << "pi/2 =    " << M_PI_2 << ", g(M_PI_2) =  " << g_ << std::endl;
  double theta2;
  double g1_th2;
  if ((alpha >= 1 && !isnan(g_) && g_ > 1) ||
      (alpha < 1 && !isnan(g_) && g_ < 1)) {
    theta2 = M_PI_2-1e-06;
    g1_th2 = theta2;
    g1(&g1_th2,1,&param);
  }
  else {
    if ((alpha < 1 && g_mth0 > 1) ||
        (alpha >= 1 && g_mth0 < 1)){
      theta2 = -theta0  +1e-06;
      g1_th2=theta2;
      g1(&g1_th2,1,&param);
    } else {
      double l_th = -theta0;
      double u_th = M_PI_2;
      g_solve g_s(1.,&param,FALSE,verbose);
      std::pair<double, double> ur1_pair;
      // toms748_solve passes lower, upper and max_iter by reference and changes them.
      max_iter = 200;
      double upper=u_th;
      double lower=l_th;
      ur1_pair = boost::math::tools::toms748_solve(g_s, lower, upper, tol, max_iter);
      double ur1 = (ur1_pair.first+ur1_pair.second)/2;
/*      g_solve log_g_s(0,&param,TRUE,verbose);
      std::pair<double, double> ur2_pair;
      // The R code used a try catch for the next line
      max_iter=200;
      upper=u_th;
      lower=l_th;
      ur2_pair = boost::math::tools::toms748_solve(log_g_s, lower, upper, tol, max_iter);
      double ur2 = (ur2_pair.first+ur2_pair.second)/2;
*/
      double g_1 = ur1;
      g1(&g_1,1,&param);
/*      double g_2 = ur2;
      g1(&g_2,1,&param);
      if (fabs(g_1-g_2)>64 * Machine_eps)
        Rcout << std::endl << "root1 pair = " << ur1_pair.first
                  << ", " << ur1_pair.second << std::endl
                  << "root2.pair = " << ur2_pair.first
                  << ", " <<  ur2_pair.second << std::endl
                  << "root1 = " << ur1
                  << ", g1(root1) = " << g_1 << std::endl
                  << "root2 = " << ur2
                  << ", g1(root2) = " << g_2 << std::endl;
      if (g_1 >= g_2) {
*/
        theta2 = ur1;
        g1_th2 = g_1;
/*      }
      else {
        theta2 = ur2;
        g1_th2 = g_2;
      }
*/
}
  }
  double eps = 1e-10;
  double g1_th0= -theta0;
  g1(&g1_th0,1,&param);
  g1_solve g1_s(eps,&param,verbose);
  bool do1;
  double th1;
  if (verbose>=2)
    Rcout << std::endl << "theta2 = " << theta2
             << ", g1(theta2) = " << g1_th2 << std::endl;
  if ((do1 = g1_th2 > eps && g1_th0 < eps)) {
    std::pair<double,double> th1_pair;
    max_iter=200;
    double lower =-theta0;
    double upper = theta2;
    th1_pair = boost::math::tools::toms748_solve(g1_s, lower, upper, tol, max_iter);
    th1=(th1_pair.first+th1_pair.second)/2;
    if (verbose>=2){
      double g1_th1 = th1;
      g1(&g1_th1,1,&param);
      Rcout << "theta1 = " << th1
                << ", g1(theta1) = " << g1_th1 << std::endl;
      Rcout << "Region 1 is " << (th1+theta0) << " long. " << std::endl
            << "Region 2 is " << (theta2-th1) << " long. " << std::endl;
    }
  }
  double g1_pi2=M_PI_2;
  g1(&g1_pi2,1,&param);
  double th3;
  bool do4;
  if ((do4 = g1_th2 > eps && g1_pi2 < eps)) {
    std::pair<double,double> th3_pair;
    max_iter=200;
    double lower = theta2;
    double upper = M_PI_2;
    th3_pair = boost::math::tools::toms748_solve(g1_s, lower, upper, tol, max_iter);
    th3=(th3_pair.first+th3_pair.second)/2;
    if (verbose>=2){
      double g1_th3 = th3;
      g1(&g1_th3,1,&param);
      Rcout << "theta3 = " << th3
                << ", g1(theta3) = " << g1_th3 << std::endl;
      Rcout << "Region 3 is " << (th3-theta2) << " long. " << std::endl
            << "Region 4 is " << (M_PI_2-th3) << " long. " << std::endl;
    }
  }
  double r1, r2, r3, r4;
  double rerr=0;
  Int_g1 int_g1(&param, dbltol, subdivisions, verbose);
  if (do1) {
    r1 = int_g1(-theta0, th1);
    rerr=rerr+int_g1.abserr;
    r2 = int_g1(th1, theta2);
    rerr=rerr+int_g1.abserr;
  }
  else {
    r1 = 0;
    r2 = int_g1(-theta0, theta2);
    rerr=rerr+int_g1.abserr;
  }
  if (do4) {
    r3 = int_g1(theta2, th3);
    rerr=rerr+int_g1.abserr;
    r4 = int_g1(th3, M_PI_2);
    rerr=rerr+int_g1.abserr;
  }
  else {
    r3 = int_g1(theta2, M_PI_2);
    rerr=rerr+int_g1.abserr;
    r4 = 0;
  }
  if (verbose) {
    if (verbose)
      Rcout << ".fct1(" << x << " , "
                <<  zeta << ",..): c2*sum(r[1:4])= "
                << c2 << "*"
                << "(" << r1
                << " + " << r2
                << " + " << r3
                << " + " << r4
                << ") = " << c2*(r1+r2+r3+r4)
              << ", abs.err = " << c2*rerr << std::endl;
  }
  if (rerr > .25*fabs((r1+r2+r3+r4)))
    return (log_flag) ? R_NegInf : 0;  ///this will cause the use of dPareto.
  if (log_flag)
    return log(c2) + log(r1 + r2 + r3 + r4);
  else ;
    return c2 * (r1 + r2 + r3 + r4);
} //fct1

//' dstable() for very small alpha > 0
//' ok only for  x > zeta := - beta * tan(pi/2 *alpha)
double dstable_smallA(double x, double alpha, double beta, bool log_flag) {
  double r = log(alpha)+log1p(beta)-(1+log(2*x+M_PI*alpha*beta));
  if (log_flag)
    return r;
  else
    return exp(r);
}

double ga1(double u, ga1_param* param) {
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

// Function to Integrate; u is a non-sorted vector!
void g2(double*u, int n, void* ext) {
  ga1_param* param = (ga1_param*) ext;
  // g2 = g(.) exp(-g(.))
  for (int i=0; i<n; i++)
    u[i]=x_exp_m_x( ga1(u[i],param));
}

// ------------------------------------------------------------------------------

// Functor passed to toms748_solve to find x such that g1(x) == value
class g2_solve {
private:
  double value;
  ga1_param *param;
  int verbose;
public:
  g2_solve(double value_in,ga1_param *param_in, int verbose_in) {
    value=value_in;
    param=param_in;
    verbose=verbose_in;
  }
  double operator()(const double u) {
    double g2_=u;
    g2(&g2_,1,param);
    g2_=fmin(g2_,1e100);
    if (verbose >=3)
      Rcout << "u = " << u
            << ", g2(u) = " << g2_ << std::endl;
      return g2_ - value;
  }
};

class Int_g2 {
private:
  ga1_param *param;
  double abs_tol;
  double rel_tol;
  int subdivisions;
  int verbose;
  static std::string msgs[];

public:
  Int_g2(ga1_param* param_in, double tol_in, int subdivisions_in, int verbose_in){
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
                << "Rdqags(g2" << "with param.u0 = " << param->u0
                << ", param.p2b = " << param->p2b
                << ", param.ea = " << param->ea
                << "lower = " << lower
                << ", upper = " << upper << ", " <<std::endl
                << "abs_tol " << abs_tol
                << ", rel_tol " << rel_tol << std::endl;

      Rdqags(g2, (void *) param, &lower, &upper, &abs_tol, &rel_tol,
             &result, &abserr, &neval, &ier,
             &subdivisions, &lenw, &last, iwork, work);
             if (ier>0) {
               if (ier<6){
                 double endpt[] = {lower, upper};
                 g2(endpt,2,(void *) param);
                 double r_min = std::min(endpt[0],endpt[1])*(upper-lower);
                 double r_max = std::max(endpt[0],endpt[1])*(upper-lower);
                 result=std::min(std::max(result,r_min),r_max);
                 warning(msgs[ier]);
               } else {
                 stop(msgs[ier]);
               }
             }
             if (verbose>=3)
               Rcout << "Integral of g2 from theta = " << lower
                         << " to theta = " << upper
                         << " = " << result
                         << " with absolute error - " << abserr
                         << " subintervals = " << last << std::endl;

               return result;
  };
};

std::string Int_g2::msgs[]={"OK","Maximum subdivisions reached","Roundoff error detected",
                            "Bad integrand behavior","Roundoff error in the extrapolation table",
                            "Integral probably divergent","Input is invalid"};

//' Auxiliary for dstable()  only used when alpha == 1 :
//' @param x  numeric *scalar*, >= 0
//' @param beta  0 < |beta| <= 1
//' @param tol
//' @param subdivisions
double fct2(double x, double beta, bool log_flag,
                   double dbltol, int subdivisions, int verbose)
{
// the standard boost eps_tolerance is relative, which doesn't work with
// near infinite numbers
  boost::uintmax_t max_iter;
  abs_eps_tolerance tol(dbltol);

  double i2b = 1/(2*beta);
  double p2b = M_PI*i2b; // = pi/(2 beta)
  double ea = -p2b*x;
  if(!isfinite(ea)) {
    Rcout << "ea is infinite returning 0" << std::endl;
    return (log_flag? R_NegInf : 0);
  }

//' g() is strictly monotone;
//' g(u) := original_g(u*pi/2)
//'  for beta > 0: increasing from g(-1) = 0   to  g(+1) = Inf
//'  for beta < 0: decreasing from g(-1) = Inf to  g(+1) = 0
//t0 = -sign(beta)*pi2# g(t0) == 0  mathematically, but not always numerically
  double u0 = -((beta>=0)-(beta<=0)); // g(u0) == 0  mathematically, but not always numerically
  ga1_param param;
  param.u0=u0;
  param.p2b=p2b;
  param.ea=ea;
  double g2_m1= -1;
  g2(&g2_m1,1,&param);
  double g2_p1= 1;
  g2(&g2_p1,1,&param);
  if (verbose >=2)
    Rcout << "ga1(-1) = " << ga1(-1,&param) << ", g2(-1) = " << g2_m1 << std::endl
          << "ga1(+1) = " << ga1(+1,&param) << ", g2(+1) = " << g2_p1 << std::endl;
  // We know that the maximum of g2(.) is = exp(-1) = 0.3679  "at" g(.) == 1
  // find that using toms748 :
  ga1_solve ga1_s(1,&param,FALSE,verbose);
  max_iter = 200;
  double lower=-1;
  double upper = 1;
  if (verbose >= 2)
    Rcout << "Calling toms748_solve(ga1_s, -1, 1, tol, max_iter = " << max_iter << std::endl;
  std::pair<double,double> ur_pair =  boost::math::tools::toms748_solve(ga1_s, lower, upper, tol, max_iter);
  double u2 = (ur_pair.first+ur_pair.second)/2;
  double g2_u2 = u2;
  g2(&g2_u2,1,&param);
  if (verbose >=2)
    Rcout << "ur_pair.first = " << ur_pair.first
              << ", ur_pair.second = " << ur_pair.second << std::endl
              << "u2 = " << u2 << ", g2(u2) = " << g2_u2 << std::endl;

  double eps = 1e-10;
  g2_solve g2_s(eps,&param,verbose);
  bool do1;
  double u1;
  if (verbose>=2)
    Rcout << std::endl << "u2 = " << u2
          << ", g2(u2) = " << g2_u2 << std::endl;
  if ((do1 = g2_u2 > eps && g2_m1 < eps)) {
    std::pair<double,double> u1_pair;
    max_iter=200;
    double lower = -1;
    double upper = u2;
    u1_pair = boost::math::tools::toms748_solve(g2_s, lower, upper, tol, max_iter);
    u1=(u1_pair.first+u1_pair.second)/2;
    if (verbose>=2){
      double g2_u1 = u1;
      g2(&g2_u1,1,&param);
      Rcout << "u1 = " << u1
            << ", g2(u1) = " << g2_u1 << std::endl;
      Rcout << "Region 1 is "  << M_PI_2*(u1+1) << " long. " << std::endl
            << "Region 2 is " << M_PI_2*(u2-u1) << " long. " << std::endl;
  }
  }
  bool do4;
  double u3;
  if ((do4 = g2_u2 > eps && g2_p1 < eps)) {
    std::pair<double,double> u3_pair;
    max_iter=200;
    double lower = u2;
    double upper = 1;
    u3_pair = boost::math::tools::toms748_solve(g2_s, lower, upper, tol, max_iter);
    u3=(u3_pair.first+u3_pair.second)/2;
    if (verbose>=2){
      double g2_u3 = u3;
      g2(&g2_u3,1,&param);
      Rcout << "u3 = " << u3
            << ", g2(u3) = " << g2_u3 << std::endl;
      Rcout << "Region 3 is " << M_PI_2*(u3-u2) << " long. " << std::endl
            << "Region 4 is " << M_PI_2*(1-u3) << " long. " << std::endl;
    }
  }
  Int_g2 int_g2(&param, dbltol, subdivisions, verbose);
  double r1, r2, r3, r4;
  double rerr=0;
  if (do1) {
    r1=int_g2(-1,u1);
    rerr+=int_g2.abserr;
    r2=int_g2(u1,u2);
    rerr+=int_g2.abserr;
  }
  else {
    r1=0;
    r2=int_g2(-1,u2);
    rerr+=int_g2.abserr;
  }
  if (do4){
    r3=int_g2(u2,u3);
    rerr+=int_g2.abserr;
    r4=int_g2(u3,1);
    rerr+=int_g2.abserr;
  }
  else{
    r3=int_g2(u2,1);
    rerr+=int_g2.abserr;
    r4=0;
  }
  if(verbose) {
    double cc = M_PI_2*fabs(i2b);
    Rprintf(".fct2(%g, %g,..): c*sum(r1+r2+r3+r4)= %g*(%g + %g + %g + %g)= %g\n with abs.err = %g\n",
                    x,beta, cc, r1, r2, r3, r4, cc*(r1+r2+r3+r4), cc*rerr);
  }
  if(log_flag)
    return log(M_PI_2) + log(fabs(i2b)) + log(r1 + r2 + r3 + r4);
  else
    return M_PI_2*fabs(i2b)*(r1 + r2 + r3 + r4);
  } // {.fct2}

double C_stable_tail(double alpha, bool log_flag) {
  if (!(0 <= alpha && alpha <= 2))
    stop("alpha is not between 0 and 2 inclusive\n");
  double r = alpha;
  if (alpha == 0)
    r = (log_flag) ? -log(2) : 0.5;
  else if (alpha == 2)
    r = (log_flag) ? R_NegInf : 0;
  else
    r = (log_flag) ? (lgamma(alpha) - log(M_PI) + log(sin(alpha * M_PI_2)))
              : tgamma(alpha)/M_PI * sin(alpha * M_PI_2);
  return r;
}

double dPareto(double x, double alpha, double beta, bool log_flag) {
  if (x < 0) {
    x = -x;
    beta = -beta;
  }
  if (log_flag)
    return log(alpha) + log1p(beta) + C_stable_tail(alpha, TRUE) -
      (1 + alpha) * log(x);
  else
    return alpha * (1 + beta) * C_stable_tail(alpha,log_flag) * pow(x,-(1 +
                  alpha));
}

double sdstable1(double x, double alpha, double beta, int log_flag,
                         double tol, double zeta_tol, int subdivisions, int verbose)
  {
    Rcout.precision(20);
    Rcout.setf(std::ios::scientific, std::ios::floatfield);
    double ret;
  if (alpha != 1) { // 0 < alpha < 2	&  |beta| <= 1 from above
// General Case
    double tanpa2 = tan( M_PI_2*alpha);
    double betan = beta * tanpa2;
    double zeta = -betan;
    double theta0 = std::min(std::max(-M_PI_2, atan(betan)/alpha), M_PI_2);
    if (verbose)
      Rcout << std::endl << "dstable(., alpha=" << alpha <<  ", beta=" << beta
                << ",..): --> theta0 = " << theta0
                <<", zeta = " << zeta << std::endl;
    if (zeta_tol == 0) {
        if (betan == 0)
          zeta_tol = 4e-16;
        else if (1 - fabs(beta) < 0.01 || alpha < 0.01)
          zeta_tol = 2e-15;
        else
          zeta_tol = 5e-05;
        if (verbose)
          Rcout << " --> zeta.tol = " << zeta_tol << std::endl;
      }
      ret = fct1(x, zeta, alpha, beta, theta0,
                 log_flag, tol, subdivisions,
                 zeta_tol, verbose);
  } // alpha != 1
  else {
    if (verbose)
      Rcout << std::endl << "dstable(., alpha=" << alpha <<  ", beta=" << beta
            << ",..):" << std::endl;
    if (x >= 0) {
      ret = fct2(x, beta, log_flag, tol,
           subdivisions, verbose);
    }
    else {
      ret = fct2(-x, -beta, log_flag , tol,
           subdivisions, verbose);}

  } //alpha == 1
  if (ret == (log_flag ? R_NegInf : 0)){
     ret = dPareto(x, alpha, beta, log_flag);
     if (verbose)
       Rcout<< "dstable1(x = " << x << "alpha = " << alpha
            << ", beta = " << beta << ", log = " << log_flag << ")" << std::endl
            << "fct1 or fct2 returned 0, using dPareto = " << ret << std::endl;
  }
  return ret;

} //sdstaple1
