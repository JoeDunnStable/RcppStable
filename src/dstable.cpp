/* If using sourceCpp run Sys.setenv("PKG_CXXFLAGS"="-I /opt/local/include") */

#include <Rcpp.h>
#include "Stable.h"
#include <Rmath.h>
#include <R_ext/Applic.h> // For integr_fn and Rdqags
#include <limits>
#include <iostream>
#include <boost/math/tools/toms748_solve.hpp>

/* This is an ancillary functions used by sdstable1 for small alpha */
double dstable_smallA(double x, double alpha, double beta, bool log_flag=FALSE);

static const double Machine_eps = std::numeric_limits<double>::epsilon();
static const double alpha_small_dstable = 1e-17;
static const double large_exp_arg = 708.39641853226407875e0;

double g_class::g_l(double th_l) {
  // Similar to g except the input variable is th_l = th+theta0
  double costh = sin(th_l-add_l);
  double att = alpha*th_l;
  double x_sin = x_m_zet/sin(att);
  double pow1 = pow(x_sin,alpha);
  double pow2 = pow(cat0*costh*pow1,(1/(alpha-1)));
  double catt_m_t = -sin((alpha-1)*th_l+add_l);
  if (fabs(beta)==1){
    double g0 = pow(cat0*pow(x_m_zet/alpha,alpha),(1/(alpha-1)))*fabs(1-alpha);
    if (alpha<1){
      if (th_l==0)
        return g0;
      else
        return pow2*catt_m_t;
    } else if (theta0>0 && th_l==th_max)
        return g0;
    else
      return pow2*catt_m_t;
  } else {//abs(beta) != 1
    if ((alpha < 1 && th_l==0) || (alpha > 1 && th_l==th_max))
      return 0;
    else if ((alpha < 1 && th_l==th_max) || (alpha >1 && th_l==0))
      return R_PosInf;
    else
    return pow2*catt_m_t;
  }
}

double g_class::g_r(double th_r) {
  // Similar to g except the input variable is th_r = M_PI_2 - th
  double att = alpha*th_r+add_r;
  double costh = sin(th_r);
  double pow1 = pow(x_m_zet/sin(att),alpha);
  double pow2 = pow(cat0 * costh * pow1,1/(alpha-1));
  double catt_m_t = sin((alpha-1)*th_r+add_r);
  if (fabs(beta)==1){
    double g0 = pow(cat0*pow(x_m_zet/alpha,alpha),(1/(alpha-1)))*fabs(1-alpha);
    if (alpha<1){
      if (th_r==th_max)
        return g0;
      else if (th_r==0)
        return R_PosInf;
      else
        return pow2*catt_m_t;
    } else if (th_r==0)
      if (add_r==0)
        return g0;
      else
        return 0;
    else if (th_r==th_max)
      return R_PosInf;
    else
      return pow2*catt_m_t;
  } else {//abs(beta) != 1
    if ((alpha>1 && th_r==th_max) || (alpha<1 && th_r==0) )
      return R_PosInf;
    else if ((alpha>1 && th_r==0) || (alpha <1 && th_r==th_max))
      return 0;
    else
      return pow2*catt_m_t;
  }
}

double g_class::ga1_r(double u_r) {

  double h = p2b+M_PI_2-u_r*M_PI_2;
  double h2b = h/p2b;
  double tanth;
  double h_tanth;
  if(u_r==2){
    tanth = R_NegInf;
    if (beta == 1)
      h_tanth = -1;
    else
      h_tanth = R_NegInf;
  } else if (u_r==0) {
    tanth = R_PosInf;
    if (beta == -1)
      h_tanth = -1;
    else
      h_tanth = R_PosInf;
  } else {
    tanth = cospi(u_r/2)/sinpi(u_r/2);
    h_tanth = h*tanth;
  }
  double exp_ea_p_h_tan_th = exp(ea+h_tanth);
  double costh = sinpi(u_r/2);
  if (u_r==2) {
    if (beta>0) {
      if (beta==1)
        return exp_ea_p_h_tan_th/M_PI_2;
      else
        return 0;
    } else {
      return R_PosInf;
    }
  } else if (u_r==0){
    if (beta<0){
      if (beta==-1){
        return exp_ea_p_h_tan_th/M_PI_2;
      } else
        return 0;
    } else
      return R_PosInf;
  } else if (exp_ea_p_h_tan_th ==0)
    return 0;
  else {
    return h2b*exp_ea_p_h_tan_th/costh;
  }
}

std::ostream& operator<<(std::ostream& os, const g_class& gc) {
  os << "g_class: " << std::endl
     << " alpha = " << gc.alpha << std::endl
     << " beta_input = " << gc.beta_input << std::endl
     << " x_input = " << gc.x_input << std::endl
     << " beta = " << gc.beta << std::endl;
  if (gc.alpha!=1) {
    os << " zeta = " << gc.zeta << std::endl
       << " theta0_x_gt_zeta = " << gc.theta0_x_gt_zeta << std::endl
       << " cos(alpha*theta0) = " << gc.cat0 << std::endl
       << " theta0 = " << gc.theta0 << std::endl
       << " x_m_zet = " << gc.x_m_zet << std::endl
       << " add_l = " << gc.add_l << std::endl
       << " add_r = " << gc.add_r << std::endl;
  }
  os << " th_max = " << gc.th_max << std::endl         //Both th_l and th_r range from 0 to th_max
     << " c2 = " << gc.c2 << std::endl
     << " fun_type = " << gc.fun_type << std::endl;

  if (gc.alpha==1) {  // variable used when alpha = 1
    os << " abs_x = " << gc.abs_x << std::endl
       << " i2b = " << gc.i2b << std::endl
       << " p2b = " << gc.p2b << std::endl
       << " ea = " << gc.ea << std::endl
       << " u0 = " << gc.u0 << std::endl;
  }
  return os;
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

// Function to integrate: dstable(..)= f(..) = c2 * \int_{-\theta_0}^{\M_PI_2} g1(u) du
void g_exp_m_g(double * th, int n, void * ext) {
  g_class * param = (g_class *) ext;
  int i;

  for (i=0; i<n; i++){
    double g_th = param->g(th[i]);
    // g1 :=  g(.) exp(-g(.))
    th[i] = x_exp_m_x(g_th);
  }
}

// This class contains everything needed to integrate g_exp_m_g except the endpoints
// of the interval.  The operator() double(a,b) actually performs the integration
// using the R version of quadpack, Rdqa gs, which is exported from the R interface
class Int_g_exp_m_g {
private:
  g_class *param;
  double abs_tol;
  double rel_tol;
  int subdivisions;
  int verbose;
  static std::string msgs[];

public:
  Int_g_exp_m_g(g_class* param, double tol, int subdivisions, int verbose) :
    param(param), abs_tol(tol), rel_tol(tol), subdivisions(subdivisions), verbose(verbose)
    {}

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
                << " beta = " << param->beta
                << ", param.cat0 = " << param->cat0
                << ", param.x_m_zet = " << param->x_m_zet << ", " << std::endl
                << "lower = " << lower
                << ", upper = " << upper << ", " <<std::endl
                << "abs_tol " << abs_tol
                << ", rel_tol " << rel_tol << std::endl;

    Rdqags(g_exp_m_g, (void *) param, &lower, &upper, &abs_tol, &rel_tol,
           &result, &abserr, &neval, &ier,
           &subdivisions, &lenw, &last, iwork, work);
    if (ier>0) {
      if (ier<6){
         warning(msgs[ier]);
      } else {
        stop(msgs[ier]);
      }
    }
    if (verbose>=3){
      if (ier > 0)
        Rcout << msgs[ier] << ":" << std::endl;
      Rcout << "Integral of g1 from theta = " << lower
                << " to theta = " << upper
                << " = " << result
                << ", with absolute error = " << abserr
                << ", subintervals = " << last << std::endl;
    }
    return result;
  };
};

std::string Int_g_exp_m_g::msgs[]={"OK","Maximum subdivisions reached","Roundoff error detected",
        "Bad integrand behavior","Roundoff error in the extrapolation table",
        "Integral probably divergent","Input is invalid"};

double g_class::integrate_g_exp_m_g(bool log_flag, double dbltol, int subdivisions,
                    double zeta_tol, int verbose) {
  // --- dstable(x, alpha, beta, ..)  for alpha < 2 ---
  // For  x = zeta, have special case formula [Nolan(1997)];
  // need to use it also for x ~= zeta, i.e., x.m.zet := |x - zeta| < delta :
  if (!isfinite(zeta)) stop("!is_finite(zeta)");
  boost::uintmax_t max_iter;
  rel_eps_tolerance tol(dbltol);

  //' g() is strictly monotone -- Nolan(1997) ["3. Numerical Considerations"]
  //'     alpha >= 1  <==>  g() is falling, ie. from Inf --> 0;  otherwise growing from 0 to +Inf

  bool do1=TRUE, do2=TRUE, do3=TRUE, do4=TRUE;
  double g_hi = g(th_max);
  double g_lo = g(0);
  if (verbose)
    Rcout << "g_hi = " << g_hi << ", g_lo = " << g_lo << std::endl;
  double theta2, g1_th2;
  if (g_hi>=g_lo && g_lo>1){
    if (verbose)
      Rcout << "Theta2 is at 0" << std::endl;
    theta2=0;
    g1_th2 = theta2;
    g_exp_m_g(&g1_th2,1,this);
    do1=FALSE;
    do2=FALSE;
  } else if (g_lo>g_hi && g_hi>1){
    if (verbose)
      Rcout << "Theta2 is at th_max" << std::endl;
    theta2=th_max;
    g1_th2 = theta2;
    g_exp_m_g(&g1_th2,1,this);
    do3=FALSE;
    do4=FALSE;
  } else {
    if (verbose)
      Rcout << "Theta2 is in the interior" << std::endl;

    g_solve g_s(0., this, verbose, TRUE);
    std::pair<double, double> ur1_pair;
    // toms748_solve passes lower, upper and max_iter by reference and changes them.
    max_iter = 1000;
    double upper=th_max;
    double lower=0;
    ur1_pair = boost::math::tools::toms748_solve(g_s, lower, upper, tol, max_iter);

    theta2 = (ur1_pair.first+ur1_pair.second)/2;
    g1_th2 = theta2;
    g_exp_m_g(&g1_th2,1,this);
  }
  double th1, th3;
  if (verbose>=2)
    Rcout << std::endl << "theta2 = " << theta2
          << ", g1(theta2) = " << g1_th2 << ", iterations = " << max_iter << std::endl;
  if (do2) {
    double g1_lo= 0;
    g_exp_m_g(&g1_lo,1,this);
    double eps = 1e-10;
    if ((do1 = g1_th2 > eps && g1_lo < eps)) {
      g_exp_m_g_solve g1_s(eps,this,verbose);
      std::pair<double,double> th1_pair;
      max_iter=1000;
      double lower = 0;
      double upper = theta2;
      th1_pair = boost::math::tools::toms748_solve(g1_s, lower, upper, tol, max_iter);
      th1=(th1_pair.first+th1_pair.second)/2;
      if (verbose>=2){
        double g1_th1 = th1;
        g_exp_m_g(&g1_th1,1,this);
        Rcout << "theta1 = " << th1
              << ", g1(theta1) = " << g1_th1 << ", Iterations = " << max_iter << std::endl;
        Rcout << "Region 1 is " << (th1) << " long. " << std::endl
              << "Region 2 is " << (theta2-th1) << " long. " << std::endl;
      }
    } else {
      if (verbose>=2)
        Rcout << "There is no Region 1" << std::endl
               << "Region 2 is " << theta2 << " long." << std::endl;
    }
  }
  if (do3){
    double g1_hi=th_max;
    g_exp_m_g(&g1_hi,1,this);
    double eps = 1e-10;
    if ((do4 = g1_th2 > eps && g1_hi < eps)) {
      g_exp_m_g_solve g1_s(eps,this,verbose);
      std::pair<double,double> th3_pair;
      max_iter=1000;
      double lower = theta2;
      double upper = th_max;
      th3_pair = boost::math::tools::toms748_solve(g1_s, lower, upper, tol, max_iter);
      th3=(th3_pair.first+th3_pair.second)/2;
      th3=fmin(th3,theta2+.01*(th_max-theta2));  // The point is to get the sharp peak
      if (verbose>=2){
        double g1_th3 = th3;
        g_exp_m_g(&g1_th3,1,this);
        Rcout << "theta3 = " << th3
              << ", g1(theta3) = " << g1_th3 << ", Iterations = " << max_iter << std::endl;
        Rcout << "Region 3 is " << (th3-theta2) << " long. " << std::endl
              << "Region 4 is " << (th_max-th3) << " long. " << std::endl;
      }
    } else {
      if (verbose>=2)
        Rcout << "There is no Region 4" << std::endl
              << "Region 3 is " << th_max-theta2 << " long." << std::endl;
    }
  }
  double r1, r2, r3, r4;
  double rerr=0;
  Int_g_exp_m_g int_g1(this, dbltol*fmin(1,1/c2), subdivisions, verbose);
  if (do1) {
    r1 = int_g1(0, th1);
    rerr=rerr+int_g1.abserr;
  }
  else {
    r1 = 0;
  }
  if (do2) {
    r2 = int_g1(do1?th1:0, theta2);
    rerr=rerr+int_g1.abserr;
  } else {
    r2 = 0;
  }
  if (do3) {
    r3 = int_g1(theta2, do4?th3:th_max);
    rerr=rerr+int_g1.abserr;
  } else {
    r3 = 0;
  }
  if (do4) {
    r4 = int_g1(th3, th_max);
    rerr=rerr+int_g1.abserr;
  }
  else {
    r4 = 0;
  }
  if (verbose) {
    if (verbose)
      Rcout << "integrate_g_exp_m_g(" << x_input << " , "
                <<  zeta << ",..): c2*sum(r[1:4])= "
                << c2 << "*"
                << "(" << r1
                << " + " << r2
                << " + " << r3
                << " + " << r4
                << ") = " << c2*(r1+r2+r3+r4)
              << ", abs.err = " << c2*rerr << std::endl;
  }
/*  if (rerr > .25*fabs((r1+r2+r3+r4)) && )
    return (log_flag) ? R_NegInf : 0;  ///this will cause the use of dPareto.
*/
  if (log_flag)
    return log(c2) + log(r1 + r2 + r3 + r4);
  else ;
    return c2 * (r1 + r2 + r3 + r4);
} //integrate_g_exp_m_g

//' dstable() for very small alpha > 0
//' ok only for  x > zeta := - beta * tan(M_PI_2 *alpha)
double dstable_smallA(double x, double alpha, double beta, bool log_flag) {
  double r = log(alpha)+log1p(beta)-(1+log(2*x+M_PI*alpha*beta));
  if (log_flag)
    return r;
  else
    return exp(r);
}


// ------------------------------------------------------------------------------

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

double g_class::sdstable1(double x, int log_flag, double tol, double zeta_tol, int subdivisions, int verbose)
  {
    Rcout.precision(20);
    Rcout.setf(std::ios::scientific, std::ios::floatfield);
    if (!isfinite(x)) {
      if (verbose)
        Rcout << "dstable1: x is not finite.  Returning 0" << std::endl;
      return log_flag ? R_NegInf: 0;
    }
    set_x(x);
    if (verbose)
      Rcout << "sdstable1" << std::endl << *this;
    double ret;
  if (alpha != 1) { // 0 < alpha < 2	&  |beta| <= 1 from above
    // General Case
    if (verbose)
      Rcout << std::endl << "dstable(., alpha=" << alpha <<  ", beta=" << beta
                << ",..): --> theta0 = " << theta0
                <<", zeta = " << zeta << std::endl;
    bool finSupp = (fabs(beta) == 1 && alpha < 1);
    if(finSupp) {
      // has *finite* support    [zeta, Inf)  if beta ==  1
      //                         (-Inf, zeta] if beta == -1
      if((beta_input == 1 && x <= zeta) || (beta_input == -1 && x >= zeta))
        return log_flag ? R_NegInf : 0;
    }
    if (zeta_tol == 0) {
/*      if (zeta == 0)
        zeta_tol = 4e-16;
      else if (1 - fabs(beta) < 0.01 || alpha < 0.01)
*/        zeta_tol = 2e-15;
/*      else
        zeta_tol = 5e-05;
*/      if (verbose)
        Rcout << " --> zeta.tol = " << zeta_tol << std::endl;
    }
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
        Rcout << "sdstable1(" << x
              << "~=" << zeta
              << " Using f.zeta()" << std::endl;
        return f_zeta;
    }
    // the real check should be about the feasibility of g() below, or its integration
    bool smallAlpha = alpha < alpha_small_dstable;

    if (smallAlpha) {
      if (x<zeta) {
        ret = dstable_smallA(-x, alpha, -beta, log_flag);
      } else {
        ret = dstable_smallA(x, alpha, beta, log_flag);
      }
      if (verbose)
        Rcout << "sdstable1(" << x << " , "
              << zeta
              << ",..): small alpha=" << alpha
              << "sdstable = " << ret <<"\n";
      return ret;
    }

  } // alpha != 1
  ret = integrate_g_exp_m_g(log_flag, tol, subdivisions, zeta_tol, verbose);
  if (ret == (log_flag ? R_NegInf : 0)  && fabs(x)>1000){
    ret = dPareto(x_input, alpha, beta_input, log_flag);
     if (verbose)
       Rcout<< "dstable1(x = " << x << "alpha = " << alpha
            << ", beta = " << beta << ", log = " << log_flag << ")" << std::endl
            << "integrate_g_exp_m_g returned 0, using dPareto = " << ret << std::endl;
  }
  return ret;

} //sdstaple1

// [[Rcpp::export]]
double test_g(double th, double x, double alpha, double beta) {
  g_class param(alpha, beta);
  param.set_x(x);
  Rcout<<param;
  return param.g(th);
}
