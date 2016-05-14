/* If using sourceCpp run Sys.setenv("PKG_CXXFLAGS"="-I /opt/local/include") */

#include "stable_double.h"
#include <iostream>
#include <iomanip>
#include "dqagp_double.h"
#include <boost/math/tools/toms748_solve.hpp>

using std::string;
using std::vector;
using std::sort;
using std::setw;
using std::right;
using std::setprecision;
using std::pair;
using boost::math::tools::toms748_solve;

double exp_m_x (double x, g_double_class* param) {return exp(-x);}

double one_m_exp_m_x (double x, g_double_class* param) {return -expm1(-x);}

// ------------------------------------------------------------------------------
//' Auxiliary for pstable()
double g_double_class::integrate_pstable() {
    // --- pstable()  for alpha < 2 ---
    // Returns the integral that preserves the most relative accuracy
    // pstable1 adjusts for whether this is the cdf or the complement of the cdf
    // or the difference between F.zeta and F.
    neval=0;
    if (alpha !=1 && !isfinite(zeta)) {
        throw std::range_error("integrate_pstable: alpha != 1 && !is_finite(zeta)");
    }
    if(!good_theta2) {
        if (verbose)
            cerr << "integrate_pstable(alpha =" << alpha << ", beta =" << beta_input << ", x=" << x_input << "): Can't find theta2.  Proceeding without it." << endl;
    }
    double r;
    if (verbose) cout << "integrate_pstable" << endl;
    // When x_m_zeta is large use the version that gets small as abs(x) becomes large
    // When x_m_zeta is small use the version that gets small as x approaches zeta
    bool use_one_m_exp_m_x = ((alpha<1) && (fun_type==fun_g_r))
                          || ((alpha>1) && (fun_type==fun_g_l))
                          || ((alpha==1) && (beta>0));
    Int_f_of_g_double int_g1 = use_one_m_exp_m_x ? Int_f_of_g_double(&one_m_exp_m_x, this, ctl_final)
                                          : Int_f_of_g_double(&exp_m_x, this, ctl_final);
    r=int_g1();
    double c1;
    if (alpha != 1){
        c1 = 1/pi;
    } else{ // alpha == 1
        c1 = .5;
    }
    abserr=c1*int_g1.abserr;
    neval+=int_g1.neval;
    ier=int_g1.ier;
    if (ier > 0 && ier != 2) {
        cout << "integrate_pstable(alpha=" << alpha << ", beta=" << beta_input << ", x =" << x_input <<"): " << int_g1.msg() << endl;
    }
    if (verbose) {
       cout << "integrate_pstable(" << x_input << " , "
            <<  zeta << ",..):" << endl
            << "c1*sum(r)= " << c1 << " * " << r << endl
            << " = " << c1  *(r) << endl
            << "abs.err = " << c1*int_g1.abserr << endl;
    }
    return c1*r;
} // g_double_class::integrate_pstable

// ------------------------------------------------------------------------------


inline double retValue(double F, int useF, int log_p) {
  return (useF) ? ((log_p) ? log(F)    : F)
                    : ((log_p) ? log1p(-F) : 1 - F);

}

double g_double_class::spstable1(double x, int lower_tail, int log_p) {
  set_x(x);
    if (verbose)
        cout << "spstable1" << endl << *this;
  if(x==PosInf)
    return retValue(1, lower_tail, log_p);
  else if (x==NegInf)
    return retValue(0, lower_tail, log_p);
  switch (dist_type) {
    case Cauchy :
      neval = 0;
      abserr = 0;
      return lower_tail ? (log_p ? log(atan(x)/pi + .5) : atan(x)/pi + .5)
                        : (log_p ? log(atan(-x)/pi + .5) : atan(-x)/pi + .5);
    case normal :
      neval = 0;
      abserr = 0;
      return lower_tail ? (log_p ? log(1+erf(x/2)) - log(2)
                                 : (1 + erf(x/2)) / 2)
                        : (log_p ? log(1+erf(-x/2)) - log(2)
                                 : (1 + erf(-x/2))/2);
    case fin_support :
      if(beta_input == 1 && x <= zeta) {
        neval = 0;
        abserr = 0;
        return retValue(0., lower_tail, log_p);
      } else if(beta_input == -1 && x >= zeta) {
        neval = 0;
        abserr = 0;
        return retValue(1., lower_tail, log_p);
      }
    case other :
      if (alpha !=1) {
        // FIXME: also provide "very small alpha" case, as in .fct1()
        if (!isfinite(zeta)) {
          throw std::range_error("spstable1: Zeta is not finite");
        }
        if(fun_type==fun_g_l) { // We're close to zeta.
          double F_zeta = (alpha<1 && fabs(beta)==1)
                              ? (lower_tail != (beta_input<1)) ? 0 : 1
                              :(lower_tail) ? (.5 - theta0_x_gt_zeta/pi) : (.5 + theta0_x_gt_zeta/pi);
          double F = F_zeta -((lower_tail != x>zeta) ? 1 : -1) * integrate_pstable();
          return (log_p) ? log(F) : F;
        } else {
          bool useF = !((x > zeta && lower_tail) || (x < zeta && !lower_tail));
          double F = fmin(1,fmax(0,integrate_pstable()));
          return retValue(F, useF, log_p);
        }
      } else { // alpha = 1
        bool useF = (x>=0) != lower_tail;
        double F = integrate_pstable();
        return retValue(F,useF,log_p);
      }
  } // switch on dist_type
}

