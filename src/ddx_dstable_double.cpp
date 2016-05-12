/* If using sourceCpp run Sys.setenv("PKG_CXXFLAGS"="-I /opt/local/include") */

#include "stable_double.h"
#include <limits>
#include <iostream>
#include <iomanip>
#include <boost/math/tools/toms748_solve.hpp>

using std::string;
using std::array;
using std::vector;
using std::sort;
using std::setw;
using std::right;
using std::setprecision;
using std::ios;
using std::pair;
using boost::math::tools::toms748_solve;

double x_1_m_alpha_x_exp_m_x(double x, g_double_class* param) {
    if (exp(-x) == 0)
        return 0;
    else
        return x * (1 - param->alpha*x) * exp(-x);
}

double g_double_class::integrate_ddx_dstable() {
    // --- dstable(x, alpha, beta, ..)  for alpha < 2 ---
    // For  x = zeta, have special case formula [Nolan(1997)];
    // need to use it also for x ~= zeta, i.e., x.m.zet := |x - zeta| < delta :
    neval=0;
    if (!isfinite(zeta)) {
        cerr << "!is_finite(zeta)" << endl;
        return NAN;
    }
    if(!good_theta2) {
        if (verbose)
          cout << "integrate_ddx_dstable(alpha =" << alpha <<", beta=" << beta_input << ", x=" << x_input << "): Can't find theta2" << endl;
        return 0;  // Causes use of dPareto or f.zeta
    }
    double r;
    if (verbose) cout << "integrate_ddx_dstable" << endl;
    integr_ctl_double ctl_ddx(*ctl_final);
    ctl_ddx.epsabs=ctl_ddx.epsrel;
    Int_f_of_g_double int_g1(x_1_m_alpha_x_exp_m_x,this, &ctl_ddx);
    r=int_g1();
    abserr=c2*int_g1.abserr;
    neval+=int_g1.neval;
    ier=int_g1.ier;
    if (ier > 0 && ier !=2) {
        cout << "integrate_ddx_dstable(alpha=" << alpha << ", beta=" << beta_input << ", x =" << x_input <<"): " << int_g1.msg() << endl;
    }
    if (verbose) {
        if (verbose)
            cout << "integrate_ddx_dstable(" << x_input << " , "
            <<  zeta << ",..):" << endl
            << "c_ddx*sum(r)= "
            << c_ddx << " * " << r << endl
            << "= " << c_ddx*(r) << endl
            << "abs.err = " << c_ddx*int_g1.abserr << endl;
    }
    return c_ddx * (r);
} //integrate_ddx_dstable


double g_double_class::ddx_sdstable1(double x)
{
  cout.precision(20);
  cout.setf(ios::scientific, ios::floatfield);
  set_x(x);
  if (verbose)
    cout << "ddx_sdstable1" << endl << *this;
  if (!isfinite(x)) {
    return 0;
  }
  switch (dist_type) {
    case Cauchy :
      neval = 0;
      abserr = 0;
      return  -2*x / (pi * pow(1 + x*x,2));
    case normal :
      neval = 0;
      abserr = 0;
      return  -x*exp(-x*x/4)/(4*sqrt(pi));
    case fin_support :
      neval = 0;
      abserr = 0;
      return 0;
    case other :
      double ret;
      double dfdx_zeta;
          
      if (alpha != 1) { // 0 < alpha < 2	&  |beta| <= 1 from above
        // General Case
        if (verbose)
          cout << endl << "ddx_dstable(., alpha=" << alpha <<  ", beta=" << beta
               << ",..): --> theta0 = " << theta0
               <<", zeta = " << zeta << endl;
        if (zeta_tol == 0) {
          zeta_tol = 200 * Machine_eps;
          if (verbose)
            cout << " --> zeta.tol = " << zeta_tol << endl;
        }
        // For  x = zeta, have special case formula [Nolan(1997)];
        // need to use it also for x ~= zeta, i.e., x.m.zet := |x - zeta| < delta :
        dfdx_zeta = tgamma(1 + 2/alpha) * sin(2*theta0_x_gt_zeta)
                         /(2*pi * pow(1 + pow(zeta,2),(1/(alpha))));
        if (x_m_zet <= zeta_tol * (zeta_tol + std::max(fabs(x), fabs(zeta)))) {
          if (verbose)
            cout << "ddx_sdstable1(" << x
                 << "~=" << zeta
                 << " Using dfdx_zeta()" << endl;
          return dfdx_zeta;
        }
        // the real check should be about the feasibility of g() below, or its integration
              
      } // alpha != 1
      ret = integrate_ddx_dstable();
      if (ret == 0) {
        if (fun_type>1){
          ret = -(1+alpha)*dPareto(x_input, alpha, beta_input, false)/x;
          if (verbose){
            cout << "ddx_sdstable1(x = " << x << "alpha = " << alpha
                 << ", beta = " << beta << ")" << endl
                 << "integrate_ddx_dstable returned 0, using dPareto = " << ret << endl;
          }
        } else {
          ret = dfdx_zeta;
          if (verbose)
            cout << "ddx_sdstable1(x = " << x << "alpha = " << alpha
                 << ", beta = " << beta << ")" << endl
                 << "integrate_ddx_dstable returned 0, using dfdx_zeta = " << ret << endl;
        }
      }
          
      return ret;
  } // switch on dist_type
} //ddx_sdstaple1

