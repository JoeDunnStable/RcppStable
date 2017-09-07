/// \file stable_distribution_pdf.cpp
/// \author Joseph Dunn
/// \copyright 2016 Joseph Dunn
/// \copyright Distributed under the terms of the GNU General Public License version 3

#include "stable_distribution.h"
#include <limits>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <boost/math/tools/toms748_solve.hpp>

namespace stable_distribution {

using std::string;
using std::vector;
using std::sort;
using std::setw;
using std::right;
using std::setprecision;
using std::ios;
using std::scientific;
using std::pair;
using std::max;
using std::min;
using boost::math::tools::toms748_solve;

/// x*exp(-x)  numerically stable, with correct limit 0 for x --> Inf
inline myFloat x_exp_m_x(myFloat x, StandardStableDistribution* ext) {
  myFloat r;
  if(isnan(x))
    r = NAN;
  else if(x > Machine::large_exp_arg) // e.g. x == Inf
    r = 0;
  else
    r = x*exp(-x);
  return r;
}

myFloat dPareto(myFloat x, myFloat alpha, myFloat beta, bool log_flag) {
  if (x < 0) {
    x = -x;
    beta = -beta;
  }
  if (log_flag)
    return log(alpha) + log1p(beta) + C_stable_tail(alpha, true) -
      (1 + alpha) * log(x);
  else
    return alpha * (1 + beta) * C_stable_tail(alpha,log_flag) * pow(x,-(1 +
                  alpha));
}

myFloat StandardStableDistribution::integrate_pdf(bool log_flag) {
  // --- pdf(x, alpha, beta, ..)  for alpha < 2 ---
  myFloat r;
  if (verbose) cout << "integrate_pdf:" << endl
                    << "Integrand is g * exp(-g)" << endl;
  
  Integral_f_of_g int_g1(&x_exp_m_x, this);
  r=int_g1();
  abserr=c2*int_g1.abserr;
  c_g_theta2_error = c2 * g_theta2_error;
  neval+=int_g1.neval;
  termination_code=int_g1.termination_code;
  last=int_g1.last;
  if (verbose) {
    cout << "  c2*sum(r)= " << c2 << " * " << r << " = " << c2*(r) << endl
         << "  abs.err = " << c2*int_g1.abserr << endl
         << "  msg = " << int_g1.msg() << endl;
  }
  if (log_flag)
    return log(c2) + log(r);
  else ;
  return c2 * (r);
} //integrate_pdf

myFloat StandardStableDistribution::pdf(myFloat x, int log_flag, Parameterization pm)
{
  cout.precision(20);
  cout.setf(ios::scientific, ios::floatfield);
  neval = 0;
  switch (pm) {
    case S0:
      x_input=x;
      set_x_m_zeta(x-zeta);
      break;
    case S1:
      x_input = x + zeta;
      set_x_m_zeta(x);
  }
  myFloat ret;
  if (verbose)
    cout << "pdf: log_flag = " << log_flag << endl << *this;
  if (!isfinite(x)) {
    abserr = 0;
    termination_code=IntegrationController::normal;
    return log_flag ? NegInf: 0;
  }
  switch (dist_type) {
    case Cauchy :
      ret =log_flag ? static_cast<myFloat>(-log(pi) - log(1 + x*x))
                    : static_cast<myFloat>(1 / (pi * (1 + x*x)));
      abserr = 0;
      if (verbose)
        cout << "  Using Cauchy Distribution = " << ret << endl;
      return ret;
    case normal :
      ret = log_flag ? static_cast<myFloat>(-x*x/4 -log(static_cast<myFloat>(2)) -log(pi)/2)
                     : exp(-x*x/4)/(2*sqrt(pi));
      abserr = 0;
      if (verbose)
        cout << "  Using Normal Distribution = " << ret << endl;
      return ret;
    case fin_support :
      ret = log_flag ? NegInf : 0;
      abserr=0;
      if (verbose)
        cout << "  Outside of support.  Returning " << ret << endl;
      return ret;
    case other :
      myFloat ret;
      myFloat f_zeta = 0.0;
      
      if (verbose)
        cout << "  General case" << endl;
      if (alpha != 1) { // 0 < alpha < 2	&  |beta| <= 1 from above
        // For  x = zeta, have special case formula [Nolan(1997)];
        if (alpha < 1 && fabs(beta)==1)
          f_zeta = log_flag ? NegInf : 0;
        else
          f_zeta = (log_flag ? static_cast<myFloat>(lgamma(1 + 1/alpha) + log(cos(theta0)) - (log(pi) + log1p(pow(zeta,2))/(2 * alpha)))
                    : static_cast<myFloat>(tgamma(1 + 1/alpha) * cos(theta0)/(pi * pow(1 + pow(zeta,2),(1/(2 * alpha))))));

        // need to use it also for x ~= zeta
        myFloat x_m_zet_tol = 0;
        switch (pm) {
          case S0:
            x_m_zet_tol = zeta_tol * (zeta_tol + std::max(fabs(x), fabs(zeta)));
            break;
          case S1:
            x_m_zet_tol = std::numeric_limits<myFloat>::min();
            break;
        }
        
        if (x_m_zet <= x_m_zet_tol) {
          ret = f_zeta;
          if (verbose)
            cout << "  " << x_input
            << " ~= " << zeta
            << " Using f_zeta() = " << ret << endl;
          return ret;
        }
        
      } // alpha != 1
      if (good_theta2) {
        ret = integrate_pdf(log_flag);
        if (verbose) cout << "pdf:" << endl;
        myFloat error = max(c_g_theta2_error, abserr);
        if (error < controller->epsrel * ret || (fun_type == 1) || (fun_type == 3 && x < 10)) {
          if (verbose)
            cout << "  Error is below threshhold or x is small. Using integral = " << ret << endl;
        } else {
          myFloat tmp = dPareto(x_m_zeta_input+zeta, alpha, beta_input, log_flag);
          if ( true /*error < pow(fabs(x), -2*(1+alpha)) */) {
            if (verbose)
              cout << "  Error is above threshhold but better than dPareto. Using integral = " << ret << endl;
          } else {
            ret = tmp;
            if (verbose)
              cout << "  Error is above threshhold and worse than dPareto.  Using dPareto = " << ret << endl;
          }
        }
      } else { // bad theta2
        abserr = NAN;
        termination_code = IntegrationController::bad_integrand;
        last = 0;
        if (fun_type>1){
          ret = dPareto(x_m_zeta_input+zeta, alpha, beta_input, log_flag);
          if (verbose){
            cout<< "  Theta2 is bad & x is large so using dPareto = " << ret << endl;
          }
        } else {
          ret = f_zeta;
          if (verbose)
            cout<< " Theta2 is bad & x is small, so using f_zeta = " << ret << endl;
        }
      } // bad theta2
      return ret;
  } // switch on dist_type

} // pdf
  
} // namespace stable_distribution
