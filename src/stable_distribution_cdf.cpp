/// \file stable_distribution_cdf.cpp
/// \author Joseph Dunn
/// \copyright 2016 Joseph Dunn
/// \copyright Distributed under the terms of the GNU General Public License version 3

#include "stable_distribution.h"
#include <iostream>
#include <iomanip>
#include <boost/math/tools/toms748_solve.hpp>

namespace stable_distribution {
  
using std::string;
using std::vector;
using std::sort;
using std::setw;
using std::right;
using std::setprecision;
using std::pair;
using boost::math::tools::toms748_solve;

myFloat exp_m_x (myFloat x, StandardStableDistribution* std_stable_dist) {
  return (x<Machine::large_exp_arg) ? exp(-x) : static_cast<myFloat>(0);}

myFloat one_m_exp_m_x (myFloat x, StandardStableDistribution* std_stable_dist) {return -expm1(-x);}

myFloat StandardStableDistribution::integrate_cdf() {
  // --- cdf()  for alpha < 2 ---
  // Returns the integral that preserves the most relative accuracy
  // cdf adjusts for whether this is the cdf or the complement of the cdf
  // or the difference between F.zeta and F.
  neval=0;
  if (verbose) cout << "integrate_cdf:" << endl;
  if(!good_theta2) {
    if (verbose)
      cout << "  Can't find theta2.  Proceeding without it." << endl;
  }
  myFloat r;

  // When x_m_zeta is large use the version that gets small as abs(x) becomes large
  // When x_m_zeta is small use the version that gets small as x approaches zeta
  bool use_one_m_exp_m_x = ((alpha<1) && (fun_type==fun_g_r))
                            || ((alpha>1) && (fun_type==fun_g_l))
                            || ((alpha==1) && (beta>0));
  if (verbose)
    cout << "Integrand is " << (use_one_m_exp_m_x ? "1 - exp(-g)" : "exp(-g)") << endl;
  Integral_f_of_g int_g1 = use_one_m_exp_m_x ? Integral_f_of_g(&one_m_exp_m_x, this)
                                             : Integral_f_of_g(&exp_m_x, this);
  r=int_g1();
  myFloat c1;
  if (alpha != 1){
    c1 = 1/pi;
  } else{ // alpha == 1
    c1 = .5;
  }
  c_g_theta2_error = g_theta2_error;
  abserr=c1*int_g1.abserr;
  neval+=int_g1.neval;
  termination_code=int_g1.termination_code;
  last=int_g1.last;
  if (verbose) {
    cout << "  c1*sum(r)= " << c1 << " * " << r
         << " = " << c1  *(r) << endl
         << "  abs.err = " << c1*int_g1.abserr << endl
         << "  msg = " << int_g1.msg() << endl;
  }
  return c1*r;
} // StandardStableDistribution::integrate_cdf

// ------------------------------------------------------------------------------


inline myFloat retValue(myFloat F, int useF, int log_p) {
  return (useF) ? ((log_p) ? ((F==0)? StandardStableDistribution::NegInf : log(F)) : F)
                : ((log_p) ? ((F==1)? StandardStableDistribution::NegInf : log1p(-F)) : 1 - F);
  
}

myFloat StandardStableDistribution::cdf(myFloat x, int lower_tail, int log_p, Parameterization pm) {
  neval = 0;
  myFloat ret;
  switch (pm) {
    case S0:
      x_input=x;
      set_x_m_zeta(x-zeta);
      break;
    case S1:
      x_input = x + zeta;
      set_x_m_zeta(x);
  }
  // default values
  c_g_theta2_error = 0;
  abserr = 0;
  if (verbose)
    cout << "cdf: lower_tail = " << lower_tail << ", log_p = " << log_p << endl << *this;
  if(x==PosInf) {
    ret = retValue(1, lower_tail, log_p);
    if (verbose)
      cout << "  x is PosInf, returning " << ret << endl;
    return ret;
  }
  else if (x==NegInf) {
    ret = retValue(0, lower_tail, log_p);
    if (verbose)
      cout << "  x is NegInf, returing " << ret << endl;
    return ret;
  }
  switch (dist_type) {
    case Cauchy :
    {
      myFloat F;
      F = (fabs(x)<=1) ? static_cast<myFloat>(atan(-fabs(x))/pi + static_cast<myFloat>(.5))
                       : static_cast<myFloat>(atan(1/fabs(x))/pi);
      ret = retValue(F, (x>0)!=lower_tail, log_p);
      if (verbose)
        cout << "  Cauchy distribution, returning " << ret << endl;
      return ret;
    }
    case normal :
    {
      myFloat F;
      F = erfc(fabs(x/2))/2;
      ret = retValue(F, (x>0)!=lower_tail, log_p);
      if (verbose)
        cout << "  Normal distribution, returning " << ret << endl;
      return ret;
    }
    case fin_support :
      if(beta_input == 1 && x_m_zeta_input <=0) {
        ret = retValue(0, lower_tail, log_p);
      } else {
        ret = retValue(1, lower_tail, log_p);
      } 
      if (verbose)
        cout << "  x is outside of support, returning " << ret << endl;
      return ret;
    case other :
      if (verbose)
        cout << "cdf: General case:" << endl;
      if (alpha !=1) {
        if(fun_type==fun_g_l) { // We're close to zeta.
          myFloat F_zeta = (alpha<1 && fabs(beta)==1)
                           ? (lower_tail != (beta_input<1)) ? 0 : 1
                           :(lower_tail) ? static_cast<myFloat>(static_cast<myFloat>(.5) - theta0_x_gt_zeta/pi)
                                         : static_cast<myFloat>(static_cast<myFloat>(.5) + theta0_x_gt_zeta/pi);
          if (verbose)
            cout << "cdf: Using difference from F_zeta, " << F_zeta << endl;
          myFloat F = F_zeta -((lower_tail != (x_m_zeta_input > 0)) ? 1 : -1) * integrate_cdf();
          ret = (log_p) ? log(F) : F;
          if (verbose)
            cout << "cdf returning " << ret << endl;
          return ret;
        } else {
          bool useF = !((x > zeta && lower_tail) || (x < zeta && !lower_tail));
          myFloat F = min(static_cast<myFloat>(1),max(static_cast<myFloat>(0),integrate_cdf()));
          ret = retValue(F, useF, log_p);
          if (verbose)
            cout << "cdf: " << endl << "  Using tail integral, returning " << ret << endl;
          return ret;
        }
      } else { // alpha = 1
        bool useF = (x>=0) != lower_tail;
        myFloat F = integrate_cdf();
        ret = retValue(F,useF,log_p);
        if (verbose)
          cout << "cdf: " << endl << "  Using tail integral, returning " << ret << endl;
        return ret;
      }
  } // switch on dist_type
} // StandardStableDistribution::cdf
  
myFloat pPareto(myFloat x, myFloat alpha, myFloat beta, bool lower_tail, bool log_p) {
  bool neg = x < 0;  // left tail
  beta = neg ? -beta : beta;
  if(log_p) {
    return (lower_tail != neg) ?
    log1p(-(1+beta)* C_stable_tail(alpha, false)* pow(fabs(x),-alpha)) :
    log1p(beta)+ C_stable_tail(alpha, true) - alpha*log(fabs(x));
  } else {
    myFloat iF = (1+beta)* C_stable_tail(alpha, false)* pow(fabs(x),-alpha);
    return (lower_tail != neg) ? 1-iF : iF;
  }
} // pPareto
  
} //namespace stable_distribution
