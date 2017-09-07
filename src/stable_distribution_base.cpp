/// @file stable_distribution_pdf.cpp
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
using std::ostringstream;
using std::pair;
using std::max;
using std::min;
using boost::math::tools::toms748_solve;


/* This is an ancillary functions used by pdf for small alpha */
myFloat pdf_smallA(myFloat x, myFloat alpha, myFloat beta, bool log_flag=false);

bool StandardStableDistribution::initialized = false;
myFloat StandardStableDistribution::pi;
myFloat StandardStableDistribution::pi2;
myFloat StandardStableDistribution::PosInf;
myFloat StandardStableDistribution::NegInf;
myFloat StandardStableDistribution::zeta_tol;
  
void StandardStableDistribution::initialize(){
  if (!Machine::initialized)
    Machine::initialize();
  pi = Machine::pi;
  pi2 = Machine::pi2;
  PosInf = Machine::PosInf;
  NegInf = Machine::NegInf;
  zeta_tol = 200*Machine::epsilon;
  initialized = true;
}

myFloat StandardStableDistribution::g_l(myFloat th_l) const{
  // Similar to g except the input variable is th_l = th+theta0
  myFloat catt_m_t = max<myFloat>(static_cast<myFloat>(0),-sin((alpha-1)*th_l+add_l));
  if ((alpha < 1 && th_l==0) || (alpha > 1 && th_l==th_max)) {
    myFloat g0 = 0;
    if (fabs(beta)==1)
      g0 = pow(cat0*pow(x_m_zet/alpha,alpha),(1/(alpha-1)))*fabs(1-alpha);
    return g0;
  }
  else if ((alpha < 1 && th_l==th_max) || (alpha >1 && th_l==0))
    return PosInf;
  else {
    myFloat costh = max<myFloat>(static_cast<myFloat>(0),sin(th_l-add_l));
    myFloat att = alpha*th_l;
    myFloat pow2;
    
    if (fabs(zeta) < 1 || fabs(x_m_zeta_input+zeta) > .1 * fabs(zeta)) {
      myFloat x_sin = x_m_zet/sin(att);
      myFloat pow1 = pow(x_sin,alpha);
      pow2 = pow(cat0*costh*pow1,(1/(alpha-1)));
    } else {
      myFloat ln_pow1 = alpha * (log(fabs(zeta))+log1p(-x_m_zeta_input/zeta-1)-log(sin(att)));
      myFloat ln_pow2 = (log(cat0) + log(costh)+ln_pow1)/(alpha-1);
      pow2 = exp(ln_pow2);
    }
    return pow2*catt_m_t;
  }

}

myFloat StandardStableDistribution::g_r(myFloat th_r) const{
  // Similar to g except the input variable is th_r = pi/2 - th
  if ((alpha>1 && th_r==th_max) || (alpha<1 && th_r==0) )
    return PosInf;
  else if ((alpha>1 && th_r==0) || (alpha <1 && th_r==th_max)){
    myFloat g0 = 0;
    if (fabs(beta)==1 && add_r==0)
      g0=pow(cat0*pow(x_m_zet/alpha,alpha),(1/(alpha-1)))*fabs(1-alpha);
    return g0;
  }
  else {
    myFloat att = alpha*th_r+add_r;
    myFloat catt_m_t = max(static_cast<myFloat>(0.),static_cast<myFloat>(sin((alpha-1)*th_r+add_r)));
    myFloat pow2;
    if (x_m_zet < 1e100) {
      myFloat costh = sin(th_r);
      myFloat pow1 = pow(x_m_zet/sin(att),alpha);
      pow2 = pow(cat0 * costh * pow1,1/(alpha-1));
    } else {
      myFloat ln_costh = log(sin(th_r));
      myFloat ln_pow1 = alpha*(log(x_m_zet) - log(sin(att)));
      pow2 = exp((log(cat0) + ln_costh + ln_pow1)/(alpha - 1));
    }
    return pow2*catt_m_t;
  }
}

myFloat StandardStableDistribution::ga1_r(myFloat u_r) const{

  myFloat h = p2b+pi2-u_r*pi2;
  myFloat h2b = h/p2b;
  myFloat tanth;
  myFloat h_tanth;
  if(u_r==2){
    tanth = NegInf;
    if (beta == 1)
      h_tanth = -1;
    else
      h_tanth = NegInf;
  } else if (u_r==0) {
    tanth = PosInf;
    if (beta == -1)
      h_tanth = -1;
    else
      h_tanth = PosInf;
  } else {
    tanth = 1/tan(pi*u_r/2);
    h_tanth = h*tanth;
  }
  myFloat exp_ea_p_h_tan_th = exp(ea+h_tanth);
  myFloat costh = sin(pi2*u_r);
  if (u_r==2) {
    if (beta>0) {
      if (beta==1)
        return exp_ea_p_h_tan_th/pi2;
      else
        return 0;
    } else {
      return PosInf;
    }
  } else if (u_r==0){
    if (beta<0){
      if (beta==-1){
        return exp_ea_p_h_tan_th/pi2;
      } else
        return 0;
    } else
      return PosInf;
  } else if (exp_ea_p_h_tan_th ==0)
    return 0;
  else {
    return h2b*exp_ea_p_h_tan_th/costh;
  }
}

myFloat StandardStableDistribution::ga1_c(myFloat u_c) const{
  
  if (fabs(u_c) > .01) {
    // We don't need the high resolution away for u_c==0
    return ga1_r(u_c-th_min);
  }
  myFloat h = h_u2-u_c*pi2;
  myFloat del_ln_h2b = log1p(-u_c/(h_u2/pi2));
  myFloat sin_u_c = sin(u_c * pi2);
  myFloat cos_u_c = cos(u_c * pi2);
  myFloat del_h_tanth;
  
  myFloat tan_u_c;
  myFloat ret;
  if(u_c==th_max){
    if (beta == 1)
      del_h_tanth = -1-h_u2*tanth_u2;
    else
      del_h_tanth = NegInf;
  } else if (u_c == th_min) {
    if (beta == -1)
      del_h_tanth = -1-h_u2*tanth_u2;
    else
      del_h_tanth = PosInf;
  } else {
    tan_u_c = sin_u_c/cos_u_c;
    if (costh_u2*costh_u2 !=0)
      del_h_tanth = -h*tan_u_c/pow(costh_u2,2)/(1+tanth_u2*tan_u_c)-pi2*u_c*tanth_u2;
    else {
      if (h*tan_u_c == 0)
        del_h_tanth = 0;
      else if (tan_u_c < 1)
        del_h_tanth = -h*tan_u_c/pow(costh_u2,2)/(1+tanth_u2*tan_u_c)-pi2*u_c*tanth_u2;
      else
        del_h_tanth = -h/pow(costh_u2,2)/(1/tan_u_c+tanth_u2)-pi2*u_c*tanth_u2;
    }
  }
  myFloat del_ln_costh = log1p(tanth_u2*sin_u_c-2*pow(sin(pi*u_c/4),2));
  if (u_c==th_max) {
    if (beta>0) {
      if (beta==1)
        ret = exp(ln_g_u2-log(h_u2/p2b)+log(costh_u2)+del_h_tanth);
      else
        ret = 0;
    } else {
      ret = PosInf;
    }
  } else if (u_c==th_min){
    if (beta<0){
      if (beta==-1){
        ret = exp(ln_g_u2-log(h_u2/p2b)+log(costh_u2)+del_h_tanth);
      } else
        ret = 0;
    } else
      ret = PosInf;
  }
  else {
    ret = exp(ln_g_u2 + del_ln_h2b + del_h_tanth - del_ln_costh);
  }
  if (isnan(ret)) {
    cerr << "ga1_c(alpha = " << alpha << ", beta = " << beta_input
         << ", x = " << x_m_zeta_input+zeta << ")" << endl
         << "u_c = " << u_c << endl;
    throw(std::range_error("ga1_c: Returning nan"));
  }
  return ret;
}

myFloat StandardStableDistribution::dlng_dth_r(myFloat th_r){
  myFloat t1 = + cos(th_r)/(sin(th_r)*(alpha-1));
  myFloat t2 = - (alpha*alpha)/(alpha-1)*cos(alpha*th_r+add_r)/sin(alpha*th_r+add_r);
  myFloat t3 = + (alpha-1) * cos((alpha-1)*th_r+add_r)/sin((alpha-1)*th_r+add_r);
  return t1+t2+t3;
}

myFloat StandardStableDistribution::dlng_dth_l(myFloat th_l) {
  myFloat t1 = cos(th_l-add_l)/sin(th_l-add_l);
  myFloat t2 = -(alpha*alpha)/(alpha-1)*cos(alpha*th_l)/sin(alpha*th_l);
  myFloat t3 = (1-alpha) * cos((1-alpha)*th_l - add_l)/sin((1-alpha)*th_l - add_l);
  return t1+t2+t3;
}

myFloat StandardStableDistribution::dlnga1_du_r(myFloat u_r) {
  myFloat t1 = -1/(1/beta + 1 - u_r);
  myFloat t2 = - pi * cos(pi2 * u_r)/sin(pi2 * u_r);
  myFloat  t3 = - pi2*pi2* (1/beta + 1 - u_r) * pow(sin(pi2 * u_r),-2);
  return t1+t2+t3;
}

ostream& operator<<(ostream& os, const StandardStableDistribution& dist) {
  os << "StandardStableDistribution: " << endl
     << " alpha = " << dist.alpha << endl
     << " beta_input = " << dist.beta_input << endl
     << " x_input = " << dist.x_input << endl
     << " x_m_zeta_input = " << dist.x_m_zeta_input << endl
     << " beta = " << dist.beta << endl;
  switch (dist.dist_type) {
    case StandardStableDistribution::Cauchy :
      os << "Cauchy distribution." << endl;
      return os;
    case StandardStableDistribution::normal :
      os << "Normal distribution with std. dev. = sqrt(2)." << endl;
      return os;
    case StandardStableDistribution::fin_support :
      os << "Outside support of distribution." << endl;
      return os;
    case StandardStableDistribution::other :
      if (dist.alpha!=1) {
        os << " zeta = " << dist.zeta << endl
           << " theta0_x_gt_zeta = " << dist.theta0_x_gt_zeta << endl
           << " cos(alpha*theta0) = " << dist.cat0 << endl
           << " theta0 = " << dist.theta0 << endl
           << " x_m_zet = " << dist.x_m_zet << endl
           << " add_l = " << dist.add_l << endl
           << " add_r = " << dist.add_r << endl;
      }
      os << " th_max = " << dist.th_max << endl         //Both th_l and th_r range from 0 to th_max
         << " c2 = " << dist.c2 << endl
         << " c_ddx = " << dist.c_ddx << endl
         << " fun_type = " << dist.fun_type << endl;
      
      if (dist.alpha==1) {  // variable used when alpha = 1
        os << " abs_x = " << dist.abs_x << endl
           << " i2b = " << dist.i2b << endl
           << " p2b = " << dist.p2b << endl
        << " ea = " << dist.ea << endl;
        if (dist.fun_type==StandardStableDistribution::fun_ga1_c) {
          os << " ln_g_u2 - " << dist.ln_g_u2 << endl
             << " costh_u2 = " << dist.costh_u2 << endl
             << " tanth_u2 = " << dist.tanth_u2 << endl
             << " h_u2 = " << dist.h_u2 << endl;
        }
      }
      os << " good_theta2 = " << dist.good_theta2 << endl
         << " g(theta2) error = " <<  dist.g_theta2_error << endl
         << " g_dd_theta2 = " << dist.g_dd_theta2 << endl;

      os << "g_map: " << endl << setw(33) << right << "theta" << setw(33) << "g(theta)" <<endl;
      for (vector<myFloat>::const_iterator ppoint=dist.points.begin(); ppoint<dist.points.end(); ppoint++) {
        os << setw(33) << setprecision(24) << scientific << *ppoint
           << setw(33) << setprecision(24) << scientific << dist.g(*ppoint) << ((*ppoint==dist.theta2)? " *":"") << endl;
      }
      return os;
  }
}

void f_of_g (myFloat *th, int n, void *ext) {
  Integral_f_of_g * int_f_of_g = (Integral_f_of_g *) ext;
  for (int i=0; i<n; i++) {
    th[i]= int_f_of_g->f_of_g(th[i]);
  }
}

myFloat Integral_f_of_g::operator() () {
  
  int verbose = controller->get_verbose();
  if (verbose==4){
    cout << endl
         << "IntegrationController::integrate(f_of_g,..)" << endl;
  }
  controller->integrate(stable_distribution::f_of_g, (void *) this, std_stable_dist->points,
                        result, abserr, neval, termination_code, last);
  
  if (verbose>=3){
    myFloat rsum=0, esum=0;
    for (int i=0; i<last; i++) {
      rsum += controller->subs.at(i).r;
      esum += controller->subs.at(i).e;
    }
    
    if (termination_code > 0)
      cout << msgs[termination_code] << ":" << endl;
    cout << "Integral of f_of_g from theta = " << std_stable_dist->points.front()
         << " to theta = " << std_stable_dist->points.back()
         << " = " << result
         << ", with absolute error = " << abserr
         << ", subintervals = " << last << endl
         << "rsum = " << rsum << ", esum = " << esum << endl;
    print_subs_summary(cout, controller->subs, last, std_stable_dist->points);
  }
  if (verbose>=4){
    print_subs(cout, controller->subs, last, std_stable_dist->points);
  }
  return result;
}

string Integral_f_of_g::msgs[]={"OK","Maximum subdivisions reached","Roundoff error detected",
        "Bad integrand behavior","Roundoff error in the extrapolation table",
        "Integral probably divergent","Input is invalid"};

void StandardStableDistribution::set_x_m_zeta(myFloat x_m_zeta_in) {
  if (x_m_zeta_in != x_m_zeta_input) {
    x_m_zeta_input=x_m_zeta_in;
    if (!isfinite(x_m_zeta_in)) {
      if (verbose)
        cout << "set_x_m_zeta: x_m_zeta is not finite." << endl;
      return;
    }
    if (alpha == 1 && beta_input == 0) {
      if (verbose) {
        cout << "set_x_m_zeta: Cauchy distribution" <<endl;
      }
      dist_type=Cauchy;
      return;
    }
    if (alpha == 2) {
      if (verbose) {
        cout << "set_x_m_zeta: normal distribution with std dev of sqrt(2)" <<endl;
      }
      dist_type=normal;
      return;
    }
    if (alpha < 1 && ((beta_input==1 && x_m_zeta_in<=0) || (beta_input==-1 && x_m_zeta_in>=0))) {
      // Not used but it's nice to have them
      beta=-beta_input;
      theta0 = -theta0_x_gt_zeta;
      
      if (verbose) {
        cout << "set_x_m_zeta: outside of distribution support" << endl;
      }
      dist_type=fin_support;
      return;
    }
    dist_type=other;
    if (alpha!=1){
      x_m_zet=fabs(x_m_zeta_in);
      if (x_m_zeta_in>=0){
        beta=beta_input;
        theta0 =theta0_x_gt_zeta;
      } else {
        beta=-beta_input;
        theta0 = -theta0_x_gt_zeta;
      }
      th_min=0;
      th_max=(pi2)+theta0;
      add_l=theta0-(pi2);
      add_r=max<myFloat>(static_cast<myFloat>(0),pi-alpha*(theta0+(pi2)));
      if (alpha<1) {
        if (beta==1)
          add_l=0;
        else if (beta==-1)
          add_l=pi;
      } else if (alpha>1) {
        if (beta==-1)
          add_r=0;
      }
      c2 = (alpha/(pi * fabs(alpha-1) * x_m_zet));
      c_ddx = c2/((alpha-1)*(x_m_zeta_in));
      if (x_m_zet<max<myFloat>(static_cast<myFloat>(1),fabs(zeta)))
        fun_type=fun_g_l;
      else
        fun_type=fun_g_r;
    } else { // alpha = 1
      abs_x=fabs(x_m_zeta_in);
      if (x_m_zeta_in >= 0) {
        beta = beta_input;
      } else {
        beta = -beta_input;
      }
      i2b=1/(2*beta);
      p2b=pi*i2b;
      ea = -p2b*abs_x;
      th_min=0;
      th_max=2;
      fun_type=fun_ga1_r;
      if (abs_x > 1e6){
        myFloat g_hi = g(th_max);
        myFloat g_lo = g(th_min);
        if (g_hi <=1 || g_lo<=1) {
          gSolve g_s(0., this, true);
          myFloat lower = th_min;
          myFloat upper = th_max;
          boost::uintmax_t max_iter = 1000;
          RelativeComparisonTolerance rel_tol(controller->epsrel);
          
          pair<myFloat,myFloat> ur1_pair = toms748_solve(g_s, lower, upper, rel_tol, max_iter);
          if (fabs(g(ur1_pair.second)-g(ur1_pair.first))<1e6){
            myFloat u2 = (ur1_pair.first+ur1_pair.second)/2;
            ln_g_u2 = g(u2);
            th_min -= u2;
            th_max -= u2;
            costh_u2 = sin(u2 * pi2);
            tanth_u2 = cos(u2 * pi2)/sin(u2 * pi2);
            h_u2 = p2b+pi2-u2*pi2;
            fun_type = fun_ga1_c;
          }
        }
        
      }
    }
    map_g();
  }
}

void StandardStableDistribution::map_g() {
  if (verbose) cout << "map_g:" << endl;
  boost::uintmax_t max_iter;
  RelativeComparisonTolerance rel_tol(controller->epsrel);
  
  //' g() is strictly monotone -- Nolan(1997) ["3. Numerical Considerations"]
  //'     alpha >= 1  <==>  g() is falling, ie. from Inf --> 0;  otherwise growing from 0 to +Inf
  
  points.resize(0);
  points.push_back(th_min);
  points.push_back(th_max);
  if (th_min==th_max) return;
  bool do_hi=true, do_lo=true;
  myFloat g_hi = g(th_max);
  myFloat g_lo = g(th_min);
  neval+=2;
  
  if (verbose)
    cout << "  g_hi = " << g_hi << ", g_lo = " << g_lo << endl;
  gSolve g_s(0., this, true);
  if (g_hi>=g_lo && g_lo>1){
    if (verbose)
      cout << "  Theta2 is at th_min" << endl;
    good_theta2=true;
    theta2=th_min;
    g_theta2_error = 0;
    max_iter=0;
    do_lo=false;
  } else if (g_lo>g_hi && g_hi>1){
    if (verbose)
      cout << "  Theta2 is at th_max" << endl;
    good_theta2=true;
    theta2=th_max;
    g_theta2_error = 0;
    max_iter=0;
    do_hi=false;
  } else {
    if (verbose)
      cout << "  Theta2 is in the interior" << endl;
    myFloat upper =th_max;
    myFloat lower = th_min;
    myFloat g_upper = g_hi;
    myFloat g_lower = g_lo;
    th_guess(1,lower, g_lower, upper, g_upper);
    
    if (verbose)
      cout << "  theta 2 range from th_guess = " << lower << " to " << upper << endl;
    // toms748_solve passes lower, upper and max_iter by reference and changes them.
    max_iter = 1000;
    
    pair<myFloat,myFloat> ur1_pair = toms748_solve(g_s, lower, upper, rel_tol, max_iter);
    theta2 = (ur1_pair.first+ur1_pair.second)/2;
    g_theta2 = g(theta2);
    neval+=max_iter+1;
    if (theta2!=th_min && theta2!=th_max) points.push_back(theta2);
    g_theta2_error = fabs(g(ur1_pair.second)-g(ur1_pair.first));
    myFloat cap_g_theta2_error = .01;
    if ( isnan(g_theta2_error) || g_theta2_error > cap_g_theta2_error ){
      if (verbose) {
        cout << endl << "  theta2 is not good, g_theta_error: "
        << g_theta2_error;
        if (isnan(g_theta2_error))
          cout << endl;
        else
          cout << " > " << cap_g_theta2_error << endl;
      }
      good_theta2=false;
//      sort(points.begin(),points.end());
//      return;
    } else
      good_theta2=true;
  }
  
  g_dd_theta2=dlng_dth(theta2);
  myFloat th;
  if (verbose>=2)
    cout << endl << "    theta2 = " << theta2
    << ", g(theta2) = " << g(theta2) << ", iterations = " << max_iter << endl
    << "    ddx_lng(theta2) = " << g_dd_theta2 <<endl;
  vector<double> ln_g_lo={-.1,-.2,-.3,-.4,-.5,-.75,
    -1,-2,-3,-4,-6,-8,-10,-12,-14,-16,
    -18,-20,-24,-28,-32,-36,-40,-50, -60,-70,-80, -100, -120, -140};
  vector<double> ln_g_hi={1, 2, 3, 4, 5, 6};
  th=theta2;
  if (do_lo) {
    vector<double> *ln_g = (isfinite(g_lo)) ? &ln_g_lo : &ln_g_hi;
    for (int j=0; j<(*ln_g).size() && (!isfinite(g_lo) || (*ln_g)[j] > log(g_lo)); j++) {
      if (verbose >=2) {
        cout << "    target g = " << exp((*ln_g)[j]) << endl;
      }
      g_s.set_value(static_cast<myFloat>((*ln_g)[j]));
      pair<myFloat,myFloat> th_pair;
      max_iter=1000;
      myFloat lower = th_min;
      myFloat upper = th;
      myFloat ln_g_th = log(g(th));
      if (!isfinite(ln_g_th)) break;
      bool bracketed = (isfinite(g_lo)) ?(log(g_lo) <= (*ln_g)[j]) && ((*ln_g)[j] < ln_g_th)
      :(ln_g_th <= (*ln_g)[j]) && ((*ln_g)[j] <= log(g_lo));
      if (!bracketed) continue;
      th_pair = toms748_solve(g_s, lower, upper, rel_tol, max_iter);
      neval+=max_iter;
      if (max_iter==1000) break;
      myFloat th_new=th_pair.second;
      if (th_new==th_min) break;
      if (rel_tol(th,th_new)) continue;
      points.push_back(th_new);
      th=th_new;
      if (verbose>=2){
        cout << "    theta = " << th
        << ", g(theta) = " << g(th) << ", Iterations = " << max_iter << endl;
      }
    }
    if (points.size() == 3) {
      // no intervals were added on the low side.  Add one if there's room
      if ( theta2 > 100 * Machine::min) {
        points.push_back((th_min+theta2)/2);
      }
    }
  }
  long npts_low = points.size();
  th=theta2;
  if (do_hi){
    vector<double> *ln_g = (isfinite(g_hi)) ? &ln_g_lo : &ln_g_hi;
    for (int j=0; j<(*ln_g).size() && (!isfinite(g_hi) || (*ln_g)[j] > log(g_hi)); j++) {
      if (verbose >=2) {
        cout << "    target g = " << exp((*ln_g)[j]) << endl;
      }
      g_s.set_value(static_cast<myFloat>((*ln_g)[j]));
      pair<myFloat,myFloat> th_pair;
      max_iter=1000;
      myFloat lower = th;
      myFloat upper = th_max;
      myFloat ln_g_th = log(g(th));
      if (!isfinite(ln_g_th)) break;
      bool bracketed = (isfinite(g_hi)) ? (log(g_hi) <= (*ln_g)[j]) && ((*ln_g)[j] < ln_g_th)
      : (ln_g_th <= (*ln_g)[j]) && ((*ln_g)[j] <= log(g_hi));
      if (!bracketed) continue;
      th_pair = toms748_solve(g_s, lower, upper, rel_tol, max_iter);
      neval+=max_iter;
      if (max_iter==1000) break;
      myFloat th_new=th_pair.first;
      if (fabs(th_new-th_max)< 200*Machine::epsilon*th_max) {
        myFloat del_th{200};
        while (fabs(th-th_max) >= 2*del_th*Machine::epsilon*th_max) {
          th_new = th_max*(1-del_th*Machine::epsilon);
          points.push_back(th_new);
          del_th = 2*del_th;
        }
        break;
      }
      if (rel_tol(th,th_new)) continue;
      points.push_back(th_new);
      th=th_new;
      if (verbose>=2){
        cout << "    theta = " << th
        << ", g(theta) = " << g(th) << ", Iterations = " << max_iter << endl;
      }
    }
    if (points.size() == npts_low) {
      // no intervals were added on the high side.  Add one if there's room
      if ( (th_max-theta2) > 100 * Machine::epsilon * th_max) {
        points.push_back(theta2 + (th_max - theta2)/2);
      }
    }
  }
  sort(points.begin(),points.end());
  return;
}

void StandardStableDistribution::th_guess(const myFloat &value,
                                  myFloat &lower, myFloat &g_lower,
                                  myFloat &upper, myFloat &g_upper) {
  myFloat theta_guess;
  switch(fun_type){
    case fun_g_l:
      theta_guess= x_m_zet*pow(value,(1-alpha)/alpha)*pow(cat0,1/alpha)*sin(-add_l)/alpha;
      break;
    case fun_g_r:
      theta_guess=sin(add_r)/cat0*pow(value,alpha-1)/pow(x_m_zet,alpha);
      break;
    case fun_ga1_r:
      theta_guess = 2 * (1 + beta) / abs_x;
    case fun_ga1_c:
      theta_guess = 0;
  }
  if (lower < theta_guess && theta_guess < upper) {
    myFloat g_theta_guess = g(theta_guess);
    if (g_theta_guess >= value) {
      if (g_upper==PosInf) {
        upper = theta_guess;
        g_upper = g_theta_guess;
      } else {
        lower = theta_guess;
        g_lower = g_theta_guess;
      }
    } else {
      if (g_upper==PosInf) {
        lower = theta_guess;
        g_lower = g_theta_guess;
      } else {
        upper = theta_guess;
        g_upper = g_theta_guess;
      }
    }
  }
}

// ------------------------------------------------------------------------------

myFloat C_stable_tail(myFloat alpha, bool log_flag) {
    if (!(0 <= alpha && alpha <= 2)) {
        throw std::range_error("C_stable_tail: alpha is not between 0 and 2 inclusive");
    }
  myFloat r = alpha;
  if (alpha == 0)
    r = (log_flag) ? -log(static_cast<myFloat>(2)) : static_cast<myFloat>(0.5);
  else if (alpha == 2)
      r = (log_flag) ? StandardStableDistribution::NegInf : static_cast<myFloat>(0);
  else
      r = (log_flag) ? static_cast<myFloat>(lgamma(alpha) - log(StandardStableDistribution::pi) + log(sin(alpha * StandardStableDistribution::pi2)))
                     : static_cast<myFloat>(tgamma(alpha)/StandardStableDistribution::pi * sin(alpha * StandardStableDistribution::pi2));
  return r;
}

 
} // namespace stable_distribution
