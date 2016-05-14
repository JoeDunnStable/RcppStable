/* If using sourceCpp run Sys.setenv("PKG_CXXFLAGS"="-I /opt/local/include") */


#include "stable_double.h"
#include <limits>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <boost/math/tools/toms748_solve.hpp>
#include <boost/math/special_functions/sin_pi.hpp>
#include <boost/math/special_functions/cos_pi.hpp>
#include <vector>

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
using boost::math::tools::toms748_solve;
using boost::math::cos_pi;
using boost::math::sin_pi;


/* This is an ancillary functions used by sdstable1 for small alpha */
double dstable_smallA(double x, double alpha, double beta, bool log_flag=false);

double g_double_class::Machine_eps=std::numeric_limits<double>::epsilon();
const double g_double_class::alpha_small_dstable = 1e-17;
double g_double_class::large_exp_arg=log(std::numeric_limits<double>::max());
double g_double_class::pi = M_PI;
double g_double_class::pi2 = M_PI_2;
double g_double_class::PosInf = std::numeric_limits<double>::infinity();
double g_double_class::NegInf = -std::numeric_limits<double>::infinity();
vector<subinterval_double> g_double_class::subs;


double g_double_class::g_l(double th_l) const{
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
      else if (th_l == th_max)
          return PosInf;
       else
          return pow2*catt_m_t;
    } else {
      if (th_l==th_max)
        if (theta0>0)
          return g0;
        else
          return 0;
      else if (th_l==0)
        return PosInf;
      else
        return pow2*catt_m_t;
    }
  } else {//abs(beta) != 1
    if ((alpha < 1 && th_l==0) || (alpha > 1 && th_l==th_max))
      return 0;
    else if ((alpha < 1 && th_l==th_max) || (alpha >1 && th_l==0))
      return PosInf;
    else
      return pow2*catt_m_t;
  }
}

double g_double_class::g_r(double th_r) const{
  // Similar to g except the input variable is th_r = pi/2 - th
  double att = alpha*th_r+add_r;
  double costh = sin(th_r);
  double pow1 = pow(x_m_zet/sin(att),alpha);
  double pow2 = pow(cat0 * costh * pow1,1/(alpha-1));
  double catt_m_t = fmax(0.,sin((alpha-1)*th_r+add_r));
  if (fabs(beta)==1){
    double g0 = pow(cat0*pow(x_m_zet/alpha,alpha),(1/(alpha-1)))*fabs(1-alpha);
    if (alpha<1){
      if (th_r==th_max)
        return g0;
      else if (th_r==0)
        return PosInf;
      else
        return pow2*catt_m_t;
    } else if (th_r==0)
      if (add_r==0)
        return g0;
      else
        return 0;
    else if (th_r==th_max)
      return PosInf;
    else
      return pow2*catt_m_t;
  } else {//abs(beta) != 1
    if ((alpha>1 && th_r==th_max) || (alpha<1 && th_r==0) )
      return PosInf;
    else if ((alpha>1 && th_r==0) || (alpha <1 && th_r==th_max))
      return 0;
    else
      return pow2*catt_m_t;
  }
}

double g_double_class::ga1_r(double u_r) const{

  double h = p2b+pi2-u_r*pi2;
  double h2b = h/p2b;
  double tanth;
  double h_tanth;
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
    tanth = cos_pi(u_r/2)/sin_pi(u_r/2);
    h_tanth = h*tanth;
  }
  double exp_ea_p_h_tan_th = exp(ea+h_tanth);
  double costh = sin(pi2*u_r);
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

double g_double_class::dlng_dth_r(double th_r){
    double t1 = + cos(th_r)/(sin(th_r)*(alpha-1));
    double t2 = - (alpha*alpha)/(alpha-1)*cos(alpha*th_r+add_r)/sin(alpha*th_r+add_r);
    double t3 = + (alpha-1) * cos((alpha-1)*th_r+add_r)/sin((alpha-1)*th_r+add_r);
    return t1+t2+t3;
}

double g_double_class::dlng_dth_l(double th_l) {
    double t1 = cos(th_l-add_l)/sin(th_l-add_l);
    double t2 = -(alpha*alpha)/(alpha-1)*cos(alpha*th_l)/sin(alpha*th_l);
    double t3 = (1-alpha) * cos((1-alpha)*th_l - add_l)/sin((1-alpha)*th_l - add_l);
    return t1+t2+t3;
}

double g_double_class::dlnga1_du_r(double u_r) {
    double t1 = -1/(1/beta + 1 - u_r);
    double t2 = - pi * cos(pi2 * u_r)/sin(pi2 * u_r);
    double  t3 = - pi2*pi2* (1/beta + 1 - u_r) * pow(sin(pi2 * u_r),-2);
    return t1+t2+t3;
}

ostream& operator<<(ostream& os, const g_double_class& gc) {
  os << "g_double_class: " << endl
     << " alpha = " << gc.alpha << endl
     << " beta_input = " << gc.beta_input << endl
     << " x_input = " << gc.x_input << endl
     << " beta = " << gc.beta << endl;
    switch (gc.dist_type) {
        case g_double_class::Cauchy :
            os << "Cauchy distribution." << endl;
            return os;
        case g_double_class::normal :
            os << "Normal distribution with std. dev. = sqrt(2)." << endl;
            return os;
        case g_double_class::fin_support :
            os << "Outside support of distribution." << endl;
            return os;
        case g_double_class::other :
            if (gc.alpha!=1) {
                os << " zeta = " << gc.zeta << endl
                << " theta0_x_gt_zeta = " << gc.theta0_x_gt_zeta << endl
                << " cos(alpha*theta0) = " << gc.cat0 << endl
                << " theta0 = " << gc.theta0 << endl
                << " x_m_zet = " << gc.x_m_zet << endl
                << " add_l = " << gc.add_l << endl
                << " add_r = " << gc.add_r << endl;
            }
            os << " th_max = " << gc.th_max << endl         //Both th_l and th_r range from 0 to th_max
            << " c2 = " << gc.c2 << endl
            << " fun_type = " << gc.fun_type << endl;

            if (gc.alpha==1) {  // variable used when alpha = 1
                os << " abs_x = " << gc.abs_x << endl
                << " i2b = " << gc.i2b << endl
                << " p2b = " << gc.p2b << endl
                << " ea = " << gc.ea << endl
                << " u0 = " << gc.u0 << endl;
            }
//            os << "ctl_final:" << endl
//            << *gc.ctl_final << endl <<endl
//            << "zeta_tol = " << gc.zeta_tol << endl
//           << "verbose = " << gc.verbose << endl;
            os << "g_map: " << endl << setw(33) << right << "theta" << setw(33) << "g(theta)" <<endl;
            for (vector<double>::const_iterator ppoint=gc.points.begin(); ppoint<gc.points.end(); ppoint++) {
                os << setw(33) << setprecision(24) << scientific << *ppoint
                   << setw(33) << setprecision(24) << scientific << gc.g(*ppoint) << ((*ppoint==gc.theta2)? " *":"") << endl;
            }
            return os;
    }
}

//' @title  x*exp(-x)  numerically stable, with correct limit 0 for x --> Inf
//' @param x  numeric
//' @return x*exp(x)
//' @author Martin Maechler
inline double x_exp_m_x(double x, g_double_class* ext) {
  double r;
  if(isnan(x))
    r = NAN;
  else if(x > g_double_class::large_exp_arg) // e.g. x == Inf
    r = 0;
  else
    r = x*exp(-x);
  return r;
}

void f_of_g (double *th, int n, void *ext) {
    Int_f_of_g_double * int_f_of_g = (Int_f_of_g_double *) ext;
    for (int i=0; i<n; i++) {
        th[i]= int_f_of_g->f_of_g(th[i]);
    }
}

double Int_f_of_g_double::operator() () {

    if (subinterval_double::ctl->verbose>=4){
        cout << endl
        << "dqagpe(f_of_g, " << "with param" << *param;
        cout << "ctl:" << endl << *subinterval_double::ctl << endl;
    }
    dqagpe(::f_of_g, (void *) this, param->points,
                  result, abserr, neval, ier,
                  g_double_class::subs, last);

    if (subinterval_double::ctl->verbose>=3){
        double rsum=0, esum=0;
        for (int i=0; i<last; i++) {
            rsum += g_double_class::subs.at(i).r;
            esum += g_double_class::subs.at(i).e;
        }

        if (ier > 0)
            cout << msgs[ier] << ":" << endl;
        cout << "Integral of f_of_g from theta = " << param->points.front()
        << " to theta = " << param->points.back()
        << " = " << result
        << ", with absolute error = " << abserr
        << ", subintervals = " << last << endl
        << "rsum = " << rsum << ", esum = " << esum << endl;
    }
    if (subinterval_double::ctl->verbose>=4){
        print_subs(g_double_class::subs, last, param->points);
    }
    return result;
}

string Int_f_of_g_double::msgs[]={"OK","Maximum subdivisions reached","Roundoff error detected",
        "Bad integrand behavior","Roundoff error in the extrapolation table",
        "Integral probably divergent","Input is invalid"};

double g_double_class::integrate_dstable(bool log_flag) {
  // --- dstable(x, alpha, beta, ..)  for alpha < 2 ---
  // For  x = zeta, have special case formula [Nolan(1997)];
  // need to use it also for x ~= zeta, i.e., x.m.zet := |x - zeta| < delta :
  neval=0;
  if (alpha != 1 && !isfinite(zeta)) {
      throw std::range_error("integrate_dstable: alpha !=1 && !is_finite(zeta)");
  }
  if(!good_theta2) {
    if (verbose)
      cout << "integrate_dstable(alpha =" << alpha <<", beta=" << beta_input << ", x=" << x_input << "): Can't find theta2" << endl;
    return log_flag ? NegInf : 0;  // Causes use of dPareto or f.zeta
  }
  double r;
  if (verbose) cout << "integrate_dstable" << endl;
  Int_f_of_g_double int_g1(&x_exp_m_x, this, ctl_final);
  r=int_g1();
  abserr=c2*int_g1.abserr;
  neval+=int_g1.neval;
  ier=int_g1.ier;
  if (ier > 0 && ier != 2) {
      cout << "integrate_dstable(alpha=" << alpha << ", beta=" << beta_input << ", x =" << x_input <<"): " << int_g1.msg() << endl;
  }
  if (verbose) {
    if (verbose)
      cout << "integrate_dstable(" << x_input << " , "
                <<  zeta << ",..):" << endl
                << "c2*sum(r)= "
                << c2 << " * " << r << endl
                << "= " << c2*(r) << endl
                << "abs.err = " << c2*int_g1.abserr << endl;
  }
  if (log_flag)
    return log(c2) + log(r);
  else ;
    return c2 * (r);
} //integrate_g_exp_m_g

void g_double_class::map_g() {
    boost::uintmax_t max_iter;
    rel_eps_tolerance_double rel_tol(ctl_final->epsrel);

    //' g() is strictly monotone -- Nolan(1997) ["3. Numerical Considerations"]
    //'     alpha >= 1  <==>  g() is falling, ie. from Inf --> 0;  otherwise growing from 0 to +Inf

    points.resize(0);
    points.push_back(0);
    points.push_back(th_max);
    if (th_max==0) return;
    bool do_hi=true, do_lo=true;
    double g_hi = g(th_max);
    double g_lo = g(0);
    neval+=2;

    if (verbose)
        cout << "g_hi = " << g_hi << ", g_lo = " << g_lo << endl;
    g_double_solve g_s(0., this, true);
    if (g_hi>=g_lo && g_lo>1){
        if (verbose)
            cout << "Theta2 is at 0" << endl;
        good_theta2=true;
        theta2=0;
        max_iter=0;
        do_lo=false;
    } else if (g_lo>g_hi && g_hi>1){
        if (verbose)
            cout << "Theta2 is at th_max" << endl;
        good_theta2=true;
        theta2=th_max;
        max_iter=0;
        do_hi=false;
    } else {
        if (verbose)
            cout << "Theta2 is in the interior" << endl;
        double upper =th_max;
        double lower = 0;
        double g_upper = g_hi;
        double g_lower = g_lo;
        th_guess(1,lower, g_lower, upper, g_upper);

        if (verbose)
           cout << "theta 2 range from th_guess = " << lower << " to " << upper << endl;
        // toms748_solve passes lower, upper and max_iter by reference and changes them.
        max_iter = 1000;

        pair<double,double> ur1_pair = toms748_solve(g_s, lower, upper, rel_tol, max_iter);
        theta2 = (ur1_pair.first+ur1_pair.second)/2;
        g_theta2 = g(theta2);
         neval+=max_iter+1;
        if (theta2!=0 && theta2!=th_max) points.push_back(theta2);
        if (max_iter==1000 || fabs(g(theta2)-1)>.01) {
            good_theta2=false;
            sort(points.begin(),points.end());
             return;
        } else
            good_theta2=true;
    }

    g_dd_theta2=good_theta2?dlng_dth(theta2):NAN;
    double th;
    if (verbose>=2)
        cout << endl << "theta2 = " << theta2
        << ", g(theta2) = " << g(theta2) << ", iterations = " << max_iter << endl;
    vector<double> ln_g_lo={-.1,-.2,-.3,-.4,-.5,-.75,
                       -1,-2,-3,-4,-6,-8,-10,-12,-14,-16,
                      -18,-20,-24,-28,-32,-36,-40,-50, -60,-70,-80, -100};
    vector<double> ln_g_hi={1, 2, 3, 4, 5};
    th=theta2;
    if (do_lo) {
        vector<double> *ln_g = (isfinite(g_lo)) ? &ln_g_lo : &ln_g_hi;
        for (int j=0; j<(*ln_g).size() && (!isfinite(g_lo) || (*ln_g)[j] > log(g_lo)); j++) {
            if (verbose >=2) {
                cout << "target g = " << exp((*ln_g)[j]) << endl;
            }
            g_s.set_value((*ln_g)[j]);
            pair<double,double> th_pair;
            max_iter=1000;
            double lower = 0;
            double upper = th;
            double ln_g_th = log(g(th));
            if (!isfinite(ln_g_th)) break;
            if ((ln_g_th-(*ln_g)[j])*(log(g_lo) - (*ln_g)[j]) >=0) continue;
            th_pair = toms748_solve(g_s, lower, upper, rel_tol, max_iter);
            neval+=max_iter;
            if (max_iter==1000) break;
            double th_new=(th_pair.first+th_pair.second)/2;
            if (th_new==0) break;
            if (rel_tol(th,th_new)) continue;
            points.push_back(th_new);
            th=th_new;
            if (verbose>=2){
                cout << "theta = " << th
                << ", g(theta) = " << g(th) << ", Iterations = " << max_iter << endl;
            }
        }
        if (points.size() == 1) {
            // no intervals were added on the low side.  Add one if there's room
            if ( theta2 > 100 * Machine_eps) {
                points.push_back(theta2 - 50*Machine_eps);
            }
        }
    }
    long npts_low = points.size();
    th=theta2;
    if (do_hi){
        vector<double> *ln_g = (isfinite(g_hi)) ? &ln_g_lo : &ln_g_hi;
        for (int j=0; j<(*ln_g).size() && (!isfinite(g_hi) || (*ln_g)[j] > log(g_hi)); j++) {
            if (verbose >=2) {
                cout << "target g = " << exp((*ln_g)[j]) << endl;
            }
            g_s.set_value((*ln_g)[j]);
            pair<double,double> th_pair;
            max_iter=1000;
            double lower = th;
            double upper = th_max;
            double ln_g_th = log(g(th));
            if (!isfinite(ln_g_th)) break;
            if ((ln_g_th-(*ln_g)[j])*(log(g_hi) - (*ln_g)[j]) >=0) continue;
            th_pair = toms748_solve(g_s, lower, upper, rel_tol, max_iter);
            neval+=max_iter;
            if (max_iter==1000) break;
            double th_new=(th_pair.first+th_pair.second)/2;
            if (fabs(th_new-th_max)<100*Machine_eps*(th_new+th_max)) break;
            if (rel_tol(th,th_new)) continue;
            points.push_back(th_new);
            th=th_new;
            if (verbose>=2){
                cout << "theta = " << th
                << ", g(theta) = " << g(th) << ", Iterations = " << max_iter << endl;
            }
        }
        if (points.size() == npts_low) {
            // no intervals were added on the high side.  Add one if there's room
            if ( (th_max-theta2) > 100 * Machine_eps) {
                points.push_back(theta2 + 50*Machine_eps);
            }
        }
    }
    sort(points.begin(),points.end());
     return;
}

void g_double_class::th_guess(const double &value,
                         double &lower, double &g_lower,
                         double &upper, double &g_upper) {
    double theta_guess;
    switch(fun_type){
        case fun_g_l:
            theta_guess= x_m_zet*pow(value,(1-alpha)/alpha)*pow(cat0,1/alpha)*sin(-add_l)/alpha;
            break;
        case fun_g_r:
            theta_guess=sin(add_r)/cat0*pow(value,alpha-1)/pow(x_m_zet,alpha);
            break;
        case fun_ga1_r:
            theta_guess = 2 * (1 + beta) / abs_x;
    }
    if (lower < theta_guess && theta_guess < upper) {
        double g_theta_guess = g(theta_guess);
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

//' dstable() for very small alpha > 0
//' ok only for  x > zeta := - beta * tan(pi2 *alpha)
double dstable_smallA(double x, double alpha, double beta, bool log_flag) {
    double r = log(alpha)+log1p(beta)-(1+log(2*x+g_double_class::pi*alpha*beta));
  if (log_flag)
    return r;
  else
    return exp(r);
}


// ------------------------------------------------------------------------------

double C_stable_tail(double alpha, bool log_flag) {
    if (!(0 <= alpha && alpha <= 2)) {
        throw std::range_error("C_stable_tail: alpha is not between 0 and 2 inclusive");
    }
  double r = alpha;
  if (alpha == 0)
    r = (log_flag) ? -log(2) : 0.5;
  else if (alpha == 2)
      r = (log_flag) ? g_double_class::NegInf : 0;
  else
      r = (log_flag) ? (lgamma(alpha) - log(g_double_class::pi) + log(sin(alpha * g_double_class::pi2)))
      : tgamma(alpha)/g_double_class::pi * sin(alpha * g_double_class::pi2);
  return r;
}

double dPareto(double x, double alpha, double beta, bool log_flag) {
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

double pPareto(double x, double alpha, double beta, bool lower_tail, bool log_p) {
    bool neg = x < 0;  // left tail
    beta = neg ? -beta : beta;
    if(log_p) {
        return (lower_tail != neg) ?
                log1p((1+beta)* C_stable_tail(alpha, false)* pow(fabs(x),-alpha)) :
                log1p(beta)+ C_stable_tail(alpha, true) - alpha*log(fabs(x));
    } else {
        double iF = (1+beta)* C_stable_tail(alpha, false)* pow(fabs(x),-alpha);
        return (lower_tail != neg) ? 1-iF : iF;
    }
}

double g_double_class::sdstable1(double x, int log_flag)
{
  cout.precision(20);
  cout.setf(ios::scientific, ios::floatfield);
  set_x(x);
  if (verbose)
    cout << "sdstable1" << endl << *this;
  if (!isfinite(x)) {
    neval = 0;
    abserr = 0;
    return log_flag ? NegInf: 0;
  }
  switch (dist_type) {
      case Cauchy :
        neval = 0;
        abserr = 0;
        return log_flag ? -log(pi) - log(1 + x*x)
                        : 1 / (pi * (1 + x*x));
      case normal :
          neval = 0;
          abserr = 0;
          return log_flag ? -x*x/4 -log(2) -log(pi)/2
                          : exp(-x*x/4)/(2*sqrt(pi));
      case fin_support :
          neval = 0;
          abserr=0;
          return log_flag ? NegInf : 0;
      case other :
          double ret;
          double f_zeta = 0.0;

          if (alpha != 1) { // 0 < alpha < 2	&  |beta| <= 1 from above
              // General Case
              if (verbose)
                  cout << endl << "dstable(., alpha=" << alpha <<  ", beta=" << beta
                  << ",..): --> theta0 = " << theta0
                  <<", zeta = " << zeta << endl;
              if (zeta_tol == 0) {
                  zeta_tol = 200*Machine_eps;
                  if (verbose)
                      cout << " --> zeta.tol = " << zeta_tol << endl;
              }
              // For  x = zeta, have special case formula [Nolan(1997)];
              // need to use it also for x ~= zeta, i.e., x.m.zet := |x - zeta| < delta :
              if (alpha < 1 && fabs(beta)==1)
                  f_zeta = log_flag ? NegInf : 0;
              else
                  f_zeta = (log_flag?
                            lgamma(1 + 1/alpha) + log(cos(theta0)) - (log(pi) + log1p(pow(zeta,2))/(2 * alpha)):
                            tgamma(1 + 1/alpha) * cos(theta0)/(pi * pow(1 + pow(zeta,2),(1/(2 * alpha)))));
              if (x_m_zet <= zeta_tol * (zeta_tol + fmax(fabs(x), fabs(zeta)))) {
                  if (verbose)
                      cout << "sdstable1(" << x
                      << "~=" << zeta
                      << " Using f.zeta()" << endl;
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
                      cout << "sdstable1(" << x << " , "
                      << zeta
                      << ",..): small alpha=" << alpha
                      << "sdstable = " << ret <<"\n";
                  return ret;
              }

          } // alpha != 1
          ret = integrate_dstable(log_flag);
          if (ret == (log_flag ? NegInf : 0)) {
              if (fun_type>1){
                  ret = dPareto(x_input, alpha, beta_input, log_flag);
                  if (verbose){
                      cout<< "dstable1(x = " << x << "alpha = " << alpha
                      << ", beta = " << beta << ", log = " << log_flag << ")" << endl
                      << "integrate_dstable returned 0, using dPareto = " << ret << endl;
                  }
              } else {
                  ret = f_zeta;
                  if (verbose)
                      cout<< "dstable1(x = " << x << "alpha = " << alpha
                      << ", beta = " << beta << ", log = " << log_flag << ")" << endl
                      << "integrate_dstable returned 0, using f_zeta = " << ret << endl;
              }
          }
          return ret;

    }

} //sdstaple1

