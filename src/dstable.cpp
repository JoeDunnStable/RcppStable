/* If using sourceCpp run Sys.setenv("PKG_CXXFLAGS"="-I /opt/local/include") */

#include "Stable.h"
#include "dqagp_v2.h"
#include <limits>
#include <array>
#include <vector>
#include <iostream>
#include <iomanip>
#include <boost/math/special_functions/cos_pi.hpp>
#include <boost/math/special_functions/sin_pi.hpp>
#include <boost/math/tools/toms748_solve.hpp>

using std::array;
using std::vector;
using std::ostream;
using std::setw;
using std::setprecision;
using std::right;
using std::ios;
using std::pair;
using boost::math::cos_pi;
using boost::math::sin_pi;
using boost::math::tools::toms748_solve;
using boost::uintmax_t;

/* This is an ancillary functions used by sdstable1 for small alpha */
double dstable_smallA(double x, double alpha, double beta, bool log_flag=false);

static const double Machine_eps = std::numeric_limits<double>::epsilon();
static const double alpha_small_dstable = 1e-17;
static const double large_exp_arg = 708.39641853226407875e0;

static const double PosInf = std::numeric_limits<double>::infinity();
static const double NegInf = -PosInf;

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
      return PosInf;
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

double g_class::ga1_r(double u_r) {

  double h = p2b+M_PI_2-u_r*M_PI_2;
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
  double costh = sin_pi(u_r/2);
  if (u_r==2) {
    if (beta>0) {
      if (beta==1)
        return exp_ea_p_h_tan_th/M_PI_2;
      else
        return 0;
    } else {
      return PosInf;
    }
  } else if (u_r==0){
    if (beta<0){
      if (beta==-1){
        return exp_ea_p_h_tan_th/M_PI_2;
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

ostream& operator<<(ostream& os, const g_class& gc) {
  os << "g_class: " << endl
     << " alpha = " << gc.alpha << endl
     << " beta_input = " << gc.beta_input << endl
     << " x_input = " << gc.x_input << endl
     << " beta = " << gc.beta << endl;
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

template<typename T>
class sort_data {
public:
  T a;
  int index;
  sort_data(T a, int index) : a(a), index(index){}
};

template<typename T>
bool sd_lt(const sort_data<T> &lhs,const sort_data<T>& rhs) {
  return lhs.a < rhs.a;
}

// This class contains everything needed to integrate g_exp_m_g except the endpoints
// of the interval.  The operator() double(a,b) actually performs the integration
// using the R version of quadpack, Rdqa gs, which is exported from the R interface
class Int_g_exp_m_g {
private:
  g_class *param;
  double abs_tol;
  double rel_tol;
  int limit;
  int verbose;
  static string msgs[];

public:
  Int_g_exp_m_g(g_class* param, double tol, int subdivisions, int verbose) :
    param(param), abs_tol(0.), rel_tol(tol), limit(subdivisions), verbose(verbose)
    {}

  double result;
  double abserr;
  int neval;
  int ier;
  double operator() (double lower, double upper, int npts2, double pts_in[]) {
    array<double, 100> points;
    for (int i=0; i<npts2-2; i++)
      points[i]=pts_in[i];
    array<subinterval, LIMIT> subs;
    int last;

    if (verbose>=4){
      Rcout << endl
            << "dqagpe(g_exp_m_g, " << "with param" << *param
            << "lower = " << lower
            << ", upper = " << upper << ", " <<endl;
      for (int i=0; i<npts2-2; i++)
        Rcout << "points[" << i << "] = " << points[i] << endl;
      Rcout << "abs_tol = " << abs_tol
            << ", rel_tol = " << rel_tol << endl;
    }
    dqagpe_v2(g_exp_m_g, (void *) param, lower, upper, npts2, points, abs_tol, rel_tol, limit,
              result, abserr, neval, ier,
              subs, last);

    if (ier>0) {
      if (ier<6){
         warning(msgs[ier]);
      } else {
        stop(msgs[ier]);
      }
    }
    if (verbose>=3){
      double rsum=0, esum=0;
      for (int i=0; i<last; i++) {
        rsum += subs.at(i).r;
        esum += subs.at(i).e;
      }

      if (ier > 0)
        Rcout << msgs[ier] << ":" << endl;
      Rcout << "Integral of g_exp_m_g from theta = " << lower
                << " to theta = " << upper
                << " = " << result
                << ", with absolute error = " << abserr
                << ", subintervals = " << last << endl
                << "rsum = " << rsum << ", esum = " << esum << endl;
    }
    if (verbose>=4){
      vector<sort_data<double> > srt_a;
      for (int i=0; i<last; i++) {
        sort_data<double> elem(subs.at(i).a,i);
        srt_a.push_back(elem);
      }
      sort(srt_a.begin(), srt_a.end(),sd_lt<double>);

      vector<sort_data<double> > srt_eord;
      for (int i=0; i<last; i++) {
        sort_data<double> elem(subs.at(i).e,i);
        srt_eord.push_back(elem);
      }
      sort(srt_eord.begin(), srt_eord.end(),sd_lt<double>);

      vector<sort_data<int> > srt_iord;
      for (int i=0; i<last; i++) {
        sort_data<int> elem(srt_eord.at(i).index,last-1-i);
        srt_iord.push_back(elem);
      }
      sort(srt_iord.begin(), srt_iord.end(),sd_lt<int>);

      Rcout << " "
            << setw(13) << right << "a"
            << setw(13) << right << "b"
            << setw(13) << right << "length"
            << setw(13) << right << "r"
            << setw(13) << right << "average"
            << setw(13) << right << "e"
            << setw(5) << right << "rank"
            << setw(6) << right << "level" << endl;
      for (int i=0; i<last;i++){
        int j = srt_a[i].index;
        bool ispt = i==0;
        for (int ipt = 0; ipt<npts2-2; ipt++) {
            ispt=ispt || subs.at(j).a==points.at(ipt);
            if (ispt) break;
        }
        if (ispt)
          Rcout << "*";
        else
          Rcout << " ";
        Rcout << setw(13) << setprecision(5) << subs.at(j).a
              << setw(13) << setprecision(5) << subs.at(j).b
              << setw(13) << setprecision(5) << subs.at(j).b-subs.at(j).a
              << setw(13) << setprecision(5) << subs.at(j).r
              << setw(13) << setprecision(5) << subs.at(j).r/(subs.at(j).b-subs.at(j).a)
              << setw(13) << setprecision(5) << subs.at(j).e
              << setw(5) << srt_iord.at(j).index+1
              << setw(6)  << subs.at(j).level << endl;
      }
    }
    return result;
  }
};


string Int_g_exp_m_g::msgs[]={"OK","Maximum subdivisions reached","Roundoff error detected",
        "Bad integrand behavior","Roundoff error in the extrapolation table",
        "Integral probably divergent","Input is invalid"};

double g_class::integrate_g_exp_m_g(bool log_flag, double dbltol, int subdivisions,
                    double zeta_tol, int verbose) {
  // --- dstable(x, alpha, beta, ..)  for alpha < 2 ---
  // For  x = zeta, have special case formula [Nolan(1997)];
  // need to use it also for x ~= zeta, i.e., x.m.zet := |x - zeta| < delta :
  if (!isfinite(zeta)) stop("!is_finite(zeta)");
  uintmax_t max_iter;
  rel_eps_tolerance tol(dbltol);

  //' g() is strictly monotone -- Nolan(1997) ["3. Numerical Considerations"]
  //'     alpha >= 1  <==>  g() is falling, ie. from Inf --> 0;  otherwise growing from 0 to +Inf

  bool do_hi=true, do_lo=true;
  double g_hi = g(th_max);
  double g_lo = g(0);
  double points[102];
  int npts2 = 2;

  if (verbose)
    Rcout << "g_hi = " << g_hi << ", g_lo = " << g_lo << endl;
  g_solve g_s(0., this, verbose, true);
  double theta2, g1_th2;
  if (g_hi>=g_lo && g_lo>1){
    if (verbose)
      Rcout << "Theta2 is at 0" << endl;
    theta2=0;
    g1_th2 = theta2;
    g_exp_m_g(&g1_th2,1,this);
    do_lo=false;
  } else if (g_lo>g_hi && g_hi>1){
    if (verbose)
      Rcout << "Theta2 is at th_max" << endl;
    theta2=th_max;
    g1_th2 = theta2;
    g_exp_m_g(&g1_th2,1,this);
    do_hi=false;
  } else {
    if (verbose)
      Rcout << "Theta2 is in the interior" << endl;

    pair<double, double> ur1_pair;
    // toms748_solve passes lower, upper and max_iter by reference and changes them.
    max_iter = 1000;
    double upper=th_max;
    double lower=0;
    ur1_pair = toms748_solve(g_s, lower, upper, tol, max_iter);
    theta2 = (ur1_pair.first+ur1_pair.second)/2;
    if (max_iter==1000 || fabs(g(theta2)-1)>.01) {
      return log_flag ? NegInf : 0;  // Causes use of dPareto or
    }
    g1_th2 = theta2;
    g_exp_m_g(&g1_th2,1,this);
    points[npts2-2]=theta2;
    npts2++;
  }
  double th;
  if (verbose>=2)
    Rcout << endl << "theta2 = " << theta2
          << ", g1(theta2) = " << g1_th2 << ", iterations = " << max_iter << endl;
  double ln_g_lo[]={-2,-4,-6,-8,-10,-12,-14,-16,-18,-20,-22,-24,-26,-28,-30,-35,-40,-50,-60,-80,-100};
  double ln_g_hi[]={1, 2, 3, 4, 5};
  th=theta2;
  if (do_lo) {
    double *ln_g;
    int npts;
    if (isfinite(g_lo)){
      ln_g=ln_g_lo;
      npts=20;
    }else{
      ln_g=ln_g_hi;
      npts=5;
    }
    for (int j=0; j<npts && (!isfinite(g_lo) || ln_g[j] > log(g_lo)); j++) {
      if (verbose >=2) {
        Rcout << "target g = " << exp(ln_g[j]) << endl;
      }
      g_s.set_value(ln_g[j]);
      pair<double,double> th_pair;
      max_iter=1000;
      double lower = 0;
      double upper = th;
      if ((log(g(th))-ln_g[j])*(log(g_lo) - ln_g[j]) >=0) continue;
      th_pair = toms748_solve(g_s, lower, upper, tol, max_iter);
      if (max_iter==1000) break;
      double th_new=(th_pair.first+th_pair.second)/2;
      if (th_new==0) break;
      if (fabs((th-th_new))< 100*Machine_eps*(th+th_new)) continue;
      points[npts2-2]=th=th_new;
      npts2++;
      if (verbose>=2){
        Rcout << "theta = " << th
              << ", g(theta1) = " << g(th) << ", Iterations = " << max_iter << endl;
      }
    }
  }
  th=theta2;
  if (do_hi){
    double *ln_g;
    int npts;
    if (isfinite(g_hi)){
      ln_g=ln_g_lo;
      npts=20;
    }else{
      ln_g=ln_g_hi;
      npts=5;
    }
    for (int j=0; j<npts && (!isfinite(g_hi) || ln_g[j] > log(g_hi)); j++) {
      if (verbose >=2) {
        Rcout << "target g = " << exp(ln_g[j]) << endl;
      }
      g_s.set_value(ln_g[j]);
      pair<double,double> th_pair;
      max_iter=1000;
      double lower = th;
      double upper = th_max;
      if ((log(g(th))-ln_g[j])*(log(g_hi) - ln_g[j]) >=0) continue;
      th_pair = toms748_solve(g_s, lower, upper, tol, max_iter);
      if (max_iter==1000) break;
      double th_new=(th_pair.first+th_pair.second)/2;
      if (fabs(th_new-th_max)<100*Machine_eps*(th_new+th_max)) break;
      if (fabs((th-th_new))< 100*Machine_eps*(th+th_new)) continue;
      points[npts2-2]=th=th_new;
      npts2++;
      if (verbose>=2){
        Rcout << "theta = " << th
              << ", g(theta) = " << g(th) << ", Iterations = " << max_iter << endl;
      }
    }
  }
  double r;
  Int_g_exp_m_g int_g1(this, dbltol, subdivisions, verbose);
  r=int_g1(0, th_max, npts2, points);
  if (verbose) {
    if (verbose)
      Rcout << "integrate_g_exp_m_g(" << x_input << " , "
                <<  zeta << ",..): c2*sum(r)= "
                << c2 << " * " << r
                << "= " << c2*(r)
                << ", abs.err = " << c2*int_g1.abserr << endl;
  }
  if (log_flag)
    return log(c2) + log(r);
  else ;
    return c2 * (r);
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
    r = (log_flag) ? NegInf : 0;
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
    return log(alpha) + log1p(beta) + C_stable_tail(alpha, true) -
      (1 + alpha) * log(x);
  else
    return alpha * (1 + beta) * C_stable_tail(alpha,log_flag) * pow(x,-(1 +
                  alpha));
}

double g_class::sdstable1(double x, int log_flag, double tol, double zeta_tol, int subdivisions, int verbose)
{
  Rcout.precision(20);
  Rcout.setf(ios::scientific, ios::floatfield);
  if (!isfinite(x)) {
    if (verbose)
      Rcout << "dstable1: x is not finite.  Returning 0" << endl;
    return log_flag ? NegInf: 0;
  }
  set_x(x);
  if (verbose)
    Rcout << "sdstable1" << endl << *this;
  double ret;
  double f_zeta;

  if (alpha != 1) { // 0 < alpha < 2	&  |beta| <= 1 from above
    // General Case
    if (verbose)
      Rcout << endl << "dstable(., alpha=" << alpha <<  ", beta=" << beta
                << ",..): --> theta0 = " << theta0
                <<", zeta = " << zeta << endl;
    bool finSupp = (fabs(beta) == 1 && alpha < 1);
    if(finSupp) {
      // has *finite* support    [zeta, Inf)  if beta ==  1
      //                         (-Inf, zeta] if beta == -1
      if((beta_input == 1 && x <= zeta) || (beta_input == -1 && x >= zeta))
        return log_flag ? NegInf : 0;
    }
    if (zeta_tol == 0) {
      zeta_tol = 2e-15;
      if (verbose)
        Rcout << " --> zeta.tol = " << zeta_tol << endl;
    }
    // For  x = zeta, have special case formula [Nolan(1997)];
    // need to use it also for x ~= zeta, i.e., x.m.zet := |x - zeta| < delta :
    f_zeta = (log_flag?
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
        Rcout << "sdstable1(" << x << " , "
              << zeta
              << ",..): small alpha=" << alpha
              << "sdstable = " << ret <<"\n";
      return ret;
    }

  } // alpha != 1
  ret = integrate_g_exp_m_g(log_flag, tol, subdivisions, zeta_tol, verbose);
  if (ret == (log_flag ? NegInf : 0)) {
    if (fun_type>1){
      ret = dPareto(x_input, alpha, beta_input, log_flag);
      if (verbose){
       Rcout<< "dstable1(x = " << x << "alpha = " << alpha
           << ", beta = " << beta << ", log = " << log_flag << ")" << endl
           << "integrate_g_exp_m_g returned 0, using dPareto = " << ret << endl;
       }
    } else {
      ret = f_zeta;
      if (verbose)
        Rcout<< "dstable1(x = " << x << "alpha = " << alpha
             << ", beta = " << beta << ", log = " << log_flag << ")" << endl
             << "integrate_g_exp_m_g returned 0, using f_zeta = " << ret << endl;
    }
  }
  return ret;

} //sdstaple1

