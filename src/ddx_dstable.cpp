/* If using sourceCpp run Sys.setenv("PKG_CXXFLAGS"="-I /opt/local/include") */

//#include <Rcpp.h>
#include "Stable.h"
//#include <Rmath.h>
#include "dqagp_v2.h"
#include <limits>
#include <iostream>
#include <iomanip>
#include <boost/math/tools/toms748_solve.hpp>

static const double Machine_eps = std::numeric_limits<double>::epsilon();
static const double alpha_small_dstable = 1e-17;
static const double large_exp_arg = 708.39641853226407875e0;

// Function to integrate: dstable(..)= f(..) = c2 * \int_{-\theta_0}^{\M_PI_2} g1(u) du
void g_1_m_alpha_g_exp_m_g(double * th, int n, void * ext) {
  g_class * param = (g_class *) ext;
  int i;

  for (i=0; i<n; i++){
    double g_th = param->g(th[i]);
    // g1 :=  g(.)(1-alpha g) exp(-g(.))
    if (exp(-g_th)==0)
      th[i]=0;
    else
      th[i] = g_th*(1-param->alpha*g_th)*exp(-g_th);
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

// This class contains everything needed to integrate g_1_m_alpha_g_exp_m_g except the endpoints
// of the interval.  The operator() double(a,b) actually performs the integration
// using the R version of quadpack, Rdqa gs, which is exported from the R interface
class Int_g_1_m_alpha_g_exp_m_g {
private:
  g_class *param;
  double abs_tol;
  double rel_tol;
  int limit;
  int verbose;
  static std::string msgs[];

public:
  Int_g_1_m_alpha_g_exp_m_g(g_class* param, double tol, int subdivisions, int verbose) :
    param(param), abs_tol(tol), rel_tol(tol), limit(subdivisions), verbose(verbose)
    {}

  double result;
  double abserr;
  int neval;
  int ier;
  double operator() (double lower, double upper, int npts2, double pts_in[]) {
    std::array<double, 100> points;
    for (int i=0; i<npts2-2; i++)
      points[i]=pts_in[i];
    std::array<subinterval, LIMIT> subs;
    int last;

    if (verbose>=4){
      Rcout << std::endl
            << "dqagpe(g_exp_m_g, " << "with param" << *param
            << "lower = " << lower
            << ", upper = " << upper << ", " <<std::endl;
      for (int i=0; i<npts2-2; i++)
        Rcout << "points[" << i << "] = " << points[i] << std::endl;
      Rcout << "abs_tol = " << abs_tol
            << ", rel_tol = " << rel_tol << std::endl;
    }
    dqagpe_v2(g_1_m_alpha_g_exp_m_g, (void *) param, lower, upper, npts2, points, abs_tol, rel_tol, limit,
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
        Rcout << msgs[ier] << ":" << std::endl;
      Rcout << "Integral of g_1_m_alpha_g_exp_m_g from theta = " << lower
                << " to theta = " << upper
                << " = " << result
                << ", with absolute error = " << abserr
                << ", subintervals = " << last << std::endl
                << "rsum = " << rsum << ", esum = " << esum << std::endl;
    }
    if (verbose>=4){
      std::vector<sort_data<double> > srt_a;
      for (int i=0; i<last; i++) {
        sort_data<double> elem(subs.at(i).a,i);
        srt_a.push_back(elem);
      }
      std::sort(srt_a.begin(), srt_a.end(),sd_lt<double>);

      std::vector<sort_data<double> > srt_eord;
      for (int i=0; i<last; i++) {
        sort_data<double> elem(subs.at(i).e,i);
        srt_eord.push_back(elem);
      }
      std::sort(srt_eord.begin(), srt_eord.end(),sd_lt<double>);

      std::vector<sort_data<int> > srt_iord;
      for (int i=0; i<last; i++) {
        sort_data<int> elem(srt_eord.at(i).index,last-1-i);
        srt_iord.push_back(elem);
      }
      std::sort(srt_iord.begin(), srt_iord.end(),sd_lt<int>);

      Rcout << " "
            << std::setw(13) << std::right << "a"
            << std::setw(13) << std::right << "b"
            << std::setw(13) << std::right << "length"
            << std::setw(13) << std::right << "r"
            << std::setw(13) << std::right << "average"
            << std::setw(13) << std::right << "e"
            << std::setw(5) << std::right << "rank"
            << std::setw(6) << std::right << "level" << std::endl;
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
        Rcout << std::setw(13) << std::setprecision(5) << subs.at(j).a
              << std::setw(13) << std::setprecision(5) << subs.at(j).b
              << std::setw(13) << std::setprecision(5) << subs.at(j).b-subs.at(j).a
              << std::setw(13) << std::setprecision(5) << subs.at(j).r
              << std::setw(13) << std::setprecision(5) << subs.at(j).r/(subs.at(j).b-subs.at(j).a)
              << std::setw(13) << std::setprecision(5) << subs.at(j).e
              << std::setw(5) << srt_iord.at(j).index+1
              << std::setw(6)  << subs.at(j).level << std::endl;
      }
    }
    return result;
  }
};


std::string Int_g_1_m_alpha_g_exp_m_g::msgs[]={"OK","Maximum subdivisions reached","Roundoff error detected",
        "Bad integrand behavior","Roundoff error in the extrapolation table",
        "Integral probably divergent","Input is invalid"};

double g_class::integrate_g_1_m_alpha_g_exp_m_g(double dbltol, int subdivisions,
                    double zeta_tol, int verbose) {
  // --- ddx dstable(x, alpha, beta, ..)  for alpha < 2 ---
  // For  x = zeta, have special case formula [Nolan(1997)];
  // need to use it also for x ~= zeta, i.e., x.m.zet := |x - zeta| < delta :
  if (!isfinite(zeta)) stop("!is_finite(zeta)");
  boost::uintmax_t max_iter;
  rel_eps_tolerance tol(dbltol);

  //' g() is strictly monotone -- Nolan(1997) ["3. Numerical Considerations"]
  //'     alpha >= 1  <==>  g() is falling, ie. from Inf --> 0;  otherwise growing from 0 to +Inf

  bool do_hi=true, do_lo=true;
  double g_hi = g(th_max);
  double g_lo = g(0);
  double points[102];
  int npts2 = 2;

  if (verbose)
    Rcout << "g_hi = " << g_hi << ", g_lo = " << g_lo << std::endl;
  g_solve g_s(0., this, verbose, true);
  double theta2, g1_th2;
  if (g_hi>=g_lo && g_lo>1){
    if (verbose)
      Rcout << "Theta2 is at 0" << std::endl;
    theta2=0;
    g1_th2 = theta2;
    g_exp_m_g(&g1_th2,1,this);
    do_lo=false;
  } else if (g_lo>g_hi && g_hi>1){
    if (verbose)
      Rcout << "Theta2 is at th_max" << std::endl;
    theta2=th_max;
    g1_th2 = theta2;
    g_exp_m_g(&g1_th2,1,this);
    do_hi=false;
  } else {
    if (verbose)
      Rcout << "Theta2 is in the interior" << std::endl;

    std::pair<double, double> ur1_pair;
    // toms748_solve passes lower, upper and max_iter by reference and changes them.
    max_iter = 1000;
    double upper=th_max;
    double lower=0;
    ur1_pair = boost::math::tools::toms748_solve(g_s, lower, upper, tol, max_iter);
    theta2 = (ur1_pair.first+ur1_pair.second)/2;
    if (max_iter==1000 || fabs(g(theta2)-1)>.01) {
      return 0;  // Causes use of dPareto or
    }
    g1_th2 = theta2;
    g_exp_m_g(&g1_th2,1,this);
    points[npts2-2]=theta2;
    npts2++;
  }
  double th;
  if (verbose>=2)
    Rcout << std::endl << "theta2 = " << theta2
          << ", g1(theta2) = " << g1_th2 << ", iterations = " << max_iter << std::endl;
  double ln_g_lo[]={-2,-4,-6,-8,-10,-12,-14,-16,-18,-20,-22,-24,-26,-28,-30,-32,-34,-36,-38,-40};
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
        Rcout << "target g = " << exp(ln_g[j]) << std::endl;
      }
      g_s.set_value(ln_g[j]);
      std::pair<double,double> th_pair;
      max_iter=1000;
      double lower = 0;
      double upper = th;
      if ((log(g(th))-ln_g[j])*(log(g_lo) - ln_g[j]) >=0) continue;
      th_pair = boost::math::tools::toms748_solve(g_s, lower, upper, tol, max_iter);
      if (max_iter==1000) break;
      double th_new=(th_pair.first+th_pair.second)/2;
      if (th_new==0) break;
      if (fabs((th-th_new))< 100*std::numeric_limits<double>::epsilon()*(th+th_new)) continue;
      points[npts2-2]=th=th_new;
      npts2++;
      if (verbose>=2){
        Rcout << "theta = " << th
              << ", g(theta1) = " << g(th) << ", Iterations = " << max_iter << std::endl;
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
        Rcout << "target g = " << exp(ln_g[j]) << std::endl;
      }
      g_s.set_value(ln_g[j]);
      std::pair<double,double> th_pair;
      max_iter=1000;
      double lower = th;
      double upper = th_max;
      if ((log(g(th))-ln_g[j])*(log(g_hi) - ln_g[j]) >=0) continue;
      th_pair = boost::math::tools::toms748_solve(g_s, lower, upper, tol, max_iter);
      if (max_iter==1000) break;
      double th_new=(th_pair.first+th_pair.second)/2;
      if (fabs(th_new-th_max)<100*std::numeric_limits<double>::epsilon()*(th_new+th_max)) break;
      if (fabs((th-th_new))< 100*std::numeric_limits<double>::epsilon()*(th+th_new)) continue;
      points[npts2-2]=th=th_new;
      npts2++;
      if (verbose>=2){
        Rcout << "theta = " << th
              << ", g(theta) = " << g(th) << ", Iterations = " << max_iter << std::endl;
      }
    }
  }
  double r;
  Int_g_1_m_alpha_g_exp_m_g int_g1(this, dbltol, subdivisions, verbose);
  r=int_g1(0, th_max, npts2, points);
  if (verbose) {
    if (verbose)
      Rcout << "integrate_g_1_m_alpha_g_exp_m_g(" << x_input << " , "
            <<  zeta << ",..): c_ddx*sum(r)= "
                << c_ddx << " * " << r
                << "= " << c_ddx*(r)
                << ", abs.err = " << c_ddx*int_g1.abserr << std::endl;
  }
  return c_ddx * (r);
} //integrate_g_1_m_alpha_g_exp_m_g


double g_class::ddx_sdstable1(double x, double tol, double zeta_tol, int subdivisions, int verbose)
{
  Rcout.precision(20);
  Rcout.setf(std::ios::scientific, std::ios::floatfield);
  if (!isfinite(x)) {
    if (verbose)
      Rcout << "ddx_dstable1: x is not finite.  Returning 0" << std::endl;
    return 0;
  }
  set_x(x);
  if (verbose)
    Rcout << "ddx_sdstable1" << std::endl << *this;
  double ret;
  double dfdx_zeta;

  if (alpha != 1) { // 0 < alpha < 2	&  |beta| <= 1 from above
    // General Case
    if (verbose)
      Rcout << std::endl << "ddx_dstable(., alpha=" << alpha <<  ", beta=" << beta
                << ",..): --> theta0 = " << theta0
                <<", zeta = " << zeta << std::endl;
    bool finSupp = (fabs(beta) == 1 && alpha < 1);
    if(finSupp) {
      // has *finite* support    [zeta, Inf)  if beta ==  1
      //                         (-Inf, zeta] if beta == -1
      if((beta_input == 1 && x <= zeta) || (beta_input == -1 && x >= zeta))
        return 0;
    }
    if (zeta_tol == 0) {
      zeta_tol = 2e-15;
      if (verbose)
        Rcout << " --> zeta.tol = " << zeta_tol << std::endl;
    }
    // For  x = zeta, have special case formula [Nolan(1997)];
    // need to use it also for x ~= zeta, i.e., x.m.zet := |x - zeta| < delta :
    dfdx_zeta = tgamma(1 + 2/alpha) * sin(2*theta0_x_gt_zeta)
                /(2*M_PI * pow(1 + pow(zeta,2),(1/(alpha))));
    if (isfinite(x) && x_m_zet <= zeta_tol * (zeta_tol + std::max(fabs(x),
                                              fabs(zeta)))) {
      if (verbose)
        Rcout << "ddx_sdstable1(" << x
              << "~=" << zeta
              << " Using dfdx_zeta()" << std::endl;
        return dfdx_zeta;
    }
    // the real check should be about the feasibility of g() below, or its integration

  } // alpha != 1
  ret = integrate_g_1_m_alpha_g_exp_m_g(tol, subdivisions, zeta_tol, verbose);
  return ret;

} //ddx_sdstaple1

