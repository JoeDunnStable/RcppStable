/* If using sourceCpp run Sys.setenv("PKG_CXXFLAGS"="-I /opt/local/include") */

#include "Stable.h"
#include <iostream>
#include <boost/math/tools/toms748_solve.hpp>
#include "dqagp_v2.h"

static const double Machine_eps = std::numeric_limits<double>::epsilon();
static const double PosInf = std::numeric_limits<double>::infinity();
static const double NegInf = -PosInf;

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

// Function to integrate: pstable(..)= f(..) = c2 * \int_{-\theta_0}^{\pi/2} g1_p(u) du
void exp_m_g(double * th, int n, void * ext) {
  g_class * param = (g_class *) ext;
  int i;

  for (i=0; i<n; i++){
    double g_th = param->g(th[i]);
    // g1_p :=  exp(-g(.))
    th[i] = exp(-g_th);
  }
}

// Function to integrate: pstable(..)= f(..) = c2 * \int_{-\theta_0}^{\pi/2} g1_p(u) du
void one_m_exp_m_g(double * th, int n, void * ext) {
  g_class * param = (g_class *) ext;
  int i;

  for (i=0; i<n; i++){
    double g_th = param->g(th[i]);
    // g1_p :=  exp(-g(.))
    th[i] = -expm1(-g_th);
  }
}

class Int_exp_m_g {
private:
  g_class *param;
  double abs_tol;
  double rel_tol;
  int limit;
  int verbose;
  static std::string msgs[];

public:
  Int_exp_m_g(g_class* param, double tol, int subdivisions, int verbose) :
    param(param), abs_tol(tol), rel_tol(tol), limit(subdivisions), verbose(verbose) {};

  double result;
  double abserr;
  int neval;
  int ier;
  double operator() (double lower, double upper, int npts2, double pts_in[]) {
    std::array<double, 100> points;
    for (int i=0; i<npts2-2; i++)
      points[i]=pts_in[i];
    double result;
    double abserr;
    int neval;
    int ier;
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
    dqagpe_v2(exp_m_g, (void *) param, lower, upper, npts2, points, abs_tol, rel_tol, limit,
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
      if (ier > 0)
        Rcout << msgs[ier] << ":" << std::endl;
      Rcout << "Integral of g_exp_m_g from theta = " << lower
            << " to theta = " << upper
            << " = " << result
            << ", with absolute error = " << abserr
            << ", subintervals = " << last << std::endl;
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
  };
};

std::string Int_exp_m_g::msgs[]={"OK","Maximum subdivisions reached","Roundoff error detected",
                            "Bad integrand behavior","Roundoff error in the extrapolation table",
                            "Integral probably divergent","Input is invalid"};


// ------------------------------------------------------------------------------
//' Auxiliary for pstable()  (for alpha != 1)
double g_class::integrate_exp_m_g(bool giveI, double dbltol, int subdivisions, int verbose)
{
  boost::uintmax_t max_iter;
  double lower;
  double upper;

  if(verbose)
    Rcout << "integrate_exp_m_g(" << *this << std::endl << "giveI = " << giveI << std::endl;

// as g() is montone, the integrand  exp(-g(.)) is too ==> maximum is at the boundary
// however, integration can be inaccuracte when g(.) quickly jumps from Inf to 0
  double eps=1e-10;
  double th_eps, th_1_m_eps;
  double exp_m_g_lo=0;
  exp_m_g(&exp_m_g_lo, 1, (void*)this);
  double exp_m_g_hi=th_max;
  exp_m_g(&exp_m_g_hi, 1, (void*)this);
  if (verbose)
    Rcout << "integrate_exp_m_g: exp(-g(0) = " << exp_m_g_lo << std::endl
          << ", exp(-g(" << th_max << ") = " << exp_m_g_hi << std::endl;
  bool have_eps = false, have_1_m_eps = false;
  if ((exp_m_g_hi-eps)*(exp_m_g_lo-eps)<0){
    g_solve g_s(log(-log(eps)), this, verbose, true);
    lower = 0;
    upper = th_max;
    max_iter = 1000;
    rel_eps_tolerance tol(1e-8);
    std::pair<double,double> ur;
    ur = boost::math::tools::toms748_solve(g_s, lower, upper, tol, max_iter);
    th_eps = (ur.first+ur.second)/2;
    have_eps=true;
    if(verbose)
      Rcout << " exp(-g)(" << th_eps <<") = " << exp(-g(th_eps)) << ", iterations = " << max_iter << std::endl;
  }
  if ((exp_m_g_hi-(1-eps))*(exp_m_g_lo-(1-eps))<0){
    g_solve g_s(log(-log(1-eps)), this, verbose, true);
    lower = (exp_m_g_lo==0) ? th_eps : 0;
    upper = (exp_m_g_hi==0) ? th_eps : th_max;
    max_iter = 1000;
    rel_eps_tolerance tol(1e-8);
    std::pair<double,double> ur;
    ur = boost::math::tools::toms748_solve(g_s, lower, upper, tol, max_iter);
    th_1_m_eps = (ur.first+ur.second)/2;
    have_1_m_eps=true;
    if(verbose)
      Rcout << "exp(-g)(" << th_1_m_eps <<") = " << 1-eps << ", iterations = " << max_iter << std::endl;
  }
  double points[22];
  int npts2 = 2;
  if (have_1_m_eps){
    npts2 += 2;
    if (exp_m_g_hi==0) {
      points[0]=th_1_m_eps;
      points[1]=th_eps;
    } else {
      points[0]=th_eps;
      points[1]=th_1_m_eps;
    }
  } else if (have_eps){
    npts2 += 1;
    points[0]=th_eps;
  }
  Int_exp_m_g int_exp_m_g(this, dbltol, subdivisions, verbose);
  double r;
  r=int_exp_m_g(0, th_max, npts2, points);
  if(verbose)
    Rcout << "--> Int r= " << r << std::endl;
  double c0, c1;
  if (alpha < 1){
    c0 = .5 - theta0/M_PI;
    c1 = 1/M_PI;
  } else if (alpha > 1){
    c0=1;
    c1=-1/M_PI;
  } else{ // alpha == 1
    c0 = 1;
    c1 = -.5;
  }
  if(giveI) {
    // { ==> alpha > 1 ==> c0 = 1; c1 = -1/pi}
    // return (1 - F) = 1 - (1 -1/pi * r) = r/pi :
    return -(r)*c1;
  } else {
    return c0 + c1* (r);
  }
} // g_class::integrate_exp_m_g

// ------------------------------------------------------------------------------


inline double retValue(double F, int useLower, int log_p) {
  return (useLower) ? ((log_p) ? log(F)    : F)
                    : ((log_p) ? log1p(-F) : 1 - F);

}

double g_class::spstable1(double x,
                 int lower_tail, int log_p,
                 double tol, int subdivisions, int verbose) {
  if(x==PosInf)
    return retValue(1, lower_tail, log_p);
  else if (x==NegInf)
    return retValue(0, lower_tail, log_p);
  set_x(x);
  if (alpha !=1) {
    //-------->>>  identically as in .fct1() for dstable() above: <<=----------
    // FIXME: also provide "very small alpha" case, as in .fct1()
    if (!isfinite(zeta)) stop("Zeta is not finite");
    bool finSupp = (fabs(beta) == 1 && alpha < 1);
    if(finSupp) {
      // has *finite* support    [zeta, Inf)  if beta ==  1
      //                         (-Inf, zeta] if beta == -1
      if(beta_input == 1 && x <= zeta)
        return retValue(0., lower_tail, log_p);
      else if(beta_input == -1 && x >= zeta)
        return retValue(1., lower_tail, log_p);
    // else .. one of the cases below
    }

    if(fabs(x - zeta) < 2 * Machine_eps) {
      // FIXME? same problem as dstable
      double r = (lower_tail) ? (.5 - theta0_x_gt_zeta/M_PI) : (.5 + theta0_x_gt_zeta/M_PI);
      if (verbose)
        Rcout << "x ~ zeta: " << ".5 "
                  << ((lower_tail) ? " - " : " + ") << "th0/pi = " << r << std::endl;
      return (log_p) ? log(r) : r;
    }
    int useLower = ((x > zeta && lower_tail) ||
                 (x < zeta && !lower_tail));
    // FIXME: for alpha > 1 -- the following computes F1 = 1 -c3*r(x)
    // and suffers from cancellation when 1-F1 is used below:
    bool giveI = !useLower && alpha > 1; // if true, int_exp_m_g() return 1-F
    double _F1 = integrate_exp_m_g(giveI, tol, subdivisions, verbose);
    if(giveI)
      return (log_p) ? log(_F1) : _F1;
    else
      return retValue(_F1, useLower, log_p);
  } else { // alpha = 1
    bool useL;
    if(beta_input>=0)
      useL = !lower_tail;
    else {
      useL=lower_tail;
    }

    bool giveI = !useL && !log_p;
    if(giveI)
      useL = true;
    double _F2 = integrate_exp_m_g(giveI, tol, subdivisions, verbose);
    return retValue(_F2,useL,log_p);
  }
}

