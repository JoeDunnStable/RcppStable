//
//  sdpqstable.cpp
//  lib_stable_double
//
//  Created by Joseph Dunn on 4/27/16.
//  Copyright © 2016 Joseph Dunn. All rights reserved.
//

#include <iostream>
using std::cout;
using std::endl;

#include <iomanip>
using std::setw;
using std::setprecision;
using std::fixed;
using std::right;
using std::scientific;

#include <algorithm>
using std::sort;

#include <fstream>
using std::ofstream;

#include <sstream>
using std::stringstream;

#include <string>
using std::string;


//#define BOOST_MATH_INSTRUMENT
#include <boost/math/tools/toms748_solve.hpp>
using boost::math::tools::toms748_solve;

#include <ctime>

#include "stable_double.h"
#include "sdpqstable_double.h"

#include "neldermeadsolver.h"
using cppoptlib::Problem;
using cppoptlib::NelderMeadSolver;

Vec sdstable(const Vec x, const double alpha, const double beta, const int log_flag,
             integr_ctl_double& ctl, double zeta_tol, int verbose){
    Vec ret(x.size());
    int i;
    g_double_class param(alpha, beta, ctl, zeta_tol, verbose);
    for (i=0; i<x.size(); i++)
        ret(i)=param.sdstable1(x(i), log_flag);
    return ret;
}


namespace boost { namespace math {

    template <class RealType, class Policy>
    inline RealType mycdf(const students_t_distribution<RealType, Policy>& dist, const RealType& x)
    {
        if((boost::math::isinf)(x))
        {
            if(x < 0) return 0; // -infinity
            return 1; // + infinity
        }
        return cdf(dist, x);
    }

    template <class RealType, class Policy>
    inline RealType mycdf(const complemented2_type<students_t_distribution<RealType, Policy>, RealType>& c)
    {
        return mycdf(c.dist, -c.param);
    }

} } // Namespace boost::math

using boost::math::mycdf;
g_quick_double::g_quick_double(g_double_class *param): param(param), dist_t(param->alpha){
    // the splines don't work well near the mode so we'll use dstable for everything near the mode
    double x_mode=stable_mode(param->alpha, param->beta_input, 1e-9, *(param->ctl_final), param->zeta_tol, param->verbose);
    double p_break_0 = 1e-6;
    x_break_0 = param->sqstable1(p_break_0, true, false, 1e-9);
    double x_break_1 = fmax(x_break_0,x_mode-10);
    double x_break_2 = fmax(x_break_1,x_mode-1);
    x_break_3 = fmax(x_break_2, x_mode-.1);
    //        double p_break_7 = fmin(p_high, 1e-6);
    double p_break_7 = 1e-6;
    //        double x_break_7 = (p_break_7 < p_high) ? param.sqstable1(p_break_7, false, false, 1e-9)
    //                                                : x_high;
    x_break_7 = param->sqstable1(p_break_7, false, false, 1e-9);
    double x_break_6 = fmin(x_break_7, x_mode+10);
    double x_break_5 = fmin(x_break_6, x_mode+1);
    x_break_4 = fmin(x_break_5, x_mode+.1);

    double pt_break_0 = mycdf(dist_t, x_break_0);
    double pt_break_1 = mycdf(dist_t,x_break_1);
    double pt_break_2 = mycdf(dist_t,x_break_2);
    double pt_break_3 = mycdf(dist_t,x_break_3);
    unsigned int n_knots_1 = (pt_break_0 < pt_break_1) ? 10 : 0;
    unsigned int n_knots_2 = (pt_break_1 < pt_break_2) ? 50 : 0;
    unsigned int n_knots_3 = (pt_break_2 < pt_break_3) ? 30 : 0;
    unsigned int n_knots = n_knots_1 + n_knots_2 + n_knots_3 + 1;
    Vec pt_knot(n_knots);
    int j=0;
    double del_pt = (n_knots_1 > 0) ? (pt_break_1-pt_break_0)/n_knots_1 : 0;
    pt_knot(j) = pt_break_0;
    j = 1;
    for (int i=0; i<n_knots_1; i++) {
        pt_knot(j) = pt_knot(j-1)+del_pt;
        j++;
    }
    pt_knot(j-1) = pt_break_1; // Otherwise it's slightly off because of rounding
    double del_x = (n_knots_2 > 0) ? (x_break_2-x_break_1)/n_knots_2: 0.;
    double x_last=x_break_1;
    for (int i=0; i<n_knots_2; i++) {
        x_last += del_x;
        pt_knot(j) = mycdf(dist_t,x_last);
        j++;
    }
    pt_knot(j-1) = pt_break_2;  //It's otherwise slightly off because of rounding
    del_x = (n_knots_3 > 0) ? (x_break_3 - x_break_2)/n_knots_3 : 0;
    x_last=x_break_2;
    for (int i=0; i<n_knots_3; i++) {
        x_last += del_x;
        pt_knot(j) = mycdf(dist_t, x_last);
        j++;
    }
    pt_knot(n_knots-1) = pt_break_3;  //It's otherwise slightly off because of rounding
    Vec y_knot(n_knots);
    j=0;
    for (unsigned int i=0; i<n_knots; ++i) {
        if (pt_knot(i) != 0 && pt_knot(i) != 1) {
            double x_tmp = quantile(dist_t,pt_knot(i));
            double y_tmp = param->sdstable1(x_tmp,true);
            if (isfinite(y_tmp)) {
                pt_knot(j)=pt_knot(i);
                y_knot(j) = y_tmp - log(pdf(dist_t, x_tmp));
                ++j;
            }
        }
    }
    pt_knot.conservativeResize(j);
    y_knot.conservativeResize(j);
    if (j>=2) {
        x_break_0 = quantile(dist_t,pt_knot(0));
        x_break_3 = quantile(dist_t,pt_knot(j-1));
    } else {
        x_break_0 = x_break_3;  // Do not use splinelow
    }
    Vec der(2);  der << 0,0;
    spline_low = cubicspline(pt_knot, y_knot, false, der);

    double pt_break_7 = mycdf(complement(dist_t, x_break_7));
    double pt_break_6 = mycdf(complement(dist_t,x_break_6));
    double pt_break_5 = mycdf(complement(dist_t,x_break_5));
    double pt_break_4 = mycdf(complement(dist_t,x_break_4));
    n_knots_1 = (pt_break_6 > pt_break_7) ? 10 : 0;
    n_knots_2 = (pt_break_5 > pt_break_6) ? 50 : 0;
    n_knots_3 = (pt_break_4 > pt_break_5) ? 30 : 0;
    n_knots = n_knots_1 + n_knots_2 + n_knots_3 + 1;
    pt_knot.resize(n_knots);
    y_knot.resize(n_knots);
    j=0;
    del_pt = (n_knots_1 > 0) ? (pt_break_6-pt_break_7)/n_knots_1 : 0;
    pt_knot(j) = pt_break_7;
    j = 1;
    for (int i=0; i<n_knots_1; i++) {
        pt_knot(j) = pt_knot(j-1)+del_pt;
        j++;
    }
    pt_knot(j-1) = pt_break_6; // Otherwise it's slightly off because of rounding
    del_x = (n_knots_2 > 0) ? (x_break_5-x_break_6)/n_knots_2: 0.;
    x_last=x_break_6;
    for (int i=0; i<n_knots_2; i++) {
        x_last += del_x;
        pt_knot(j) = mycdf(complement(dist_t,x_last));
        j++;
    }
    pt_knot(j-1) = pt_break_5;  //It's otherwise slightly off because of rounding
    del_x = (n_knots_3 > 0) ? (x_break_4 - x_break_5)/n_knots_3 : 0;
    x_last=x_break_5;
    for (int i=0; i<n_knots_3; i++) {
        x_last += del_x;
        pt_knot(j) = mycdf(complement(dist_t, x_last));
        j++;
    }
    pt_knot(n_knots-1) = pt_break_4;  //It's otherwise slightly off because of rounding
    j=0;
    for (unsigned int i=0; i<n_knots; ++i) {
        if (pt_knot(i) != 0 && pt_knot(i) != 1) {
            double x_tmp = quantile(complement(dist_t,pt_knot(i)));
            double y_tmp = param->sdstable1(x_tmp,true);
            if (isfinite(y_tmp)) {
                pt_knot(j)=pt_knot(i);
                y_knot(j) = y_tmp - log(pdf(dist_t,x_tmp));
                ++j;
            };
        }
    }
    pt_knot.conservativeResize(j);
    y_knot.conservativeResize(j);
    if (j>=2) {
        x_break_7 = quantile(complement(dist_t,pt_knot(0)));
        x_break_4 = quantile(complement(dist_t,pt_knot(j-1)));
    } else {
        x_break_7 = x_break_4;
    }
    spline_high = cubicspline(pt_knot, y_knot, false, der);
}

    // This function returns an approximation to the log likelihood of the observations x
Vec g_quick_double::operator() (const Vec& x)
{
    unsigned int n = static_cast<unsigned int>(x.size());
    Vec ret(n);
    if (param->alpha < 2 - g_double_class::Machine_eps*10){
        pMat idx(sort_indexes(x));
        Vec xx = idx.inverse() * x;  //xx contains the values from x sorted in ascending order
        unsigned int n0;
        for (n0=0; xx(n0)<x_break_0; ++n0)
        {
            ret(n0) = param->sdstable1(xx(n0), true);
        }
        unsigned int n3;
        for (n3=n; n3 > 1 && xx(n3-1) > x_break_7; n3--){
            ret(n3-1) = param->sdstable1(xx(n3-1),true);
        }
        unsigned int n1;
        for (n1=n0; n1<n-1 && xx(n1)<x_break_3; ++n1);
        Vec pt_x_low(n1-n0);
        Vec ln_dt_x_low(n1-n0);
        for (unsigned int i = n0; i<n1; ++i) {
            pt_x_low(i-n0) = mycdf(dist_t,xx(i));
            ln_dt_x_low(i-n0) = log(pdf(dist_t, xx(i)));
        }
        ret.segment(n0,n1-n0) = spline_low(pt_x_low) + ln_dt_x_low;
        unsigned int n2;
        for (n2=n1; n2 < n3 && xx(n2)<x_break_4; ++n2) {
            ret(n2) = param->sdstable1(xx(n2),true);
        }
        Vec pt_x_high(n3-n2);
        Vec ln_dt_x_high(n3-n2);
        for (unsigned int i = n2; i<n3; ++i) {
            pt_x_high(i-n2) = mycdf(complement(dist_t,xx(i)));
            ln_dt_x_high(i-n2) = log(pdf(dist_t, xx(i)));
        }
        ret.segment(n2,n3-n2)= spline_high(pt_x_high) + ln_dt_x_high;
        return idx * ret;
    } else {
        for (int i=0; i<n; ++i) {
            ret(i) = -x(i)*x(i)/4 -log(2) -log(M_PI)/2;
        }
        return ret;
    }
}

double g_quick_double::operator() (const double& x)
{
    if (param->alpha < 2 - g_double_class::Machine_eps*10){
        if (x<x_break_0) {
            return param->sdstable1(x, true);
        } else if (x >= x_break_7){
            return param->sdstable1(x,true);
        } else if (x >= x_break_0 && x < x_break_3) {
            Vec pt_x_low(1);
            double ln_dt_x_low;
            pt_x_low(0) = mycdf(dist_t,x);
            ln_dt_x_low =log(pdf(dist_t,x));
            return spline_low(pt_x_low)(0) + ln_dt_x_low;
        } else if (x >= x_break_3 && x < x_break_4) {
            return param->sdstable1(x,true);
        } else {
            Vec pt_x_high(1);
            pt_x_high(0) = mycdf(complement(dist_t,x));
            double ln_dt_x_high = log(pdf(dist_t,x));
            return spline_high(pt_x_high)(0) + ln_dt_x_high;
        }
    } else {
        return -x*x/4 -log(2) -log(M_PI)/2;
    }
}

Vec g_quick_double::pt(const Vec& x) {
    Vec ret(x.size());
    for (int i=0; i<x.size(); i++) {
        if (x(i) >=x_break_0 && x(i) < x_break_3) ret(i) = mycdf(dist_t,x(i));
        else if (x(i)>=x_break_4 && x(i) < x_break_7) ret(i) = mycdf(complement(dist_t,x(i)));
        else ret(i) = NAN;
    }
    return ret;
}

Vec sdstable_quick(const Vec x, const double alpha, const double beta, const int log_flag,
             integr_ctl_double& ctl, double zeta_tol, int verbose){
  g_double_class param(alpha, beta, ctl, zeta_tol, verbose);
  g_quick_double param_quick(&param);
  Vec ll = param_quick(x);
  return (log_flag) ? ll : ll.array().exp();
}

Vec spstable(Vec z, double alpha, double beta, int lower_tail, int log_p,
             integr_ctl_double& ctl, double zeta_tol, int verbose) {
    Vec ret(z.size());
    int i;
    g_double_class param(alpha, beta, ctl, zeta_tol, verbose);
    for (i=0; i<z.size(); i++)
        ret(i)=param.spstable1(z(i), lower_tail, log_p);
    return ret;
}

Vec sqstable(const Vec p, const double alpha, const double beta, const int lower_tail, const int log_p,
             const double dbltol, integr_ctl_double& ctl,
             const double zeta_tol, const int verbose) {
    Vec ret(p.size());
    int i;
    g_double_class param(alpha, beta, ctl, zeta_tol, verbose);
    for (i=0; i<p.size(); i++)
        ret(i)=param.sqstable1(p(i), lower_tail, log_p, dbltol);
    return ret;
}

Vec ddx_sdstable(const Vec x, const double alpha, const double beta,
             integr_ctl_double& ctl, double zeta_tol, int verbose){
  Vec ret(x.size());
  int i;
  g_double_class param(alpha, beta, ctl, zeta_tol, verbose);
  for (i=0; i<x.size(); i++)
    ret(i)=param.ddx_sdstable1(x(i));
  return ret;
}

double capped_sdstable(const Vec y, const double alpha, const double beta, const double gamma, const double delta,
                       const bool quick, integr_ctl_double &ctl, const double zeta_tol, const int verbose) {
    g_double_class param(alpha, beta, ctl, zeta_tol, verbose);
    Vec yy = (y.array() - delta)/gamma;
    double log_gamma = log(gamma);
    double out =0;
    double p_low = param.spstable1(-1e300, true, true);
    double p_high = param.spstable1(1e300, false, true);
    if (quick) {
        g_quick_double sdstable_quick(&param);
        int j = 0;
        for (int i=0; i<yy.size(); i++) {
            if (yy(i) > 1e300)
                out += p_high;
            else if (yy(i) < -1e300)
                out += p_low;
            else {
                yy(j) = yy(i);
                j++;
            }
        }
        yy.conservativeResize(j);
        out += sdstable_quick(yy).sum() - j*log_gamma;
    } else {
        for (int i=0; i<yy.size(); i++) {
            if (yy(i) > 1e300)
                out += p_high;
            else if (yy(i) < -1e300)
                out += p_low;
            else
                out += param.sdstable1(yy(i),true) - log_gamma;
        }
    }
    out = -2 * out / y.size();
    return out;
}

class McCullochFit : public Problem<double> {
public:
    double q_kurt;
    double q_skew;
    double alpha_min;
    double alpha_max;
    ostream *trace;
    double dbltol;
    integr_ctl_double &ctl;
    double zeta_tol;
    int verbose;
    double value(const Vec &par) {
        double alpha = par_to_alpha(par(0));
        double beta = atan(par(1))/M_PI_2;
        *trace << setw(18) << setprecision(8) << alpha
               << setw(18) << setprecision(8) << beta;
       (*trace).flush();
        g_double_class param(alpha, beta, ctl, zeta_tol, verbose);
        Vec probs(5); probs << .05, .25, .5, .75, .95;
        bool lower_tail = true;
        bool log_p = false;
        Vec qs = sqstable(probs, alpha, beta, lower_tail, log_p, dbltol, ctl, zeta_tol, verbose);
        double qs_kurt = (qs(4) - qs(0))/(qs(3)-qs(1));
        double qs_skew = (qs(4) -2 * qs(2) + qs(0))/(qs(4)-qs(0));
        *trace << setw(18) << setprecision(8) << qs_kurt
        << setw(18) << setprecision(8) << qs_skew;
        double out = pow(q_kurt/qs_kurt-1,2)+100*pow(qs_skew-q_skew,2);
        *trace << setw(18) << setprecision(8) << out << endl;
        return fmin(out,1e100);  //Optim doesn't like infinite numbers
    }
    double alpha_to_par(const double alpha) {
        double a = fmin(alpha_max,fmax(alpha_min,alpha));
        return tan((alpha_max-a)/(alpha_max-alpha_min)*(-M_PI_2)+(a-alpha_min)/(alpha_max-alpha_min)*M_PI_2);

    }
    double par_to_alpha(const double par) {
        return alpha_min + (alpha_max-alpha_min)*(.5+atan(par)/M_PI);
    }
    double norm(const Vec& par1, const Vec& par2) {
        Vec del(2);
        del(0) = par_to_alpha(par1(0)) - par_to_alpha(par2(0));
        del(1) = (atan(par1(1))-atan(par2(1)))/M_PI_2;
        return del.array().abs().maxCoeff();
    }

    McCullochFit(const double q_kurt, const double q_skew, const double alpha_min, const double alpha_max, ostream *trace, const double dbltol, integr_ctl_double &ctl, const double zeta_tol, const int verbose )
    : q_kurt(q_kurt), q_skew(q_skew), alpha_min(alpha_min), alpha_max(alpha_max), trace(trace), dbltol(dbltol), ctl(ctl), zeta_tol(zeta_tol), verbose(verbose) {};
};

class DunnFit : public Problem<double> {
public:
    double q_kurt;
    double q_mode;
    double skew;
    double alpha_min;
    double alpha_max;
    ostream *trace;
    double dbltol;
    integr_ctl_double &ctl;
    double zeta_tol;
    int verbose;
    double value(const Vec &par) {
        double alpha = par_to_alpha(par(0));
        double beta = atan(par(1))/M_PI_2;
        *trace << setw(18) << setprecision(8) << alpha
        << setw(18) << setprecision(8) << beta;
        (*trace).flush();
        g_double_class param(alpha, beta, ctl, zeta_tol, verbose);
        Vec probs(5); probs << .05, .25, .5, .75, .95;
        bool lower_tail = true;
        bool log_p = false;
        Vec qs = sqstable(probs, alpha, beta, lower_tail, log_p, dbltol, ctl, zeta_tol, verbose);
        double mode = stable_mode(alpha, beta, dbltol, ctl, zeta_tol, verbose);
        double q_delta = .5*(qs(4)-qs(0));
        double p_minus = param.spstable1(mode-q_delta,true,false);
        double p_plus = param.spstable1(mode + q_delta, false, false);
        double s_skew = (p_plus - p_minus)/(p_plus + p_minus);
        double qs_kurt = (qs(4) - qs(0))/(qs(3)-qs(1));
        *trace << setw(18) << setprecision(8) << qs_kurt
        << setw(18) << setprecision(8) << s_skew;
        double out = pow(q_kurt/qs_kurt-1,2)+pow(s_skew-skew,2);
        *trace << setw(18) << setprecision(8) << out << endl;
        return fmin(out,1e100);  //Optim doesn't like infinite numbers
    }
    double alpha_to_par(const double alpha) {
        double a = fmin(alpha_max,fmax(alpha_min,alpha));
        return tan((alpha_max-a)/(alpha_max-alpha_min)*(-M_PI_2)+(a-alpha_min)/(alpha_max-alpha_min)*M_PI_2);

    }
    double par_to_alpha(const double par) {
        return alpha_min + (alpha_max-alpha_min)*(.5+atan(par)/M_PI);
    }
    double norm(const Vec& par1, const Vec& par2) {
        Vec del(2);
        del(0) = par_to_alpha(par1(0)) - par_to_alpha(par2(0));
        del(1) = (atan(par1(1))-atan(par2(1)))/M_PI_2;
        return del.array().abs().maxCoeff();
    }

    DunnFit(const double q_kurt, const double q_mode, const double skew, const double alpha_min, const double alpha_max, ostream *trace, const double dbltol, integr_ctl_double &ctl, const double zeta_tol, const int verbose )
    : q_kurt(q_kurt), q_mode(q_mode), skew(skew), alpha_min(alpha_min), alpha_max(alpha_max), trace(trace), dbltol(dbltol), ctl(ctl), zeta_tol(zeta_tol), verbose(verbose) {};
};

class MLEFit : public Problem<double> {
public:
    Vec y;
    double alpha_min;
    double alpha_max;
    bool quick;
    ostream *trace;
    integr_ctl_double &ctl;
    double zeta_tol;
    int verbose;
    double value(const Vec &par) {
        double alpha = par_to_alpha(par(0));
        double beta = atan(par(1))/M_PI_2;
        double gamma = exp(par(2));
        double delta = par(3);
        *trace << setw(18) << setprecision(8) << alpha
        << setw(18) << setprecision(8) << beta
        << setw(18) << setprecision(8) << gamma
        << setw(18) << setprecision(8) << delta;
        (*trace).flush();
        double out = 0;
        Vec yy = (y.array()-delta)/gamma;
        g_double_class param(alpha, beta, ctl, zeta_tol, verbose);
        double p_high = param.spstable1(1e300, false, true);
        double p_low = param.spstable1(-1e300, true, true);
        double log_gamma = log(gamma);
        if (quick) {
            g_quick_double sdstable_quick(&param);
            int j = 0;
            for (int i=0; i<yy.size(); i++) {
                if (yy(i) > 1e300)
                    out += p_high;
                else if (yy(i) < -1e300)
                    out += p_low;
                else {
                    yy(j) = yy(i);
                    j++;
                }
            }
            yy.conservativeResize(j);
            out += sdstable_quick(yy).sum() - j*log_gamma;
        } else {
            for (int i=0; i<yy.size(); i++) {
                if (yy(i) > 1e300)
                    out += p_high;
                else if (yy(i) < -1e300)
                    out += p_low;
                else
                    out += param.sdstable1(yy(i),true) - log_gamma;
            }
        }
        out = -2 * out / y.size();
        *trace << setw(18) << setprecision(8) << out << endl;
        return fmax(fmin(out,1e100),-1e+100);  //Optim doesn't like infinite numbers
    }
    double alpha_to_par(const double alpha) {
        double a = fmin(alpha_max,fmax(alpha_min,alpha));
        return tan((alpha_max-a)/(alpha_max-alpha_min)*(-M_PI_2)+(a-alpha_min)/(alpha_max-alpha_min)*M_PI_2);

    }
    double par_to_alpha(const double par) {
        return alpha_min + (alpha_max-alpha_min)*(.5+atan(par)/M_PI);
    }
    double norm(const Vec& par1, const Vec& par2) {
        Vec del(4);
        del(0) = par_to_alpha(par1(0)) - par_to_alpha(par2(0));
        del(1) = (atan(par1(1))-atan(par2(1)))/M_PI_2;
        del(2) = exp(par1(2))-exp(par2(2));
        del(3) = par1(3)-par2(3);
        return del.array().abs().maxCoeff();
    }

    MLEFit(const Vec &y, const double alpha_min, const double alpha_max, bool quick, ostream *trace, integr_ctl_double &ctl, const double zeta_tol, const int verbose )
    : y(y), alpha_min(alpha_min), alpha_max(alpha_max), quick(quick), trace(trace), ctl(ctl), zeta_tol(zeta_tol), verbose(verbose) {};
};


class QMLEFit : public Problem<double> {
public:
    Vec y;
    g_double_class param;
    bool quick;
    ostream *trace;
    g_quick_double sdstable_quick;

    double value(const Vec &par) {
        double gamma = exp(par(0));
        double delta = par(1);
        *trace << setw(18) << setprecision(8) << param.alpha
        << setw(18) << setprecision(8) << param.beta
        << setw(18) << setprecision(8) << gamma
        << setw(18) << setprecision(8) << delta;
        double out = 0;
        Vec yy = (y.array()-delta)/gamma;
        double p_high = param.spstable1(1e300, false, true);
        double p_low = param.spstable1(-1e300, true, true);
        double log_gamma = log(gamma);
        if (quick) {
            int j = 0;
            for (int i=0; i<yy.size(); i++) {
                if (yy(i) > 1e300)
                    out += p_high;
                else if (yy(i) < -1e300)
                    out += p_low;
                else {
                    yy(j) = yy(i);
                    j++;
                }
            }
            yy.conservativeResize(j);
            out += sdstable_quick(yy).sum() - j*log_gamma;
        } else {
            for (int i=0; i<y.size(); i++) {
                if (yy(i) > 1e300)
                    out += p_high;
                else if (yy(i) < -1e300)
                    out += p_low;
                else
                    out += param.sdstable1(yy(i),true) - log_gamma;
            }
        }
        out = -2 * out / y.size();
        *trace << setw(18) << setprecision(8) << out << endl;
        return fmax(fmin(out,1e100),-1e+100);  //Optim doesn't like infinite numbers
    }
    double norm(const Vec& par1, const Vec& par2) {
        Vec del(2);
        del(0) = exp(par1(0))-exp(par2(0));
        del(1) = par1(1)-par2(1);
        return del.array().abs().maxCoeff();
    }

    QMLEFit(const Vec &y, const double &alpha, const double &beta, bool quick, ostream *trace, integr_ctl_double &ctl, const double zeta_tol, const int verbose )
    : y(y), param(alpha, beta, ctl, zeta_tol, verbose), quick(quick), trace(trace), sdstable_quick(&param) {};

};

inline double to_prob(double p){
    return fmax(0.,fmin(1.,p));
}

Vec quantile(Vec &x, Vec probs)
{
    double eps = 100 * std::numeric_limits<double>::epsilon();
    long np = probs.size();
    if ((probs.array() < -eps || probs.array() > 1 + eps).any())
        throw std::range_error("quantile: 'probs' outside [0,1]");
    probs.array().unaryExpr(std::ptr_fun(to_prob));
    Vec qs(np);
    long n = x.size();
    if (n > 0 && np > 0) {
        std::sort(x.data(), x.data()+n);
        for (int j=0; j<np; ++j) {
            double index = (n - 1) * probs(j);
            int lo = floor(index);
            int hi = ceil(index);
            double h = index - lo;
            qs(j) = (1-h) * x(lo) + h * x(hi);
        }
        return qs;
    } else {
        throw std::range_error("quantile: Both x and prob must be of length greater than 0");
    }
}

double p_sample(Vec x_sample, double x ) {
    long n = x_sample.size();
    long n_lower = 0;
    for (int i=0; i<n; ++i) {
        if (x_sample(i) <= x) {
            ++n_lower;
        }
    }
    return static_cast<double>(n_lower)/n;
}

double mode_sample(Vec x_sample) {
    int n = x_sample.size();
    std::sort(x_sample.data(),x_sample.data()+n);
    int m = n/100;
    int i_mode = -1;
    double x_delta = std::numeric_limits<double>::infinity();
    for (int i=0; i+m<n; ++i) {
        if (x_sample(i+m) - x_sample(i) < x_delta) {
            i_mode = i;
            x_delta = x_sample(i+m) - x_sample(i);
        }
    }
    return x_sample(i_mode+m/2);
}

ostream& operator<< (ostream &os, const fit_result &fr) {
    os << setw(10) << fr.method
    << setw(14) << setprecision(5) << scientific << fr.alpha
    << setw(14) << setprecision(5) << fr.beta
    << setw(14) << setprecision(5) << fr.gamma
    << setw(14) << setprecision(5) << fr.delta
    << setw(14) << setprecision(5) << fr.two_ll_n
    << setw(7) << fr.n
    << setw(14) << setprecision(5) << fr.q_kurt
    << setw(14) << setprecision(5) << fr.q_skew
    << setw(14) << setprecision(5) << fr.q_scale
    << setw(14) << setprecision(5) << fr.q_location
    << setw(6) << fr.convergence
    << setw(6)  << fr.iterations
    << setw(6) << setprecision(1) << fixed << fr.cpu_time << endl;
    return os;
}

void result_heading(ostream &os) {
    os << setw(10) << right << "method"
    << setw(14) << right << "alpha"
    << setw(14) << right << "beta"
    << setw(14) << right << "gamma"
    << setw(14) << right << "delta"
    << setw(14) << right << "two_ll_n"
    << setw(7) << right << "n"
    << setw(14) << right << "q_kurt"
    << setw(14) << right << "q_skew"
    << setw(14) << right << "q_scale"
    << setw(14) << right << "q_location"
    << setw(6) << right << "conv?"
    << setw(6)  << right << "iter"
    << setw(6) << right << "time" << endl;
}

std::vector<fit_result> stable_fit(const Vec yy, integr_ctl_double ctl, const double dbltol,
                                   const string type,const bool quick, const int verbose) {
    const double zeta_tol=0;
    int n=static_cast<int>(yy.size());
    ofstream trace("stable_fit_trace.txt");

    // First McCulloch's method

    Vec probs(5); probs << .05, .25, .5, .75, .95;
    Vec y(yy);
    Vec q = quantile(y,probs);  // y is sorted as a side effect
    double alpha;
    double beta;
    string convergence;
    clock_t phase1_times;
    int iterations;
    double q_kurt=(q(4)-q(0))/(q(3)-q(1));
    double q_skew=(q(4)+q(0)-2*q(2))/(q(4)-q(0));
    double q_mode, skew;
    NelderMeadSolver<double> solver;
    NelderMeadSolver<double>::Info &myctrl = solver.ctrl();
    myctrl.iterations = 1000;
    myctrl.obj_spread = 1e-10;
    myctrl.x_spread = 1e10;   // Effectively not used.
    const NelderMeadSolver<double>::Info &myinfo = solver.info();
    if (q_kurt < 1000) {
        phase1_times = -clock();
        Vec par_McCulloch(2); par_McCulloch << 1., q_skew;
        McCullochFit mcculloch_fit(q_kurt, q_skew, .1, 2., &trace, dbltol, ctl, zeta_tol, verbose);
        solver.minimize(mcculloch_fit, par_McCulloch);
        iterations = myinfo.iterations;

        phase1_times+= clock();
        if (myinfo.obj_spread <= myctrl.obj_spread && myinfo.x_spread <= myctrl.x_spread)
            convergence = "True";
        else
            convergence = "False";
        alpha = mcculloch_fit.par_to_alpha(par_McCulloch(0));
        beta = atan(par_McCulloch(1))/M_PI_2;
    } else {
        q_mode = mode_sample(y);
        double p_upper = 1-p_sample(y,q_mode+.5*(q(4)-q(0)));
        double p_lower = p_sample(y,q_mode-.5*((q(4)-q(0))));
        skew = (p_upper-p_lower)/(p_upper + p_lower);
        phase1_times = -clock();
        Vec par_Dunn(2); par_Dunn << 1., skew;
        DunnFit dunn_fit(q_kurt, q_mode, skew, .1, 2., &trace, dbltol, ctl, zeta_tol, verbose);
        solver.minimize(dunn_fit, par_Dunn);
        iterations = myinfo.iterations;

        phase1_times += clock();
        if (myinfo.obj_spread <= myctrl.obj_spread && myinfo.x_spread <= myctrl.x_spread)
            convergence = "True";
        else
            convergence = "False";
        alpha = dunn_fit.par_to_alpha(par_Dunn(0));
        beta = atan(par_Dunn(1))/M_PI_2;
    }
    Vec ps(5);
    ps << .05, .25, .5, .75, .95;
    Vec qs = sqstable(ps, alpha, beta, true, false, dbltol, ctl, zeta_tol, verbose);
    double gamma = (q(3)-q(1))/(qs(3)-qs(1));
    double delta;
    if (q_kurt < 1000) {
        delta = q(2) - gamma * qs(2);
    } else {
        delta = q_mode - gamma * stable_mode(alpha, beta, dbltol, ctl, zeta_tol, verbose);
    }
    double alpha_min = .1, alpha_max = 2.;
    MLEFit mle_fit(y, alpha_min, alpha_max, quick, &trace, ctl, zeta_tol, verbose);
    Vec par_mle(4);
    par_mle << mle_fit.alpha_to_par(alpha), tan(M_PI_2*beta), log(gamma), delta;

    QMLEFit q_mle_fit(y, alpha, beta, quick, &trace, ctl, zeta_tol, verbose);
    Vec par_q_mle(2);
    par_q_mle << log(gamma), delta;

    qs= gamma * qs.array() + delta;
    fit_result fr_phase1((q_kurt < 1000) ? "McCulloch" : "Dunn", alpha, beta, gamma, delta, -mle_fit.value(par_mle), n, qs, convergence, iterations,
                         static_cast<double>(phase1_times)/CLOCKS_PER_SEC);

    std::vector<fit_result> results;
    results.push_back(fr_phase1);

    iterations = 0;
    if (type=="mle") {
        clock_t mle_times = -clock();
        myctrl.iterations = 1000;
        solver.minimize(mle_fit, par_mle);
        iterations += myinfo.iterations;
        if (myinfo.obj_spread <= myctrl.obj_spread && myinfo.x_spread <= myctrl.x_spread)
            convergence = "True";
        else
            convergence = "False";
        alpha = mle_fit.par_to_alpha(par_mle(0));
        beta = atan(par_mle(1))/M_PI_2;
        gamma = exp(par_mle(2));
        delta = par_mle(3);
        mle_times += clock();
        qs=gamma*sqstable(ps, alpha, beta, true, false, dbltol, ctl, zeta_tol, verbose).array() + delta;
        fit_result fr_mle("mle", alpha, beta, gamma, delta,-mle_fit.value(par_mle),n,qs,convergence,iterations,
                          static_cast<double>(mle_times)/CLOCKS_PER_SEC);
        results.push_back(fr_mle);
    }
    else if (type=="q_mle") {
        // alpha and beta are from McCulloch's method
        clock_t q_mle_times = -clock();
        myctrl.iterations=1000;
        solver.minimize(q_mle_fit, par_q_mle);
        iterations += myinfo.iterations;
        string convergence;
        if (myinfo.obj_spread <= myctrl.obj_spread && myinfo.x_spread <= myctrl.x_spread)
            convergence = "True";
        else
            convergence = "False";
        par_mle(2)=par_q_mle(0);
        gamma = exp(par_q_mle(0));
        par_mle(3) = delta = par_q_mle(1);
        q_mle_times += clock();
        qs=gamma*sqstable(ps, alpha, beta, true, false, dbltol, ctl, zeta_tol, verbose).array()+delta;
        fit_result fr_q_mle("q_mle", alpha, beta, gamma, delta, -mle_fit(par_mle), n, qs, convergence, iterations,
                            static_cast<double>(q_mle_times)/CLOCKS_PER_SEC);
        results.push_back(fr_q_mle);
    }
    return results;
}

Vec srstable(const double alpha, const double beta, const Vec &u1, const Vec &u2) {
  if (u1.size() != u2.size()) {
    throw std::range_error("rstable: u1.size() != u2.size()");
  }
  size_t n = u1.size();
  Vec ret(n);
  for (size_t i = 1; i<n; ++i) ret(i)=srstable1(alpha, beta, u1(i), u2(i));
  return ret;
}

