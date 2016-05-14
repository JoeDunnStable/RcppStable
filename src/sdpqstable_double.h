//
//  sdpqstable_double.h
//  lib_stable_double
//
//  Created by Joseph Dunn on 4/27/16.
//  Copyright © 2016 Joseph Dunn. All rights reserved.
//

#ifndef sdpqstable_double_h
#define sdpqstable_double_h

#include <Eigen/Dense>
#include "cubicspline.h"

//#define BOOST_MATH_INSTRUMENT
#include <boost/math/distributions/students_t.hpp>
using boost::math::students_t;

typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vec;

Vec sdstable(const Vec x, const double alpha, const double beta, const int log_flag,
             integr_ctl_double& ctl, double zeta_tol, int verbose);

class g_quick_double {
    double x_break_0;
    double x_break_3;
    double x_break_4;
    double x_break_7;
    g_double_class *param;
    students_t dist_t;
    cubicspline spline_low;
    cubicspline spline_high;
public:
    g_quick_double(g_double_class *param);
    Vec operator() (const Vec& x);
    double operator() (const double& x);
    Vec get_knots_low() {return spline_low.get_knots();}
    Vec get_knots_high() {return spline_high.get_knots();}
    uVec region(const Vec& x) {
        uVec ret(x.size());
        for (int i =0; i<x.size(); i++) {
            if (x(i) < x_break_0) ret(i) = 0;
            else if (x(i) < x_break_3) ret(i)=1;
            else if (x(i) < x_break_4) ret(i)=2;
            else if (x(i) < x_break_7) ret(i)=3;
            else ret(i)=4;
        }
        return ret;
    }
    Vec pt(const Vec& x);
};

Vec spstable(Vec z, double alpha, double beta, int lower_tail, int log_p,
             integr_ctl_double& ctl, double zeta_tol, int verbose);
Vec sdstable_quick(const Vec x, const double alpha, const double beta, const int log_flag,
                   integr_ctl_double& ctl, double zeta_tol, int verbose);
Vec sqstable(const Vec p, const double alpha, const double beta, const int lower_tail, const int log_p,
             const double dbltol, integr_ctl_double& ctl,
             const double zeta_tol, const int verbose);
Vec ddx_sdstable(const Vec x, const double alpha, const double beta,
                 integr_ctl_double& ctl, double zeta_tol, int verbose);
double capped_sdstable(const Vec y, const double alpha, const double beta, const double gamma, const double delta,
                       const bool quick, integr_ctl_double &ctl, const double zeta_tol, const int verbose);
Vec quantile(Vec &x, Vec probs);

class fit_result {
public:
    string method;
    double alpha;
    double beta;
    double gamma;
    double delta;
    double two_ll_n;
    int n;
    double q_kurt;
    double q_skew;
    double q_scale;
    double q_location;
    string convergence;
    unsigned int iterations;
    double cpu_time;
    fit_result(string method, double alpha, double beta, double gamma, double delta,
               double two_ll_n, int n, Vec qs, string convergence, unsigned int iterations, double cpu_time)
    : method(method), alpha(alpha), beta(beta), gamma(gamma), delta(delta),
    two_ll_n(two_ll_n), n(n), convergence(convergence), iterations(iterations),
    cpu_time(cpu_time){
        q_kurt=(qs(4)-qs(0))/(qs(3)-qs(1));
        q_skew=(qs(4)+qs(0)-2*qs(2))/(qs(4)-qs(0));
        q_scale=qs(3)-qs(1);
        q_location=qs(2);
    }
    friend ostream& operator<< (ostream &os, const fit_result &fr);
};

Vec srstable(const double alpha, const double beta, const Vec &u1, const Vec &u2);

void result_heading(ostream &os);

std::vector<fit_result> stable_fit(const Vec y, integr_ctl_double ctl, const double dbltol = 1e-10,
                                   const string type="q",const bool quick=false, const int verbose=0);

#endif /* sdpqstable_double_h */
