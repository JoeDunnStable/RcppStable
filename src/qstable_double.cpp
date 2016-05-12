//
//  qstable_double.cpp
//  
//
//  Created by Joseph Dunn on 1/2/16.
//
//

#include "stable_double.h"
#include "q_guess.h"
#include <boost/math/tools/toms748_solve.hpp>
#include <boost/math/special_functions/erf.hpp>

using std::pair;
using boost::math::erf_inv;
using boost::math::erfc_inv;

// Functor passed to toms748 solve to find q for unit stable distribution
class p_double_solve {
private:
    double p;
    g_double_class* param;
    int lower_tail;
    int log_p;
    int maxiter;
    
public:
    p_double_solve(double p, g_double_class* param,
            int lower_tail, int log_p) :
    p(p), param(param), lower_tail(lower_tail) , log_p(log_p){}
    double operator()(const double q) {
        if (param->verbose)
            cout << "Calling spstable with parmerters" << endl
            << "q = " << q << ", alpha = " << param->alpha
            << ", beta = " << param->beta_input << endl
            << "targeting p = " << p;
        
        double r = param->spstable1(q, lower_tail, log_p)-p;
        if (param->verbose)
            cout << ", Resulting delta = " << r << endl;
        return r;
    }
};

namespace boost { namespace math { namespace tools {
    
    template <class F, class T, class Tol, class Policy>
    pair<T, T> bracket_and_solve_root2(F f, const T& guess, T factor, bool rising, Tol tol, boost::uintmax_t& max_iter, const Policy& pol)
    {
        BOOST_MATH_STD_USING
        static const char* function = "bracket_and_solve_root2<%1%>";
        //
        // Set up inital brackets:
        //
        T a = guess;
        T b = a;
        T fa = f(a);
        T fb = fa;
        //
        // Set up invocation count:
        //
        boost::uintmax_t count = max_iter - 1;
        
        int step = 32;
        
        if((fa < 0) == rising)
        {
            //
            // Zero is to the right of b, so walk upwards
            // until we find it:
            //
            while((boost::math::sign)(fb) == (boost::math::sign)(fa))
            {
                if(count == 0)
                    return boost::math::detail::pair_from_single(policies::raise_evaluation_error(function, "Unable to bracket root, last nearest value was %1%", b, pol));
                //
                // Heuristic: normally it's best not to increase the step sizes as we'll just end up
                // with a really wide range to search for the root.  However, if the initial guess was *really*
                // bad then we need to speed up the search otherwise we'll take forever if we're orders of
                // magnitude out.  This happens most often if the guess is a small value (say 1) and the result
                // we're looking for is close to std::numeric_limits<T>::min().
                //
                if((max_iter - count) % step == 0)
                {
                    factor *= 2;
                    if(step > 1) step /= 2;
                }
                //
                // Now go ahead and move our guess by "factor":
                //
                a = b;
                fa = fb;
                b += factor;
                fb = f(b);
                --count;
                BOOST_MATH_INSTRUMENT_CODE("a = " << a << " b = " << b << " fa = " << fa << " fb = " << fb << " count = " << count);
            }
        }
        else
        {
            //
            // Zero is to the left of a, so walk downwards
            // until we find it:
            //
            while((boost::math::sign)(fb) == (boost::math::sign)(fa))
            {
                if(count == 0)
                    return boost::math::detail::pair_from_single(policies::raise_evaluation_error(function, "Unable to bracket root, last nearest value was %1%", a, pol));
                //
                // Heuristic: normally it's best not to increase the step sizes as we'll just end up
                // with a really wide range to search for the root.  However, if the initial guess was *really*
                // bad then we need to speed up the search otherwise we'll take forever if we're orders of
                // magnitude out.  This happens most often if the guess is a small value (say 1) and the result
                // we're looking for is close to numeric_limits<T>::min().
                //
                if((max_iter - count) % step == 0)
                {
                    factor *= 2;
                    if(step > 1) step /= 2;
                }
                //
                // Now go ahead and move are guess by "factor":
                //
                b = a;
                fb = fa;
                a -= factor;
                fa = f(a);
                --count;
                BOOST_MATH_INSTRUMENT_CODE("a = " << a << " b = " << b << " fa = " << fa << " fb = " << fb << " count = " << count);
            }
        }
        max_iter -= count;
        max_iter += 1;
        pair<T, T> r = toms748_solve(
                                     f, a, b, fa, fb,
                                     tol, count, pol);
        max_iter += count;
        BOOST_MATH_INSTRUMENT_CODE("max_iter = " << max_iter << " count = " << count);
        BOOST_MATH_LOG_COUNT(max_iter)
        return r;
    }
    
    template <class F, class T, class Tol>
    inline pair<T, T> bracket_and_solve_root2(F f, const T& guess, const T& factor, bool rising, Tol tol, boost::uintmax_t& max_iter)
    {
        return bracket_and_solve_root2(f, guess, factor, rising, tol, max_iter, policies::policy<>());
    }
    
} // namespace tools
} // namespace math
} // namespace boost

double q_guess(double p, double alpha,double  beta, int lower_tail,int log_p);

// Returns a quantile for the unit stable distribution.
double g_double_class::sqstable1(double p, int lower_tail, int log_p, double dbltol) {
    
    if (lower_tail) {
        if (fabs(log_p ? exp(p) : p) < dbltol)
                return NegInf;
        else if (fabs("1"-log_p ? exp(p) : p) < dbltol)
            return PosInf;
    } else {
        if (fabs(log_p ? exp(p) : p) < dbltol)
            return PosInf;
        else if (fabs("1"-log_p ? exp(p) : p) < dbltol)
            return NegInf;
    }
    if (alpha == 1 && beta_input == 0) {
        if (verbose)
            cout << "Returning inverse cdf of the Cauchy distribution. " << endl;
        neval = 0;
        abserr = 0;
        return lower_tail ? (log_p ? tan((exp(p)-.5)*pi) : tan((p-.5)*pi))
                          : (log_p ? -tan((exp(p)-.5)*pi) : -tan((p-.5)*pi)) ;
    }
    if (alpha == 2) {
        if (verbose)
            cout << "Returning inverse cdf of the normal distribution. " << endl;
        neval = 0;
        abserr = 0;
        return lower_tail ? (log_p ? 2*erf_inv(2*exp(p)-1) : 2*erf_inv(2*p-1))
                          : (log_p ? -2*erf_inv(2*exp(p)-1) : -2*erf_inv(2*p-1));
    }
    p_double_solve p_s(p, this, lower_tail, log_p);
    pair<double,double> r;
    rel_eps_tolerance_double tol(dbltol);
    double guess = q_guess(p,alpha,beta_input,lower_tail,log_p);
    if (verbose)
        cout << "Guess for q " << guess << endl;
    double factor = 2.;
    bool rising = lower_tail;
    boost::uintmax_t maxiter = 1000;
    r=boost::math::tools::bracket_and_solve_root2(p_s,guess,factor,rising,tol,maxiter);
    if (verbose)
        cout << "r.first = " << r.first << ", r.second - " << r.second
        << " in " << maxiter << " iterations" << endl;
    return (r.first+r.second)/2;
}

