///
/// \file  adaptive_integration.h
/// \author Joseph Dunn
/// \copyright 2016 Joseph Dunn
/// \copyright Distributed under the terms of the GNU General Public License version 3

#ifndef adaptive_integration_h
#define adaptive_integration_h
#include <vector>
#include <string>
using std::string;
//#include <cstddef>
#include "myFloat.h"

namespace adaptive_integration {
  
using std::vector;
using std::ostream;
  
/// Class to hold machine constants.
/// Needed to fix problems in numeric_limits when using mpreal
class Machine {
public:
/// Initialize the static variable with machine constants
  static void initialize();
  static bool initialized;
  static int digits10;     /// the number of decimal digits
  static myFloat epsilon;    ///< the machine epsilon
  static myFloat min;     ///< the smallest absolute value
  static myFloat max;     ///< the largest absolute value
  static myFloat large_exp_arg;  ///< the maximum allowable argument to exp()
  static myFloat pi;      ///< pi in the precison of myFloat
  static myFloat pi2;     ///< pi2 in the precision of myFloat
  static myFloat PosInf;  ///< Positive infinity
  static myFloat NegInf;  ///< Negative infinity
  static void print(ostream& os);
};
  
/// vectorizing function   f(x[1:n], ...) -> x[]  {overwriting x[]}. 

typedef void Integrand(myFloat *x, ///< [in,out] points to evaluate and result of evaluation 
                       int n,     ///< [in] the number of points 
                        void *ex  ///< pointer to the environment
                        );

/// a subinterval in the Gauss Kronrod integreation process. 
class Subinterval {
public:
    myFloat a;              ///< the left endpoint 
    myFloat b;              ///< the right endpoint 
    myFloat r;              ///< the approximation to the integral of f 
    myFloat e;              ///< the estimated error of the integration 
    myFloat rabs;           ///< the appoximation to the integral of abs(f) 
    myFloat defabs;         ///< the approximation of the integral of abs(f-mean(f)) 
    int level;             ///< the level of the subinterval.  Initial = 0 
    int ndin;              ///< ndin = 1 for subintervals that are selected for initial valuation 

    /// default subinterval with everything set to 0 
    Subinterval() : a(0), b(0), r(0), e(0), ndin(0){};

    /// test whether the subinterval can be subdivided 
    bool is_divisible() const{
      myFloat& epmach = Machine::epsilon;
      myFloat& uflow = Machine::min;
      return e!=0 && (std::max<myFloat>(fabs(a), fabs(b)) >
                       (1 + 100 * epmach) * (fabs((a+b)/2) + 1000 * uflow));
    }
    /// calculate the approximations to the integrals over the subinterval 
    ///
    /// computes i = integral of f over (a,b), with error estimate
    ///          j = integral of abs(f) over (a,b)
    void integrate(Integrand f,                       ///< the function to be integrated
                   void* ex,                          ///< a pointer to the data used by f
                   const int n_gauss,                 ///< the number of Gauss nodes
                   const vector<myFloat>& w_gauss,    ///< the n_gauss weights for the Gauss integration
                   const vector<myFloat>& x_kronrod,  ///< the 2n_gauss+1 abscissae for the Kronrod integration.
                                                      ///< x_kronrod(1), x_kronrod(3), ...  abscissae of the n_gauss-point
                                                      ///< gauss rule.
                                                      ///< x_kronrod(0), x_kronrod(2), ...  abscissae which are optimally
                                                      ///< added to the n_gauss-point gauss rule.
                   const vector<myFloat>& w_kronrod  ///< the 2n_gauss+1 weights for the Kronrod integration
                   );
};

/// the data and comparison operator used to rank subintervals for subdivision 
class SubintervalComp {
public:
    int level_max;  ///< the cap on the level 

    /// operator ordering subintervals for subdivision 
    bool operator() (const Subinterval lhs, const Subinterval rhs) const {
        if ((!lhs.is_divisible()) && rhs.is_divisible()) return true;
        else if (lhs.is_divisible() && !rhs.is_divisible()) return false;
        else if (lhs.level == level_max && rhs.level < level_max) return true;
        else if (lhs.level < level_max && rhs.level == level_max) return false;
        else if (lhs.ndin < rhs.ndin) return true;
        else if (lhs.ndin > rhs.ndin) return false;
        else if (lhs.e < rhs.e)return true;
        else return false;
    }
    /// constructor for the subinterval comparison
    SubintervalComp(const int level_max ///< [in] cap on the level
                   ) : level_max(level_max) {}
};

/// print out subintervals 
void print_subs(ostream & os,                    ///< [in,out] the output stream to use
                const vector<Subinterval>& subs, ///< [in] a reference to the vector of subintervals 
                const int last,                  ///< [in] the number of subintervals to print 
                const vector<myFloat>& points     ///< [in]the original endpoints before the subdivision process 
);
  
/// print out a summary of the subintervals
void print_subs_summary(ostream & os,            ///< [in,out] the output stream to use
                const vector<Subinterval>& subs, ///< [in] a reference to the vector of subintervals
                const int last,                  ///< [in] the number of subintervals to print
                const vector<myFloat>& points    ///< [in] the original endpoints before the subdivision process
);
  
/// integration controller for adaptive Gauss Kronrod integration
/// adaptively integrate a function over a finite interval given an initial subdivision.
///
/// the routine calculates an approximation result to a given definite integral
/// i = integral of f over (a,b), hopefully satisfying following claim for accuracy
/// abs(i-result) < max(epsabs,epsrel*abs(i)). break points of the integration
/// interval, where local difficulties of the integrand may occur(e.g. singularities,
/// discontinuities),provided by user.
///
/// @author Joseph Dunn based on the original Fortran code of
/// @author piessens,robert, appl. math. & progr. div. - k.u.leuven
/// @author de doncker,elise, appl. math & progr. div. - k.u.leuven
///
class IntegrationController {
private:
  bool  noext;               ///< flag indicating no extrapolation is to be used.
  int n_gauss;               ///< the number of Gauss nodes
  vector<myFloat> w_gauss;   ///< the n_gauss weights for the Gauss integration
  vector<myFloat> x_kronrod; ///< the 2n_gauss+1 abscissae for the Kronrod integration.
                             ///< x_kronrod(1), x_kronrod(3), ...  abscissae of the n_gauss-point
                             ///< gauss rule.
                             ///< x_kronrod(0), x_kronrod(2), ...  abscissae which are optimally
                             ///< added to the n_gauss-point gauss rule.

  vector<myFloat> w_kronrod; ///< the 2n_gauss+1 weights for the Kronrod integration
  int limit;                 ///< the maximum number of subintervals
  int verbose;               ///< flag indicating verbose output
public:
  /// constructor for integration controller.  Calls kronrod to determine weights and nodes
  IntegrationController(bool noext,       ///< flag indicating no extrapolation is to used
                        int n_gauss,      ///< the number of Gauss nodes
                        myFloat epsabs,   ///< the target absolute error
                        myFloat epsrel,   ///< the target relative error
                        int limit,        ///< the maximum number of subintervals
                        int verbose       ///< flag indicating verbose output
  );
  
  
  myFloat epsabs;           ///< the target absolute error
  myFloat epsrel;           ///< the target relative error
  vector<Subinterval> subs; ///< work area for the integrator

  /// the termination codes returned by the integrator
  enum TerminationCode {normal, maximum_subdivisions, roundoff_error, bad_integrand,
                        extrapolation_roundoff_error, divergent_integral, bad_input};

  /// integrate the function f over the subinterval using the adaptive integration
  void integrate(Integrand f,                        ///< [in] the function to integrate
                 void* ex,                           ///< [in] pointer to the environment of the function
                 const vector<myFloat>& points,      ///< [in] the initial subdivsioin of the interval
                 myFloat& result,                    ///< [out] the approximation to the integral of f
                 myFloat& abserr,                    ///< [out] the estimated abs error of the approximation
                 int& neval,                         ///< [out] the number of function evaluations needed
                 TerminationCode& termination_code,  ///< [out] an error indicator.
                 int& last                           ///< [out] the number of subintervals actually used
  );
  
  /// get the verbose indicator
  int get_verbose() {return verbose;}
  
  /// print out the control information
  friend ostream& operator<< (ostream& os, const IntegrationController& ctl);
};

  /// A class with information needed to integrate a partiuclar Integrand
  class Integral {
  private:
    Integrand* f;
    void *ex;                ///< a pointer to an instance of
    const vector<myFloat>& points;      ///< the initial subdivsioin of the interval
    
    static string msgs[];                     ///< error messages from the integration
    IntegrationController& controller;  ///< the controller to use
    int verbose;             ///< the level of diagnostic output
    
  public:
    
    /// construct the functor
    Integral(Integrand* f,                  ///< pointer to the integrand
             void* ex,                      ///< a pointer to the integrands environment
             const vector<myFloat>& points, ///< A reference to the initial subdivision
             IntegrationController& ctl,     ///< [in] pointer to the integration controller
             int verbose                    ///< the level of diagnostic output
            )
        : f(f), ex(ex), points(points), controller(ctl), verbose(verbose){}
    
    /// return a human readable error message
    string msg() {return msgs[termination_code];};
    
    /// return a pointer to the instance of StandardStableDistribution
    myFloat result;        ///< the approximation of integral
    myFloat abserr;        ///< the estimated absolute error of the approximation
    int neval;            ///< the number of function evaluations used
    IntegrationController::TerminationCode termination_code; ///< the error indicator returned from IntegrationController::integrate
    int last;             ///< the number of subintervals required
    
    /// return the approximation to the integral of f(g(x))
    myFloat operator() ();
  };

} //namespace adaptive_integration

#endif // adaptive_integration_h
