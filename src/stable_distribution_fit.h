//
/// \file stable_distribution_fit.h
/// \author Joseph Dunn
/// \copyright 2016 Joseph Dunn
/// \copyright Distributed under the terms of the GNU General Public License version 3

#ifndef stable_distribution_fit_h
#define stable_distribution_fit_h

#include "cubicspline.h"
#include <boost/math/distributions/students_t.hpp>
#include "stable_distribution_Vec.h"

namespace stable_distribution {
  
using adaptive_integration::IntegrationController;
using boost::math::students_t_distribution;

/// the data and functions necessary to calculate the quick version of pdf for large vectors x
///
/// uses pdf directly for points with very large absolute value and points near the mode
/// and uses cubic spline interpolation otherwise.

class DstableQuick {
  myFloat x_break_0;        ///< the lower bound for the lower cubic spline
  myFloat x_break_3;        ///< the upper bound of the region for the lower cubic spline
  myFloat x_break_4;        ///< the lower bound for the upper cubic spline
  myFloat x_break_7;        ///< the upper bound for the upper cubic spline
  StandardStableDistribution *std_stable_dist;   ///< pointer to the integration controller
  students_t_distribution<myFloat> dist_t;       ///< the standard student t distribution
  cubicspline spline_low;  ///< the lower cubic spline
  cubicspline spline_high; ///< the upper cubic spline
public:
  
  /// constract an instance of DstableQuick from a point to StandardStableDistribution
  DstableQuick(StandardStableDistribution *std_stable_dist ///< pointer to an instance of StandardStableDistribution
  );
  
  /// apply the spline or pdf to a vector of points x
  Vec operator() (const Vec& x ///< [in] vector of points at which to calculate the pdf
  );
  
  /// apply the spline or pdf to a single point x
  myFloat operator() (const myFloat& x ///< [in] the point at which to calculate the pdf
  );
  
  /// retrieve a vector of the knots for the lower spline
  Vec get_knots_low() {return spline_low.get_knots();}
  
  /// retrieve a vector of the knots for the upper spline
  Vec get_knots_high() {return spline_high.get_knots();}
  
  /// determine the region of the each of the points x
  uVec region(const Vec& x ///< [in] a vector of point at which to determine the region
  ) {
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
  /// return the student t cdf used as an intermediate step in calculation of spline
  Vec pt(const Vec& x ///< [in] a vector of points at which to calculate the cdf
  );
}; // class DstableQuick

/// return the pdf of the stable distribution using the quick approximation
Vec pdf_quick(const Vec& x,            ///< [in] the vector of points at which to calculate the pdf
                   const myFloat alpha,     ///< [in] the structural parameter of the distribution
                   const myFloat beta,      ///< [in] the skewness parameter of the distribution
                   const Vec& gamma,     ///< [in] the scale parameter of the distribution
                   const Vec& delta,     ///< [in] the location parameter of the distribution
                   const int pm,            ///< [in] the parameterization, 0, 1 or 2
                   const int log_flag,     ///< [in] return log of the pdf, otherwise the pdf
                   IntegrationController& ctl, ///< [in,out] reference to the integration controller
                   const int verbose       ///< [in] indicator for verbose output
                  );

/// return the pdf of the stable distribution for a vector of points y treating
/// all points with large absolute value as single points +infinity or -infinity
myFloat capped_pdf(const Vec& y,            ///< [in] the vector of points at which to calculate the pdf
                        const myFloat alpha,     ///< [in] the structural parameter of the distribution
                        const myFloat beta,      ///< [in] the skewness parameter of the distribution
                        const myFloat gamma,     ///< [in] the scale parameter of the distribution
                        const myFloat delta,     ///< [in] the location parameter of the distribution
                        const bool quick,       ///< [in] use the quick version of pdf
                        IntegrationController &ctl, ///< [in,out] reference to the integration controller
                        const int verbose       ///< [in] indicator for verbose output
                       );

/// return quantiles of a sample x sorting x in process
Vec quantile(Vec &x,         ///< [in,out] the sample points
             const Vec probs ///< [in] the requested probabilities
            );

/// the result of a run of stable fit
class FitResult {
public:
  string method;            ///< the method used: McCulloch, Dunn, mle or qmle
  myFloat alpha;             ///< the structural parameter of the stable dist
  myFloat beta;              ///< the skewness parameter of the stable distribution
  myFloat gamma;             ///< the scale parameter of the stable distribution
  myFloat delta;             ///< the location parameter of the distribution
  myFloat two_ll_n;          ///< two time the average log likelihood
  int n;                    ///< the number of sample points in the fit
  myFloat q_kurt;            ///< (q(.95)-q(.05)/(q(.75)-q(.25)
  myFloat q_skew;            ///< (q(.95)-2q(.5)+q(.05))/(q(.95)-q(.05)
  myFloat q_scale;           ///< (q(.75)-q(.25(
  myFloat q_location;        ///< q(.5)
  string convergence;       ///< did the algorithm converge
  unsigned int iterations;  ///< the number of iteration used
  myFloat cpu_time;          ///< the cpu time used

  /// constructor which simply copies the data and calculates the McCulloch's statistics
  FitResult(string method,            ///< the method used: McCulloch, Dunn, mle or qmle
            myFloat alpha,             ///< the structural parameter of the stable dist
            myFloat beta,              ///< the skewness parameter of the stable distribution
            myFloat gamma,             ///< the scale parameter of the stable distribution
            myFloat delta,             ///< the location parameter of the distribution
            myFloat two_ll_n,          ///< two time the average log likelihood
            int n,                    ///< the number of sample points in the fit
            Vec qs,                   ///< quantiles .05, .25, .5, .75, and .95
            string convergence,       ///< did the algorithm converge
            unsigned int iterations,  ///< the number of iteration used
            myFloat cpu_time           ///< the cpu time used
  )
    : method(method), alpha(alpha), beta(beta), gamma(gamma), delta(delta),
    two_ll_n(two_ll_n), n(n), convergence(convergence), iterations(iterations),
    cpu_time(cpu_time) {
    q_kurt=(qs(4)-qs(0))/(qs(3)-qs(1));
    q_skew=(qs(4)+qs(0)-2*qs(2))/(qs(4)-qs(0));
    q_scale=qs(3)-qs(1);
    q_location=qs(2);
  }
  
  /// send a FitResult to an output stream
  friend ostream& operator<< (ostream &os,        ///< the output stream
                              const FitResult &fr ///< the fit result to be sent
                              );
}; // class FitResult

/// place a heading on os for display of fit results
void result_heading(ostream &os ///< reference to the ostream to use
                   );

/// return a vector of fit results for the sample points in y
std::vector<FitResult> stable_fit(const Vec& y,                 ///< [in] the sample points to fit
                                  IntegrationController &ctl,  ///< [in,out] reference to integration controller
                                  const myFloat dbltol = 1e-10, ///< [in] the tolerance to use in the fit
                                  const string type="q",       ///< [in] the method to use: q, mle or q_mle
                                  const bool quick=false,      ///< [in] use the quick version of pdf
                                  const int verbose=0          ///< [in] indicator for verbose output
                                 );

} // namespace stable_distribution

#endif /* stable_distribution_fit_hpp */

