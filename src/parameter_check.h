//
///  \file parameter_check.h
///
///  \author Joseph Dunn on 7/27/16.
///  \copyright 2016 Joseph Dunn.
//

#ifndef parameter_check_h
#define parameter_check_h

#include <Eigen/Dense>
#include "stable_distribution.h"

namespace stable_distribution {
  
typedef Eigen::Matrix<myFloat, Eigen::Dynamic, 1> Vec;

void inline parameter_check(const string fcn, const Vec& x, const myFloat alpha, const myFloat beta,
                            const Vec& gamma, const Vec& delta, const int pm) {
  if (alpha <= 0 || alpha > 2) {
    throw std::out_of_range(fcn + ": alpha < 0 or alpha > 2");
  }
  if (beta < -1 || beta > 1) {
    throw std::out_of_range(fcn + ": beta < -1 or beta > 1");
  }
  if ((gamma.array() < 0).any()) {
    throw std::out_of_range(fcn + ": gamma < 0");
  }
  if (pm < 0 || pm > 2) {
    throw std::out_of_range(fcn + ": pm is not 0, 1 or 2");
  }
  if (x.size() != delta.size() || x.size() != delta.size()) {
    throw std::out_of_range(fcn + ": x, gamma and delta must be the same length");
  }
}

void inline parameter_switch(const myFloat alpha, const myFloat beta,
                      const Vec& gamma, const Vec& delta, const int pm,
                      IntegrationController& ctl, const int verbose,
                      Vec& gamma0, Vec& delta0) {
  myFloat& pi2 = StandardStableDistribution::pi2;
  if (pm == 0) {
    gamma0 = gamma;
    delta0 = delta;
  } else if (pm == 1) {
    gamma0=gamma;
    if (alpha != round(alpha))
      delta0 = delta + beta * gamma * tan(pi2*alpha);
    else if (alpha == 1)
      delta0 = delta.array() + beta * gamma.array() * gamma.array().log()/pi2;
    else
      delta0 = delta;
  } else { // pm=2
    gamma0 = pow(alpha, -1/alpha) * gamma;
    StandardStableDistribution std_stable_dist(alpha, beta, ctl, verbose);
    delta0 = delta - gamma0 * std_stable_dist.mode(sqrt(std::numeric_limits<myFloat>::epsilon())).first;
  }
}

} //namespace stable_distribution

#endif /* parameter_check_h */
