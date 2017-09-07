//
/// \file stable_distribution_Vec.cpp
/// \author Joseph Dunn
/// \copyright 2016 Joseph Dunn
/// \copyright Distributed under the terms of the GNU General Public License version 3

#include "stable_distribution_Vec.h"
#include "parameter_check.h"

namespace stable_distribution {
  
using namespace adaptive_integration;
  
Vec std_pdf(const Vec& x, const myFloat alpha, const myFloat beta, int log_flag,
            IntegrationController& ctl, const int verbose) {
  StandardStableDistribution std_stable_dist(alpha, beta, ctl, verbose);
  Vec ret(x.size());
  for (int i=0; i<x.size(); i++)
    ret(i) = std_stable_dist.pdf(x(i), log_flag);
  return ret;
}
  
Vec pdf(const Vec& x, const myFloat alpha, const myFloat beta,
        const Vec& gamma, const Vec& delta, const int pm, const int log_flag,
             IntegrationController& ctl, const int verbose){

  parameter_check("pdf", x, alpha, beta, gamma, delta, pm);

  // Switch to S0 parameterization
  Vec gamma0(gamma.size());
  Vec delta0(delta.size());
  parameter_switch(alpha, beta, gamma, delta, pm, ctl, verbose,
                   gamma0, delta0);
  
// Shift and Scale:
  Vec x_std = (x - delta0).array()/gamma0.array();
  
  Vec ret = std_pdf(x_std, alpha, beta, log_flag, ctl, verbose);
  if (log_flag)
    return ret.array() - gamma0.array().log();
  else
    return ret.array()/gamma0.array();
}

Vec std_cdf(const Vec& x, const myFloat alpha, const myFloat beta, int lower_tail, int log_flag,
            IntegrationController& ctl, const int verbose) {
  StandardStableDistribution std_stable_dist(alpha, beta, ctl, verbose);
  Vec ret(x.size());
  for (int i=0; i<x.size(); i++)
    ret(i) = std_stable_dist.cdf(x(i), lower_tail, log_flag);
  return ret;
}
  
Vec cdf(const Vec& x, const myFloat alpha, const myFloat beta,
        const Vec& gamma, const Vec& delta, const int pm, const int lower_tail, const int log_p,
             IntegrationController& ctl, const int verbose) {

  parameter_check("pdf", x, alpha, beta, gamma, delta, pm);
  
  // Switch to S0 parameterization
  Vec gamma0(gamma.size());
  Vec delta0(delta.size());
  parameter_switch(alpha, beta, gamma, delta, pm, ctl, verbose,
                   gamma0, delta0);
  
  // Shift and Scale:
  Vec x_std = (x - delta0).array()/gamma0.array();
  return std_cdf(x_std, alpha, beta, lower_tail, log_p, ctl, verbose);
}

Vec std_quantile(const Vec& p, const myFloat alpha, const myFloat beta,
             const int lower_tail, const int log_p,
             const myFloat dbltol, IntegrationController& ctl,
             const int verbose) {
  StandardStableDistribution std_stable_dist(alpha, beta, ctl, verbose);
  Vec ret(p.size());
  for (int i=0; i<p.size(); i++)
    ret(i)=std_stable_dist.quantile(p(i), lower_tail, log_p, dbltol);
  return ret;
}

Vec quantile(const Vec& p, const myFloat alpha, const myFloat beta,
              const Vec& gamma, const Vec& delta, const int pm, const int lower_tail, const int log_p,
             const myFloat dbltol, IntegrationController& ctl,
             const int verbose) {

  parameter_check("quantile", p, alpha, beta, gamma, delta, pm);
  
  // Switch to S0 parameterization
  Vec gamma0(gamma.size());
  Vec delta0(delta.size());
  parameter_switch(alpha, beta, gamma, delta, pm, ctl, verbose,
                   gamma0, delta0);
  Vec ret = std_quantile(p, alpha, beta, lower_tail, log_p, dbltol, ctl, verbose);
  return gamma0.array()*ret.array() + delta0.array();
}

Vec std_ddx_pdf(const Vec& x, const myFloat alpha, const myFloat beta,
            IntegrationController& ctl, const int verbose){
  Vec ret(x.size());
  StandardStableDistribution std_stable_dist(alpha, beta, ctl, verbose);
  for (int i=0; i<x.size(); i++)
    ret(i)=std_stable_dist.ddx_pdf(x(i));
  return ret;
}

Vec ddx_pdf(const Vec& x, const myFloat alpha, const myFloat beta,
            const Vec& gamma, const Vec& delta, const int pm,
                 IntegrationController& ctl, const int verbose){
  parameter_check("ddx_pdf", x, alpha, beta, gamma, delta, pm);
  
  // Switch to S0 parameterization
  Vec gamma0(gamma.size());
  Vec delta0(delta.size());
  parameter_switch(alpha, beta, gamma, delta, pm, ctl, verbose,
                   gamma0, delta0);
  
  // Shift and Scale:
  Vec x_std = (x - delta0).array()/gamma0.array();
  Vec ret = std_ddx_pdf(x_std, alpha, beta, ctl, verbose);
  return ret.array()/gamma0.array();
}

Vec std_random_stable(const myFloat alpha, const myFloat beta,
                  const Vec &u1, const Vec &u2) {
  size_t n = u1.size();
  Vec ret(n);
  for (size_t i = 0; i<n; ++i)
    ret(i)=random_stable(alpha, beta, u1(i), u2(i));
  return ret;
}
  
Vec random_stable(const myFloat alpha, const myFloat beta,
                  const Vec& gamma, const Vec& delta, const int pm,
                  IntegrationController& ctl, const int verbose,
                  const Vec &u1, const Vec &u2) {
  parameter_check("random_stable", u1, alpha, beta, gamma, delta, pm);
  if (u1.size() != u2.size()) {
    throw std::out_of_range("random_stable: u1.size() != u2.size()");
  }

  // Switch to S0 parameterization
  Vec gamma0(gamma.size());
  Vec delta0(delta.size());
  
  // We need a controller to calculate the mode
  parameter_switch(alpha, beta, gamma, delta, pm, ctl, verbose,
                   gamma0, delta0);
  return gamma0.array()* std_random_stable(alpha, beta, u1, u2).array() + delta0.array();
  
}
  
} //namespace stable_distribution

