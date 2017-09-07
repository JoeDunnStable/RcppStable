///
/// @file stable_distribution_random.cpp
/// \author Joseph Dunn
/// \copyright 2016 Joseph Dunn
/// \copyright Distributed under the terms of the GNU General Public License version 3

#include "stable_distribution.h"

namespace stable_distribution {
  
  myFloat random_stable(myFloat alpha, myFloat beta, myFloat u1, myFloat u2)
  {
    // Description:
    //	 Returns one random variates for standard stable DF
    if (!StandardStableDistribution::initialized)
      StandardStableDistribution::initialize();
    auto &pi =  StandardStableDistribution::pi;
    auto &pi2 =  StandardStableDistribution::pi2;
    
    // Calculate uniform and exponential distributed random numbers:
    myFloat theta = pi * (u1-1./2.);
    myFloat w = -log(u2);
    
    myFloat result;
    if (alpha == 1 ) {
      result = (1/pi2)*((pi2+beta*theta)*tan(theta)-beta*log(pi2*w*cos(theta)/(pi2+beta*theta)));
    } else if (alpha==0) {
      if (w > 1)
        return 0;
      else if (theta + beta * pi2 > 0)
        return std::numeric_limits<myFloat>::infinity();
      else
        return -std::numeric_limits<myFloat>::infinity();
    } else {
      myFloat b_tan_pa = beta*tan(pi2*alpha);
      myFloat theta0 = min(max<myFloat>(-pi2, atan(b_tan_pa) / alpha), pi2);
      myFloat c = pow(1+pow(b_tan_pa,2),1/(2*alpha));
      myFloat a_tht = alpha*(theta+theta0);
      myFloat r = c*sin(a_tht)*exp(-log(cos(theta))/alpha+((1-alpha)/alpha)*log(cos(theta-a_tht)/w));
      // Use Parametrization 0:
      result = r - b_tan_pa;
    }
    return result;
  }
} // namespace stable_distribution

