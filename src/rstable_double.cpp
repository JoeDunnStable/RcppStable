//
//  rstable_double.cpp
//  lib_stable_double
//
//  Created by Joseph Dunn on 2/29/16.
//  Copyright © 2016 Joseph Dunn. All rights reserved.
//

#include "stable_double.h"

double srstable1(double alpha, double beta, double u1, double u2)
{
// Description:
//	 Returns one random variates for standard stable DF
    
// Calculate uniform and exponential distributed random numbers:
    double theta = M_PI * (u1-1./2.);
    double w = -log(u2);
    
    double result;
    if (alpha == 1 ) {
        result = (1/M_PI_2)*((M_PI_2+beta*theta)*tan(theta)-beta*log(M_PI_2*w*cos(theta)/(M_PI_2+beta*theta)));
    } else {
        double b_tan_pa = beta*tan(M_PI_2*alpha);
        double theta0 = fmin(fmax(-M_PI_2, atan(b_tan_pa) / alpha), M_PI_2);
        double c = pow(1+pow(b_tan_pa,2),1/(2*alpha));
        double a_tht = alpha*(theta+theta0);
        double r = ( c*sin(a_tht)/pow(cos(theta),1/alpha) ) *
                   pow(cos(theta-a_tht)/w,(1-alpha)/alpha);
// Use Parametrization 0:
        result = r - b_tan_pa;
    }
    return result;
}

