//
//  stable_mode_double.cpp
//  lib_stable_double
//
//  Created by Joseph Dunn on 1/21/16.
//  Copyright © 2016 Joseph Dunn. All rights reserved.
//

#include "stable_double.h"
#include <boost/math/tools/toms748_solve.hpp>

double stable_mode(double alpha, double beta, double dbltol,
                   integr_ctl_double& ctl, double zeta_tol, int verbose)
{
    if(alpha * beta == 0){
        return 0.;
    }
    else {
        if (fabs(alpha-1.)<1e-7) alpha=1.;  // The routine blows up otherwise
        double upper, lower;
        double outer_bound;
        double zeta;
        double eps=1.e-15;
        if ((alpha<1) && (fabs(beta)==1) ) {
            zeta = -beta*tan(M_PI*alpha/2);
            if (alpha < .1) return zeta;
            else outer_bound= zeta *(1-eps);
        
        } else
            outer_bound = (beta>0) ? -.7 : .7;
 
        if (beta>0){
            lower=outer_bound;
            upper=0.;
        } else {
            lower=0;
            upper=outer_bound;
        }
        ddx_sdstable_double_solve ddx_s(0, alpha, beta, dbltol, ctl, zeta_tol, verbose);
        rel_eps_tolerance_double tol(dbltol);
        boost::uintmax_t max_iter=1000;
        std::pair<double,double> mode=boost::math::tools::toms748_solve(ddx_s,lower,upper,tol,max_iter);
        return mode.first;
    }
}


