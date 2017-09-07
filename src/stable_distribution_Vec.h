//
/// \file  stable_distribution_Vec.h
/// \author Joseph Dunn
/// \copyright 2016 Joseph Dunn
/// \copyright Distributed under the terms of the GNU General Public License version 3

#ifndef stable_distribution_Vec_h
#define stable_distribution_Vec_h

#include <vector>
#include <Eigen/Dense>
#include "stable_distribution.h"

namespace stable_distribution {

using std::vector;
using adaptive_integration::IntegrationController;
  
/// an Eigen vector containing myFloats
typedef Eigen::Matrix<myFloat, Eigen::Dynamic, 1> Vec;

/// calculate the pdf of the standard stable distribution for vector of x's
Vec std_pdf(const Vec& x,                     ///< [in] the vector of points at which to calculate the pdf
        const myFloat alpha,         ///< [in] the structural parameter of the distribution
        const myFloat beta,          ///< [in] the skewness parameter of the distribution
        const int log_flag,          ///< [in] return log of the pdf, otherwise the pdf
        IntegrationController& ctl,  ///< [in,out] reference to the integration controller
        const int verbose            ///< [in] indicator for verbose output
);

/// calculate the pdf of the stable distribution for vector of x's
Vec pdf(const Vec& x,                     ///< [in] the vector of points at which to calculate the pdf
             const myFloat alpha,         ///< [in] the structural parameter of the distribution
             const myFloat beta,          ///< [in] the skewness parameter of the distribution
             const Vec& gamma,            ///< [in] a vector of scale parameters for distribution
             const Vec& delta,            ///< [in] a vector of location parameters for distribution
             const int pm,                ///< [in] the parameterization, 0, 1 or 2
             const int log_flag,          ///< [in] return log of the pdf, otherwise the pdf
             IntegrationController& ctl,  ///< [in,out] reference to the integration controller
             const int verbose            ///< [in] indicator for verbose output
             );

/// return the cdf of the standard stable distribution for a vector of points z
Vec std_cdf(const Vec& z,                     ///< [in] a vector of points at which to calculate the distribution
        const myFloat alpha,         ///< [in] the structural parameter of the distribution
        const myFloat beta,          ///< [in] the skewness parameter of the distribution
        const int lower_tail,        ///< [in] return the lower tail, otherwise the upper tail
        const int log_p,             ///< [in] return the log(pdf), otherwise the pdf
        IntegrationController& ctl,  ///< [in,out] reference to the integration controller
        const int verbose            ///< [in] indicator for verbose output
);

/// return the cdf of the stable distribution for a vector of points z
Vec cdf(const Vec& z,                     ///< [in] a vector of points at which to calculate the distribution
             const myFloat alpha,         ///< [in] the structural parameter of the distribution
             const myFloat beta,          ///< [in] the skewness parameter of the distribution
             const Vec& gamma,            ///< [in] a vector of scale parameters for distribution
             const Vec& delta,            ///< [in] a vector of location parameters for distribution
             const int pm,                ///< [in] the parameterization, 0, 1 or 2
             const int lower_tail,        ///< [in] return the lower tail, otherwise the upper tail
             const int log_p,             ///< [in] return the log(pdf), otherwise the pdf 
             IntegrationController& ctl,  ///< [in,out] reference to the integration controller
             const int verbose            ///< [in] indicator for verbose output
             );

/// return the inverse of the cdf of the standard stable distribution for a vector of probabilities p
Vec std_quantile(const Vec& p,               ///< [in] vector of probabilities at which to calculate the inverse
             myFloat alpha,               ///< [in] the structural parameter of the distribution
             myFloat beta,                ///< [in] the skewness parameter of the distribution
             int lower_tail,              ///< [in] return the lower tail, otherwise the upper tail
             int log_p,                   ///< [in] p is the log(pdf), otherwise the pdf
             const myFloat dbltol,        ///< [in] the tolerance for the inverse
             IntegrationController& ctl,  ///< [in,out] reference to the integration controller
             const int verbose            ///< [in] indicator for verbose output
);

/// return the inverse of the cdf of the stable distribution for a vector of probabilities p
Vec quantile(const Vec& p,               ///< [in] vector of probabilities at which to calculate the inverse
             myFloat alpha,               ///< [in] the structural parameter of the distribution
             myFloat beta,                ///< [in] the skewness parameter of the distribution
             const Vec& gamma,            ///< [in] a vector of scale parameters for distribution
             const Vec& delta,            ///< [in] a vector of location parameters for distribution
             const int pm,                ///< [in] the parameterization, 0, 1 or 2
             int lower_tail,              ///< [in] return the lower tail, otherwise the upper tail
             int log_p,                   ///< [in] p is the log(pdf), otherwise the pdf
             const myFloat dbltol,        ///< [in] the tolerance for the inverse
             IntegrationController& ctl,  ///< [in,out] reference to the integration controller
             const int verbose            ///< [in] indicator for verbose output
             );

/// return the derivative wrt x of the pdf of the standard stable distribution for vector of x
Vec std_ddx_pdf(const Vec& x,                     ///< [in] the vector of points at which to calculate the pdf
            const myFloat alpha,         ///< [in] the structural parameter of the distribution
            const myFloat beta,          ///< [in] the skewness parameter of the distribution
            IntegrationController& ctl,  ///< [in,out] reference to the integration controller
            const int verbose            ///< [in] indicator for verbose output
);

/// return the derivative wrt x of the pdf of the stable distribution for vector of x
Vec ddx_pdf(const Vec& x,                     ///< [in] the vector of points at which to calculate the pdf
                 const myFloat alpha,         ///< [in] the structural parameter of the distribution
                 const myFloat beta,          ///< [in] the skewness parameter of the distribution
                 const Vec& gamma,            ///< [in] a vector of scale parameters for distribution
                 const Vec& delta,            ///< [in] a vector of location parameters for distribution
                 const int pm,                ///< [in] the parameterization, 0, 1 or 2
                 IntegrationController& ctl,  ///< [in,out] reference to the integration controller
                 const int verbose            ///< [in] indicator for verbose output
                 );

/// return a vector of deviates from the standard stable distribution
Vec std_random_stable(const myFloat alpha,    ///< [in] the structural parameter of the distribution
                  const myFloat beta,          ///< [in] the skewness parameter of the distribution
                  const Vec &u1,               ///< [in] a vector of the first std uniform deviates used
                  const Vec &u2                ///< [in] a vector of the second std uniform deviates used
);

/// return a vector of deviates from the stable distribution
Vec random_stable(const myFloat alpha,    ///< [in] the structural parameter of the distribution
             const myFloat beta,          ///< [in] the skewness parameter of the distribution
             const Vec& gamma,            ///< [in] a vector of scale parameters for distribution
             const Vec& delta,            ///< [in] a vector of location parameters for distribution
             const int pm,                ///< [in] the parameterization, 0, 1 or 2
             IntegrationController& ctl,  ///< [in,out] reference to the integration controller
             const int verbose,           ///< [in] indicator for verbose output
             const Vec &u1,               ///< [in] a vector of the first std uniform deviates used
             const Vec &u2                ///< [in] a vector of the second std uniform deviates used
             );

} //namespace stable_distribution
  
#endif // stable_distribution_Vec_h
