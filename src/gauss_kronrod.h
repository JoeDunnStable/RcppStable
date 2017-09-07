///
/// \file  gauss_kronrod.h
/// \author Joseph Dunn
/// \copyright 2016 Joseph Dunn
/// \copyright Distributed under the terms of the GNU General Public License version 3

#ifndef gauss_kronrod_h
#define gauss_kronrod_h
#include "myFloat.h"
#include <vector>

namespace gauss_kronrod {
using std::vector;

/// Calculate integration nodes and weights given recursion coeffecients of orthogonal polynomials.
///  @see http://people.sc.fsu.edu/~jburkardt/f_src/toms726/toms726.html for the fortran version
///
void toms726(const int n_gauss,               ///< [in] the number of nodes
             const vector<myFloat>& a,  ///< [in] x-a is the coefficient of pj in the recursion
             const vector<myFloat>& b,  ///< [in] b is the coefficient of pj-1 in the recursion
             vector<myFloat>& x,        ///< [out] the nodes of the quadrature formula
             vector<myFloat>& w,        ///< [out] the weights attached to the nodes
             const int verbose          ///< [in] a flag indicating verbose output
);

/// Generate recurrence coefficients for monic Jacobi polynomials on [0,1].
///
///   r_jacobi01(n,a,b,c,d) generates the first n recurrence
///   coefficients for monic Jacobi polynomials on [0,1] with
///   parameters a and b. These are orthogonal on [0,1] relative
///   to the weight function w(t)=(1-t)^a t^b. The n alpha-
///   coefficients are stored in c, the n beta-
///   coefficients in d.
///
///   Translated to C++ by Joseph Dunn, Dec. 26, 2015
///
void r_jacobi01(const int n_gauss,         ///< [in] the number of recurrence coeficients to generate
                const myFloat a,     ///< [in] the exponent of (1-t) in the weight function
                const myFloat b,     ///<[in] the exponented of t in the weight function
                vector<myFloat>& c,  ///< [out] the alpha coefficients
                vector<myFloat>& d   ///< [out] the beta coefficients
);

/// Generate recurrence coefficients for monic Jacobi polynomials on [-1,1]
///
///   r_jacobi(n,a,b,c,d) generates the first n recurrence
///   coefficients for monic Jacobi polynomials with parameters
///   a and b. These are orthogonal on [-1,1] relative to the
///   weight function w(t)=(1-t)^a(1+t)^b. The n alpha-coefficients
///   are stored in c, the n beta-coefficients in d
///
///   Translated to C++ by Joseph Dunn, Dec. 26, 2015
///   Supplied by Dirk Laurie, 6-22-1998; edited by Walter
///   Gautschi, 4-4-2002.
///

void r_jacobi(const int n_gauss,        ///< [in] the number of recurrence coeficients to generate
              const myFloat a,    ///< [in] the exponent of (1-t) in the weight function
              const myFloat b,    ///<[in] the exponented of t in the weight function
              vector<myFloat>& c, ///< [out] the alpha coefficients
              vector<myFloat>& d  ///< [out] the beta coefficients
);

/// Generate the recurrence coeficients for the Jacobi-Kronrod matrix.
///
///   r_kronrod(n,a0,b0,a,b) produces the alpha- and beta-elements in
///   the Jacobi-Kronrod matrix of order 2n+1 for the weight
///   function (or measure) w. The input data for the weight
///   function w are the recurrence coefficients of the associated
///   orthogonal polynomials, which are stored in the vectors a0 & b0 of
///   dimension [ceiling(3*n/2)+1]. The alpha-elements are stored in
///   a, the beta-elements in b, of
///   dimension (2*n+1)
///
///   Supplied by Dirk Laurie, 6-22.1998
///   Edited by Pavel Holoborodko, n_gaussovember 7, 2011:
///
///   Translated to C++ by Joseph Dunn, Dec. 26 2015
///
///@see http://dip.sun.ac.za/~laurie/papers/kronrod/kronrod.ps

void r_kronrod(int n_gauss,                      ///< [in] the number of nodess in the Jacobi Gauss quadrature formula
               const vector<myFloat>& a0,  ///< [in] the alpha coeficients of the Jacobi Gauss recurrence
               const vector<myFloat>& b0,  ///< [in] the beta coefficients of the Jacobi Gauss recurrence
               vector<myFloat>& a,         ///< [out] the alpha coefficients of the Jacobi-Kronrod recurrence
               vector<myFloat>& b         ///< [out] the beta coefficients of the Jacobi-Kronrod recurrence
);
} // namespace gauss_kronrod

#endif // gauss_kronrod_h
