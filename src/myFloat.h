//
///  @file myFloat.h
/// \author Joseph Dunn
/// \copyright 2016 Joseph Dunn
/// \copyright Distributed under the terms of the GNU General Public License version 3

#ifndef myFloat_h
#define myFloat_h

#undef BOOST_MATH_OVERFLOW_ERROR_POLICY
#define BOOST_MATH_OVERFLOW_ERROR_POLICY ignore_error
// #include <boost/math/constants/constants.hpp>
#include <boost/lexical_cast.hpp>

using boost::lexical_cast;

#include <RcppEigen.h>
#include <Eigen/Eigenvalues>

#ifdef CPP_BIN_FLOAT

#include <boost/multiprecision/cpp_bin_float.hpp>
typedef boost::multiprecision::cpp_bin_float_50 myFloat;

#elif defined(MPFR_FLOAT_50)

#include <boost/multiprecision/mpfr.hpp>

typedef boost::multiprecision::number<boost::multiprecision::mpfr_float_backend<50>, boost::multiprecision::et_off> myFloat;

#elif defined(MPREAL)

//#include <boost/math/bindings/mpreal.hpp>
#include <mpreal.h>
typedef mpfr::mpreal myFloat;
using mpfr::digamma;
using mpfr::const_pi;
using mpfr::const_euler;

#include <boost/math/tools/real_cast.hpp>
namespace mpfr {
  template <class Policy>
  inline long long lltrunc(mpfr::mpreal const& x, const Policy& pol)
  {
    return boost::math::tools::real_cast<long long>(boost::math::trunc(x, pol));
  }

}


#else

#include <cmath>
typedef double myFloat;
using std::isfinite;
using std::isnan;
using std::isinf;
#include <boost/math/special_functions/digamma.hpp>
using boost::math::digamma;
#include <boost/math/special_functions/zeta.hpp>
using boost::math::zeta;
#endif

#if defined(CPP_BIN_FLOAT) || defined(MPFR_FLOAT_50)

#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/zeta.hpp>
#include <boost/math/special_functions/erf.hpp>
//using fmax = boost::multiprecision::max;
//using fmin = boost::multiprecision::min;

using boost::math::isfinite;
using boost::math::isnan;
using boost::math::isinf;
using boost::math::expm1;
using boost::math::log1p;
using boost::math::erf;
using boost::math::erfc;
using boost::math::tgamma;
using boost::math::lgamma;
using boost::math::digamma;
using boost::math::zeta;

inline myFloat pow(myFloat x, myFloat y) {
  if (x<0)
    return NAN;
  else if (x == 0 && y == 0)
    return NAN;
  else if (x == 0 && y < 0)
    return std::numeric_limits<myFloat>::infinity();
  else if (x == 0 && y>0)
    return 0;
  else
    return boost::multiprecision::pow(x, y);
}

#endif

#if !defined(MPREAL)

#include <boost/math/constants/constants.hpp>
inline myFloat const_pi() {
  return boost::math::constants::pi<myFloat>();
}

inline myFloat const_euler() {
  return boost::math::constants::euler<myFloat>();
}

#endif

#if defined(CPP_BIN_FLOAT) || defined(MPREAL) || defined(MPFR_FLOAT_50)

namespace Eigen {

  template<> struct NumTraits<myFloat> : GenericNumTraits<myFloat>
  {
    static inline Real dummy_precision() { return epsilon()*static_cast<myFloat>(1024); }
  };

} // namespace eigen


#endif

#include <string>
using std::string;

/// Return the type of floating point numbers used
inline string floating_type() {
#if defined(CPP_BIN_FLOAT)
  return string("cpp_bin_float");
#elif defined(MPFR_FLOAT_50)
  return string("mpfr_float_50");
#elif defined(MPREAL)
  return string("mpreal");
#else
  return string("double");
#endif

}



#endif /* myFloat_h */

