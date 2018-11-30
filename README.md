RcppStable
==========

Introduction
------------
This R package contains density, probability and quantile functions,
and random number generation for (skew) stable distributions, using the
parametrizations of Nolan. Except for some top level interface routines, the
calculation is entirely in C++ using routines from the accompanying 
stable_distribution package. The skew stable distribution is a
generalization of the normal distribution appropriate for distributions with
power law tails.  The package started out as a port of the stabledist package
from R, but the algorithms have been extensively modified to improve
numerical stability. 

Motivation
----------
Originally the goal was to improve the speed of the stabledist package
by using C++.  However, the package has evolved to the point that
it's not only faster, but it's more reliably accurate.

Installation
------------
The package can be installed by downloading to a directory and then 
using the standard R procedure for installing packages.  

Documentation
-------------
Standard R documentation is included in the package and is available 
in the man sudirectory.

Acknowledgements
----------------
This package and the accoumpanying stable_distribution package 
build on the work of a number of others.

1. John Nolan developed the integral representation of the skew stable
distribution that's used in this package
[density.pdf](http://fs2.american.edu/jpnolan/www/stable/density.pdf).
He has released a compiled version 3.12.02 of his Fortran 
program STABLE,
[stable.exe](http://academic2.american.edu/~jpnolan/stable/stable.exe).
It's not open source and he's not explicit about its license. 
2. Many of the routines in this package started out life as C++ translations of
the stabledist package for R, the source code for which is available at 
[stabledist_0.7-0](https://cran.r-project.org/src/contrib/stabledist_0.7-0.tar.gz)
3. The routines in the adaptive_integration routine started out life
as machine C++ translations of Fortran routines in QUADPACK, which is part of 
SLATEC and therefore in the public domain (http://en.wikipedia.org/wiki/QUADPACK).
The routines were then heavily modified to take advantage of the C++ language.
4. One of the modifications made to QUADPACK is the addition of the ability to calculate
the nodes and weights for the Gauss Kronrod integration on the fly.  For this purpose
Dirk Laurie's method is used [kronrod.ps](http://dip.sun.ac.za/~laurie/papers/kronrod/kronrod.ps).
The routines used here are C++ translations of Dirk Laurie's MATLAB code, which
is included in Walter Gautschi's OPQ suite [OPQ](https://www.cs.purdue.edu/archives/2002/wxg/codes/OPQ.html).
5. Some of the Boost header files are used.  They're available at http://www.boost.org.
6. The Eigen 3 package of headers is used.  It's available at http://eigen.tuxfamily.org.
7. Pavel Holoborodko's Mpfr C++ header is used to wrap the GNU mpfr library for the
mpreal float type.
The header needs to be modified to allow it to work with Boost by commenting out the 
line #define MPREAL_HAVE_DYNAMIC_STD_NUMERIC_LIMITS. 
The original is available at http://www.holoborodko.com/pavel/mpfr/.
8. The GNU mpfr package is used for the multiprecision version.  It's available at http://www.mpfr.org/mpfr-current/#download.
9. Three components, meta.h, Problem.h and NelderMeadSolver.h, of Patrick Wieschollek's CppNumericalSolver package are used 
by the stable_fit routine and are included in the present package.  The original
is available at https://github.com/PatWie/CppNumericalSolvers.
10. The documentation is prepared using Doxygen, which is available at http://doxygen.org, and 
mathjax



