#include <limits>
#include <cmath>
#include <array>
#include "dqagp_v2.h"
#include <stdio.h>
#include <iostream>
#include <iomanip>

const double subinterval::epmach = std::numeric_limits<double>::epsilon();
const double subinterval::uflow = std::numeric_limits<double>::min();
const double subinterval::wtg_level = 0;
/*
double eps_mat[50][50];
*/

void
dqelg_v2(
  int& n,
  std::array<double, 52>& epstab,
  double& result,
  double& abserr,
  std::array<double, 3>& res3la,
  int& nres)
{
  //***begin prologue  dqelg_v2
  //***refer to  dqagie,dqagoe,dqagpe_v2,dqagse
  //***revision date  830518   (yymmdd)
  //***keywords  epsilon algorithm, convergence acceleration,
  //             extrapolation
  //***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
  //           de doncker,elise,appl. math & progr. div. - k.u.leuven
  //***purpose  the routine determines the limit of a given sequence of
  //            approximations, by means of the epsilon algorithm of
  //            p.wynn. an estimate of the absolute error is also given.
  //            the condensed epsilon table is computed. only those
  //            elements needed for the computation of the next diagonal
  //            are preserved.
  //***description
  //
  //           epsilon algorithm
  //           standard fortran subroutine
  //           double precision version
  //
  //           parameters
  //              n      - integer
  //                       epstab(n) contains the new element in the
  //                       first column of the epsilon table.
  //
  //              epstab - double precision
  //                       vector of dimension 52 containing the elements
  //                       of the two lower diagonals of the triangular
  //                       epsilon table. the elements are numbered
  //                       starting at the right-hand corner of the
  //                       triangle.
  //
  //              result - double precision
  //                       resulting approximation to the integral
  //
  //              abserr - double precision
  //                       estimate of the absolute error computed from
  //                       result and the 3 previous results
  //
  //              res3la - double precision
  //                       vector of dimension 3 containing the last 3
  //                       results
  //
  //              nres   - integer
  //                       number of calls to the routine
  //                       (should be zero at first call)
  //
  //***end prologue  dqelg_v2
  //
  //           list of major variables
  //           -----------------------
  //
  //           e0     - the 4 elements on which the computation of a new
  //           e1       element in the epsilon table is based
  //           e2
  //           e3                 e0
  //                        e3    e1    new
  //                              e2
  //           newelm - number of elements to be computed in the new
  //                    diagonal
  //           error  - error = abs(e1-e0)+abs(e2-e1)+abs(new-e2)
  //           result - the element in the new diagonal with least value
  //                    of error
  //
  //           machine dependent constants
  //           ---------------------------
  //
  //           epmach is the largest relative spacing.
  //           oflow is the largest positive magnitude.
  //           limexp is the maximum number of elements the epsilon
  //           table can contain. if this number is reached, the upper
  //           diagonal of the epsilon table is deleted.
  //
  //***first executable statement  dqelg_v2
  double epmach = std::numeric_limits<double>::epsilon();
  double oflow = std::numeric_limits<double>::max();
  nres++;
  abserr = oflow;
  result = epstab.at(n-1);
  if (n < 3) {
    abserr = fmax(abserr, 0.5e+01 * epmach * fabs(result));
    return;
  }
  int limexp = 50;
  epstab.at(n+1) = epstab.at(n-1);
  int newelm = (n - 1) / 2;
  epstab.at(n-1) = oflow;
  int num = n;
  int k1 = n;
  for (int i=1; i<=newelm; i++) {
    int k2 = k1 - 1;
    int k3 = k1 - 2;
    double res = epstab.at(k1+1);
    double e0 = epstab.at(k3-1);
    double e1 = epstab.at(k2-1);
    double e2 = res;
    double e1abs = fabs(e1);
    double delta2 = e2 - e1;
    double err2 = fabs(delta2);
    double tol2 = fmax(fabs(e2), e1abs) * epmach;
    double delta3 = e1 - e0;
    double err3 = fabs(delta3);
    double tol3 = fmax(e1abs, fabs(e0)) * epmach;
    if (err2 <= tol2 && err3 <= tol3) {
      //
      //           if e0, e1 and e2 are equal to within machine
      //           accuracy, convergence is assumed.
      //           result = e2
      //           abserr = abs(e1-e0)+abs(e2-e1)
      //
      result = res;
      abserr = fmax(err2 + err3, 0.5e+01 * epmach * fabs(result));
      return;
    }
    double e3 = epstab.at(k1-1);
    epstab.at(k1-1) = e1;
    double delta1 = e1 - e3;
    double err1 = fabs(delta1);
    double tol1 = fmax(e1abs, fabs(e3)) * epmach;
    //
    //           if two elements are very close to each other, omit
    //           a part of the table by adjusting the value of n
    //
    if (err1 <= tol1 || err2 <= tol2 || err3 <= tol3) {
      n = i + i - 1;
      break;
    }
    double ss = 0.1e+01 / delta1 + 0.1e+01 / delta2 - 0.1e+01 / delta3;
    double epsinf = fabs(ss * e1);
    //
    //           test to detect irregular behaviour in the table, and
    //           eventually omit a part of the table adjusting the value
    //           of n.
    //
    if (epsinf <= 0.1e-03) {
      n = i + i - 1;
      break;
    }
    //
    //           compute a new element and eventually adjust
    //           the value of result.
    //
    res = e1 + 0.1e+01 / ss;
    epstab.at(k1-1) = res;
    k1 = k1 - 2;
    double error = err2 + fabs(res - e2) + err3;
    if (error <= abserr) {
      abserr = error;
      result = res;
    }
  }
  //
  //           shift the table.
  //
  if (n == limexp) {
    n = 2 * (limexp / 2) - 1;
  }
  int ib = 1;
  if ((num / 2) * 2 == num) {
    ib = 2;
  }
  int ie = newelm + 1;
  for (int i=1; i<=ie; i++) {
    int ib2 = ib + 2;
    epstab.at(ib-1) = epstab.at(ib2-1);
    ib = ib2;
  }
  if (num != n) {
    int indx = num - n + 1;
    for (int i=1; i<=n; i++) {
      epstab.at(i-1) = epstab.at(indx-1);
      indx++;
    }
/*

    for (int i=0; i<num-2; i++) {
      for (int j=0; j<=num/2; j++)
        if (j>i+1)
          eps_mat[i][j]=0;
        else
          eps_mat[i][j]=eps_mat[i+num-n][j];
    }
*/
  }

  if (nres >= 4) {
    abserr = fabs(result - res3la.at(2)) + fabs(result - res3la.at(1)) + fabs(result - res3la.at(0));
    res3la.at(0) = res3la.at(1);
    res3la.at(1) = res3la.at(2);
    res3la.at(2) = result;
  } else{
    res3la.at(nres-1) = result;
    abserr = oflow;
  }
  //
  //           compute error estimate
  //
  abserr = fmax(abserr, 0.5e+01 * epmach * fabs(result));
/*
  std::cout << "dqelg_v2: "
            << "result = " << std::setw(20) << std::setprecision(15) << result
            << ", abserr = " << std::setw(20) << std::setprecision(15) << abserr
            << std::endl;


  if (n==3) {
    for (int i=0; i<50; i++)
      for (int j=0; j<50; j++)
        eps_mat[i][j]=0.;
    eps_mat[0][1]=epstab[0];
    eps_mat[0][0]=epstab[1];
    eps_mat[1][0]=epstab[2];
  } else {
    int row = (n-2);  //base zero
    int col = 0;  //base zero
    for (int i=n-1; i>=0; i-=2, row--, col++) {
      eps_mat[row][col]=epstab[i];
    }
  }

  for (int row = 0; row <(n-1);row++){
    for (int col=0; col<(n+1)/2; col++){
      if (eps_mat[row][col]==0) break;
      std::cout << std::setw(11) << std::setprecision(5) << result-eps_mat[row][col];
    }
    std::cout << std::endl;
  }
*/
}

void
subinterval::dqk21_v2(
  integr_fn f,
  void* ex
  )
{
  static double wg[] = {0.066671344308688137593568809893332e0,
                        0.149451349150580593145776339657697e0,
                        0.219086362515982043995534934228163e0,
                        0.269266719309996355091226921569469e0,
                        0.295524224714752870173892994651338e0};
  static double xgk[] ={0.995657163025808080735527280689003e0,
                        0.973906528517171720077964012084452e0,
                        0.930157491355708226001207180059508e0,
                        0.865063366688984510732096688423493e0,
                        0.780817726586416897063717578345042e0,
                        0.679409568299024406234327365114874e0,
                        0.562757134668604683339000099272694e0,
                        0.433395394129247190799265943165784e0,
                        0.294392862701460198131126603103866e0,
                        0.148874338981631210884826001129720e0,
                        0.000000000000000000000000000000000e0};
  static double wgk[] ={0.011694638867371874278064396062192e0,
                        0.032558162307964727478818972459390e0,
                        0.054755896574351996031381300244580e0,
                        0.075039674810919952767043140916190e0,
                        0.093125454583697605535065465083366e0,
                        0.109387158802297641899210590325805e0,
                        0.123491976262065851077958109831074e0,
                        0.134709217311473325928054001771707e0,
                        0.142775938577060080797094273138717e0,
                        0.147739104901338491374841515972068e0,
                        0.149445554002916905664936468389821e0};

  //***begin prologue  dqk21_v2
  //***date written   800101   (yymmdd)
  //***revision date  830518   (yymmdd)
  //***category no.  h2a1a2
  //***keywords  21-point gauss-kronrod rules
  //***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
  //           de doncker,elise,appl. math. & progr. div. - k.u.leuven
  //***purpose  to compute i = integral of f over (a,b), with error
  //                           estimate
  //                       j = integral of abs(f) over (a,b)
  //***description
  //
  //           integration rules
  //           standard fortran subroutine
  //           double precision version
  //
  //           parameters
  //            on entry
  //              f      - double precision
  //                       function subprogram defining the integrand
  //                       function f(x). the actual name for f needs to be
  //                       declared e x t e r n a l in the driver program.
  //
  //              a      - double precision
  //                       lower limit of integration
  //
  //              b      - double precision
  //                       upper limit of integration
  //
  //            on return
  //              result - double precision
  //                       approximation to the integral i
  //                       result is computed by applying the 21-point
  //                       kronrod rule (resk) obtained by optimal addition
  //                       of abscissae to the 10-point gauss rule (resg).
  //
  //              abserr - double precision
  //                       estimate of the modulus of the absolute error,
  //                       which should not exceed abs(i-result)
  //
  //              resabs - double precision
  //                       approximation to the integral j
  //
  //              resasc - double precision
  //                       approximation to the integral of abs(f-i/(b-a))
  //                       over (a,b)
  //
  //***references  (none)
  //***end prologue  dqk21_v2
  //
  //           the abscissae and weights are given for the interval (-1,1).
  //           because of symmetry only the positive abscissae and their
  //           corresponding weights are given.
  //
  //           xgk    - abscissae of the 21-point kronrod rule
  //                    xgk(2), xgk(4), ...  abscissae of the 10-point
  //                    gauss rule
  //                    xgk(1), xgk(3), ...  abscissae which are optimally
  //                    added to the 10-point gauss rule
  //
  //           wgk    - weights of the 21-point kronrod rule
  //
  //           wg     - weights of the 10-point gauss rule
  //
  // gauss quadrature weights and kronron quadrature abscissae and weights
  // as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
  // bell labs, nov. 1981.
  //
  //           list of major variables
  //           -----------------------
  //
  //           centr  - mid point of the interval
  //           hlgth  - half-length of the interval
  //           absc   - abscissa
  //           fval*  - function value
  //           resg   - result of the 10-point gauss formula
  //           resk   - result of the 21-point kronrod formula
  //           reskh  - approximation to the mean value of f over (a,b),
  //                    i.e. to i/(b-a)
  //
  //           machine dependent constants
  //           ---------------------------
  //
  //           epmach is the largest relative spacing.
  //           uflow is the smallest positive magnitude.
  //
  //***first executable statement  dqk21_v2
  double epmach = std::numeric_limits<double>::epsilon();
  double uflow = std::numeric_limits<double>::min();
  //
  double centr = 0.5e+00 * (a + b);
  double hlgth = 0.5e+00 * (b - a);
  double dhlgth = fabs(hlgth);
  //
  //           compute the 21-point kronrod approximation to
  //           the integral, and estimate the absolute error.
  //
  double resg = 0.0e+00;
  double fc = centr;
  f(&fc, 1, ex);
  double resk = wgk[10] * fc;
  rabs = fabs(resk);
  int j = 0;
  int jtw = 0;
  double absc = 0.;
  double fval1 = 0.;
  double fval2 = 0.;
  double fv1[10];
  double fv2[10];
  int jtwm1 = 0;
  for (j=1; j<=5; j++) {
    jtw = 2 * j;
    absc = hlgth * xgk[jtw-1];
    fv1[jtw-1] = centr-absc;
    fv2[jtw-1] = centr+absc;
    jtwm1 = 2 * j - 1;
    absc = hlgth * xgk[jtwm1-1];
    fv1[jtwm1-1] = centr-absc;
    fv2[jtwm1-1] = centr+absc;
  }
  f(fv1, 10, ex);
  f(fv2, 10, ex);
  double fsum = 0.;
  for (j=1; j<=5; j++) {
    jtw = 2 * j;
    fval1 = fv1[jtw-1];
    fval2 = fv2[jtw-1];
    fsum = fval1+fval2;
    resg += wg[j-1] * fsum;
    resk += wgk[jtw-1] * fsum;
    rabs += wgk[jtw-1] * (fabs(fval1) + fabs(fval2));
  }
  for (j=1; j<=5; j++) {
    jtwm1 = 2 * j - 1;
    fval1 = fv1[jtwm1-1];
    fval2 = fv2[jtwm1-1];
    fsum = fval1 + fval2;
    resk += wgk[jtwm1-1] * fsum;
    rabs += wgk[jtwm1-1] * (fabs(fval1) + fabs(fval2));
  }
  double reskh = resk * 0.5e+00;
  defabs = wgk[10] * fabs(fc - reskh);
  for (j=1; j<=10; j++) {
    defabs += wgk[j-1] * (fabs(fv1[j-1] - reskh) + fabs(fv2[j-1] - reskh));
  }
  r = resk * hlgth;
  rabs = rabs * dhlgth;
  defabs = defabs * dhlgth;
  e = fabs((resk - resg) * hlgth);
  if (defabs != 0.0e+00 && e != 0.0e+00) {
    e = defabs * fmin(0.1e+01, pow((0.2e+03 * e / defabs),
      1.5e+00));
  }
  if (rabs > uflow / (0.5e+02 * epmach)) {
    e = fmax((epmach * 0.5e+02) * rabs, e);
  }
}

void
dqpsrt(
  int const& limit,
  int const& last,
  int& maxerr,
  double& ermax,
  std::array<subinterval, LIMIT>& subs,
  std::array<int, LIMIT>& iord,
  int& nrmax)
{
  //***begin prologue  dqpsrt
  //***refer to  dqage,dqagie,dqagpe,dqawse
  //***routines called  (none)
  //***revision date  810101   (yymmdd)
  //***keywords  sequential sorting
  //***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
  //           de doncker,elise,appl. math. & progr. div. - k.u.leuven
  //***purpose  this routine maintains the descending ordering in the
  //            list of the local error estimated resulting from the
  //            interval subdivision process. at each call two error
  //            estimates are inserted using the sequential search
  //            method, top-down for the largest error estimate and
  //            bottom-up for the smallest error estimate.
  //***description
  //
  //           ordering routine
  //           standard fortran subroutine
  //           double precision version
  //
  //           parameters (meaning at output)
  //              limit  - integer
  //                       maximum number of error estimates the list
  //                       can contain
  //
  //              last   - integer
  //                       number of error estimates currently in the list
  //
  //              maxerr - integer
  //                       maxerr points to the nrmax-th largest error
  //                       estimate currently in the list
  //
  //              ermax  - double precision
  //                       nrmax-th largest error estimate
  //                       ermax = elist(maxerr)
  //
  //              elist  - double precision
  //                       vector of dimension last containing
  //                       the error estimates
  //
  //              iord   - integer
  //                       vector of dimension last, the first k elements
  //                       of which contain pointers to the error
  //                       estimates, such that
  //                       elist(iord(1)),...,  elist(iord(k))
  //                       form a decreasing sequence, with
  //                       k = last if last.le.(limit/2+2), and
  //                       k = limit+1-last otherwise
  //
  //              nrmax  - integer
  //                       maxerr = iord(nrmax)
  //
  //***end prologue  dqpsrt
  //
  //           check whether the list contains more than
  //           two error estimates.
  //
  //***first executable statement  dqpsrt
  if (last <= 2) {
    iord.at(0) = 1;
    iord.at(1) = 2;
    maxerr = iord.at(nrmax-1);
    ermax = subs.at(maxerr-1).e;
    return;
  }
  double errmax = subs.at(maxerr-1).e;
  if (nrmax != 1) {
    int ido = nrmax - 1;
    int isucc;
    for (int i=1; i<=ido; i++) {
      isucc = iord.at(nrmax - 2);
      // ***jump out of do-loop
      if (errmax <= subs.at(isucc-1).e) {
        break;
      }
      //
      //           this part of the routine is only executed if, due to a
      //           difficult integrand, subdivision increased the error
      //           estimate. in the normal case the insert procedure should
      //           start after the nrmax-th largest error estimate.
      //
      iord.at(nrmax-1) = isucc;
      nrmax = nrmax - 1;
    }
  }
  //
  //           compute the number of elements in the list to be maintained
  //           in descending order. this number depends on the number of
  //           subdivisions still allowed.
  //
  int jupbn = last;
  if (last > (limit / 2 + 2)) {
    jupbn = limit + 3 - last;
  }
  double errmin = subs.at(last-1).e;
  //
  //           insert errmax by traversing the list top-down,
  //           starting comparison from the element elist(iord(nrmax+1)).
  //
  int jbnd = jupbn - 1;
  int ibeg = nrmax + 1;
  bool insert_both_at_end=true;
  int i;
  if (ibeg <= jbnd) {
    int isucc;
    for (i=ibeg; i<=jbnd; i++) {
      isucc = iord.at(i-1);
      // ***jump out of do-loop
      if (errmax >= subs.at(isucc-1).e) {
        insert_both_at_end=false;
        break;
      }
      iord.at(i - 2) = isucc;
    }
  }
  if (insert_both_at_end) {
    iord.at(jbnd-1) = maxerr;
    iord.at(jupbn-1) = last;
    maxerr = iord.at(nrmax-1);
    ermax = subs.at(maxerr-1).e;
    return;
  }
  //
  //           insert errmin by traversing the list bottom-up.
  //
  iord.at(i-2) = maxerr;
  int k = jbnd;
  bool insert_one_at_i=true;
  for (int j=i; j<=jbnd; j++) {
    int isucc = iord.at(k-1);
    // ***jump out of do-loop
    if (errmin < subs.at(isucc-1).e) {
      insert_one_at_i=false;
      break;
    }
    iord.at(k) = isucc;
    k = k - 1;
  }
  if (insert_one_at_i){
    iord.at(i-1) = last;
  } else {
    iord.at(k) = last;
  }
  maxerr = iord.at(nrmax-1);
  ermax = subs.at(maxerr-1).e;
}

void
dqagpe_v2(
  integr_fn f,
  void* ex,
  double const& a,
  double const& b,
  int const& npts2,
  const std::array<double, 100>& points,
  double const& epsabs,
  double const& epsrel,
  int const& limit,
  double& result,
  double& abserr,
  int& neval,
  int& ier,
  std::array<subinterval, LIMIT>& subs,
  int& last)
{
  //***begin prologue  dqagpe_v2
  //***date written   800101   (yymmdd)
  //***revision date  830518   (yymmdd)
  //***category no.  h2a2a1
  //***keywords  automatic integrator, general-purpose,
  //             singularities at user specified points,
  //             extrapolation, globally adaptive.
  //***author  piessens,robert ,appl. math. & progr. div. - k.u.leuven
  //           de doncker,elise,appl. math. & progr. div. - k.u.leuven
  //***purpose  the routine calculates an approximation result to a given
  //            definite integral i = integral of f over (a,b), hopefully
  //            satisfying following claim for accuracy abs(i-result).le.
  //            max(epsabs,epsrel*abs(i)). break points of the integration
  //            interval, where local difficulties of the integrand may
  //            occur(e.g. singularities,discontinuities),provided by user.
  //***description
  //
  //        computation of a definite integral
  //        standard fortran subroutine
  //        double precision version
  //
  //        parameters
  //         on entry
  //            f      - double precision
  //                     function subprogram defining the integrand
  //                     function f(x). the actual name for f needs to be
  //                     declared e x t e r n a l in the driver program.
  //
  //            a      - double precision
  //                     lower limit of integration
  //
  //            b      - double precision
  //                     upper limit of integration
  //
  //            npts2  - integer
  //                     number equal to two more than the number of
  //                     user-supplied break points within the integration
  //                     range, npts2.ge.2.
  //                     if npts2.lt.2, the routine will end with ier = 6.
  //
  //            points - double precision
  //                     vector of dimension npts2, the first (npts2-2)
  //                     elements of which are the user provided break
  //                     points. if these points do not constitute an
  //                     ascending sequence there will be an automatic
  //                     sorting.
  //
  //            epsabs - double precision
  //                     absolute accuracy requested
  //            epsrel - double precision
  //                     relative accuracy requested
  //                     if  epsabs.le.0
  //                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
  //                     the routine will end with ier = 6.
  //
  //            limit  - integer
  //                     gives an upper bound on the number of subintervals
  //                     in the partition of (a,b), limit.ge.npts2
  //                     if limit.lt.npts2, the routine will end with
  //                     ier = 6.
  //
  //         on return
  //            result - double precision
  //                     approximation to the integral
  //
  //            abserr - double precision
  //                     estimate of the modulus of the absolute error,
  //                     which should equal or exceed abs(i-result)
  //
  //            neval  - integer
  //                     number of integrand evaluations
  //
  //            ier    - integer
  //                     ier = 0 normal and reliable termination of the
  //                             routine. it is assumed that the requested
  //                             accuracy has been achieved.
  //                     ier.gt.0 abnormal termination of the routine.
  //                             the estimates for integral and error are
  //                             less reliable. it is assumed that the
  //                             requested accuracy has not been achieved.
  //            error messages
  //                     ier = 1 maximum number of subdivisions allowed
  //                             has been achieved. one can allow more
  //                             subdivisions by increasing the value of
  //                             limit (and taking the according dimension
  //                             adjustments into account). however, if
  //                             this yields no improvement it is advised
  //                             to analyze the integrand in order to
  //                             determine the integration difficulties. if
  //                             the position of a local difficulty can be
  //                             determined (i.e. singularity,
  //                             discontinuity within the interval), it
  //                             should be supplied to the routine as an
  //                             element of the vector points. if necessary
  //                             an appropriate special-purpose integrator
  //                             must be used, which is designed for
  //                             handling the type of difficulty involved.
  //                         = 2 the occurrence of roundoff error is
  //                             detected, which prevents the requested
  //                             tolerance from being achieved.
  //                             the error may be under-estimated.
  //                         = 3 extremely bad integrand behaviour occurs
  //                             at some points of the integration
  //                             interval.
  //                         = 4 the algorithm does not converge.
  //                             roundoff error is detected in the
  //                             extrapolation table. it is presumed that
  //                             the requested tolerance cannot be
  //                             achieved, and that the returned result is
  //                             the best which can be obtained.
  //                         = 5 the integral is probably divergent, or
  //                             slowly convergent. it must be noted that
  //                             divergence can occur with any other value
  //                             of ier.gt.0.
  //                         = 6 the input is invalid because
  //                             npts2.lt.2 or
  //                             break points are specified outside
  //                             the integration range or
  //                             (epsabs.le.0 and
  //                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28))
  //                             or limit.lt.npts2.
  //                             result, abserr, neval, last, rlist(1),
  //                             and elist(1) are set to zero. alist(1) and
  //                             blist(1) are set to a and b respectively.
  //
  //            subs   - As of version 2 an array of subintervals of length
  //                     limit is used to contain the work.  The first last
  //                     subinterval contain the partition of the integration
  //                     range.  Eash subinterval has the following fields:
  //
  //                     a  - double precision
  //                          the left end point of the subinterval in the
  //                          partition of the given integration range (a,b)
  //
  //                     b  - double precision
  //                          the right end point of the subinterval in the
  //                          partition of the given integration range (a,b)
  //
  //                     r  - double precision
  //                          the integral approximation on the subinterval in the
  //                          partition of the given integration range (a,b)
  //
  //                     e  - double precision
  //                          the moduli of the absolute error on the subinterval in the
  //                          partition of the given integration range (a,b)
  //
  //                     ndin - integer
  //                          subintervals with ndin = 1 are handled first. the is based
  //                          on the variation in f over the subinterval on entry.
  //
  //            pts    - double precision
  //                     array of dimension at least npts2, containing the
  //                     integration limits and the break points of the
  //                     interval in ascending sequence.
  //
  //            last   - integer
  //                     number of subintervals actually produced in the
  //                     subdivisions process
  //
  //***references  (none)
  //***routines called  dqelg_v2,dqk21_v2,dqpsrt
  //***end prologue  dqagpe_v2
  //
  //            the dimension of rlist2 is determined by the value of
  //            limexp in subroutine epsalg (rlist2 should be of dimension
  //            (limexp+2) at least).
  //
  //            list of major variables
  //            -----------------------
  //
  //           rlist2    - array of dimension at least limexp+2
  //                       containing the part of the epsilon table which
  //                       is still needed for further computations
  //           maxerr    - pointer to the interval with largest error
  //                       estimate
  //           errmax    - elist(maxerr)
  //           erlast    - error on the interval currently subdivided
  //                       (before that subdivision has taken place)
  //           area      - sum of the integrals over the subintervals
  //           errsum    - sum of the errors over the subintervals
  //           errbnd    - requested accuracy max(epsabs,epsrel*
  //                       abs(result))
  //           *****1    - variable for the left subinterval
  //           *****2    - variable for the right subinterval
  //           last      - index for subdivision
  //           nres      - number of calls to the extrapolation routine
  //           numrl2    - number of elements in rlist2. if an appropriate
  //                       approximation to the compounded integral has
  //                       been obtained, it is put in rlist2(numrl2) after
  //                       numrl2 has been increased by one.
  //           erlarg    - sum of the errors over the intervals larger
  //                       than the smallest interval considered up to now
  //           extrap    - logical variable denoting that the routine
  //                       is attempting to perform extrapolation. i.e.
  //                       before subdividing the smallest interval we
  //                       try to decrease the value of erlarg.
  //           noext     - logical variable denoting that extrapolation is
  //                       no longer allowed (true-value)
  //
  //            machine dependent constants
  //            ---------------------------
  //
  //           epmach is the largest relative spacing.
  //           uflow is the smallest positive magnitude.
  //           oflow is the largest positive magnitude.
  //
  //***first executable statement  dqagpe_v2
/*
  for (int i =0; i<50; i++)
    for (int j=0; j<50; j++)
      eps_mat[i][j]=0;
*/
double const epmach = std::numeric_limits<double>::epsilon();
  double min_length = std::numeric_limits<double>::max();
  //
  //            test on validity of parameters
  //            -----------------------------
  //
  ier = 0;
  neval = 0;
  last = 0;
  result = 0.0e+00;
  abserr = 0.0e+00;
  int npts = npts2 - 2;
  if (npts2 < 2 || limit <= npts || limit > LIMIT || (epsabs <= 0.0e+00 &&
      epsrel < fmax(0.5e+02 * epmach, 0.5e-28))) {
    ier = 6;    // Invalid Input
    return;
  }
  //
  //            if any break points are provided, sort them into an
  //            ascending sequence.
  //
  double sign = (a>b) ? -1 : 1;
  std::array<double, 102> pts;
  pts.at(0) = fmin(a, b);
  if (npts != 0) {
    for (int i=1; i<=npts; i++) {
      pts.at(i) = points.at(i-1);
    }
  }
  pts.at(npts + 1) = fmax(a, b);
  int nint = npts + 1;
  std::sort(pts.begin(),pts.begin()+npts+2);
  if (pts.at(0) != fmin(a, b) || pts.at(npts2-1) != fmax(a, b)) {
//    std::cout << "Failed endpoint check" << std::endl
//              << "a = " << a << ", b = " << b << std::endl;
//    for (int i = 0; i<npts2; i++) std::cout << pts[i] << std::endl;
    ier = 6;          // Invalid input
      return;
  }

  //
  //            compute first integral and error approximations.
  //            ------------------------------------------------
  //
  double resabs = 0.0e+00;
  int i;
  for (i=1; i<=nint; i++) {
    subs.at(i-1).a = pts.at(i-1);
    subs.at(i-1).b = pts.at(i);
    subs.at(i-1).dqk21_v2(f, ex);
    abserr += subs.at(i-1).e;
    result += subs.at(i-1).r;
    if (subs.at(i-1).e == subs.at(i-1).defabs && subs.at(i-1).e != 0.0e+00) {
      subs.at(i-1).ndin = 1;
    } else {
      subs.at(i-1).ndin = 0;
    }
    resabs += subs.at(i-1).rabs;
    subs.at(i-1).level =0;
    min_length = fmin(min_length, fabs(subs.at(i-1).b-subs.at(i-1).a));
  }
  double errsum = 0.0e+00;
  for (i=1; i<=nint; i++) {
    errsum += subs.at(i-1).e;
  }
  //
  //           test on accuracy.
  //
  last = nint;
  neval = 21 * nint;
  double dres = fabs(result);
  double errbnd = fmax(epsabs, epsrel * dres);
  if (abserr <= 0.1e+03 * epmach * resabs && abserr > errbnd) {
    ier = 2;       // Roundoff error detected
  }
  std::make_heap(subs.begin(),subs.begin()+nint);
  if (limit < npts2) {
    ier = 1;         // Maximum number of subdivisions reached
  }
  if (ier != 0 || abserr <= errbnd) {
    if (ier > 2) {
      ier = ier - 1;
    }
    result = result * sign;
    return;
  }
  //
  //           initialization
  //           --------------
  //
  std::array<double, 52> rlist2;
  rlist2.fill(0.);
  rlist2.at(0) = result;
  std::array<double, 3> res3la;
  res3la.fill(0.);
  double correc = 0.;
  double area = result;
  int nres = 0;
  int numrl2 = 1;
  int ktmin = 0;
  bool extrap = false;
  bool noext = true;
  double erlarg = errsum;
  double ertest = errbnd;
  int level_max = 0;
  int iroff1 = 0;
  int iroff2 = 0;
  int iroff3 = 0;
  int ierro = 0;
  const double uflow = std::numeric_limits<double>::min();
  const double oflow = std::numeric_limits<double>::max();
  abserr = oflow;
  int ksgn = -1;
  if (dres >= (0.1e+01 - 0.5e+02 * epmach) * resabs) {
    ksgn = 1;
  }
  bool final_check, test_for_divergence, need_sum;
  //
  //           main do-loop
  //           ------------
  //
  for (last=npts2; last<=limit; last++) {
    //           At this point last is the number of intervals after the biscection.
    //
    //           bisect the subinterval with the largest score
    //           estimate.
    //
    std::pop_heap(subs.begin(), subs.begin()+last-1);
    int maxerr = last - 1;
    double errmax = subs.at(maxerr-1).e;
    double erlast = errmax;
    double rold = subs.at(maxerr-1).r;
    int level_cur = subs.at(maxerr-1).level+1;
    level_max = std::max(level_cur,level_max);
    double a1 = subs.at(maxerr-1).a;
    double b2 = subs.at(maxerr-1).b;
    double b1 = (a1+b2)/2;
    double a2 = b1;
    subs.at(maxerr-1).b = b1;
    subs.at(last-1).a = a2;
    subs.at(last-1).b = b2;
    subs.at(maxerr-1).dqk21_v2(f, ex);
    subs.at(last-1).dqk21_v2(f, ex);
    //
    //           improve previous approximations to integral
    //           and error and test for accuracy.
    //
    neval += 42;
    double area12 = subs.at(maxerr-1).r + subs.at(last-1).r;
    double erro12 = subs.at(maxerr-1).e + subs.at(last-1).e;
    errsum += erro12 - errmax;
    area += area12 - rold;
    if (subs.at(maxerr-1).defabs != subs.at(maxerr-1).e
                && subs.at(last-1).defabs != subs.at(last-1).e
                && erro12 > epsrel*fabs(area12)) {
      if ( fabs(rold - area12) <= 0.1e-04 * fabs(
          area12) && erro12 >= 0.99e+00 * errmax) {
        if (extrap) {
          iroff2++;
        }
        if (!extrap) {
          iroff1++;
        }
      }
      if (last > 10 && erro12 > errmax) {
        iroff3++;
      }
    }
    subs.at(maxerr-1).ndin = 0;
    subs.at(maxerr-1).level=level_cur;
    subs.at(last-1).ndin = 0;
    subs.at(last-1).level=level_cur;
    std::push_heap(subs.begin(),subs.begin()+last-1);
    std::push_heap(subs.begin(),subs.begin()+last);
    errbnd = fmax(epsabs, epsrel * fabs(area));
    min_length = fmin(min_length,fabs(b2-a1)/2);
    //
    //           test for roundoff error and eventually set error flag.
    //
/*
    std::cout << "a1 = " << a1 << ", b2 = " << b2
              << ", level = " << level_cur <<std::endl
              << "rold = " << rold << ", rnew = " << area12 << std::endl
              << "errorold = " << errmax << ", errornew = " << erro12 << std::endl
              << "iroff1 = " << iroff1
              << ", iroff2 = " << iroff2
              << ", iroff3 = " << iroff3 << std::endl;
*/
    if (iroff1 + iroff2 >= 10 || iroff3 >= 20) {
      ier = 2;        // Roundoff error detected
    }
    if (iroff2 >= 5) {
      ierro = 3;
    }
    //
    //           set error flag in the case that the number of
    //           subintervals equals limit.
    //
    if (last == limit) {
      ier = 1;   // Maximum number of subdivisions reached
    }
    //
    //           set error flag in the case of bad integrand behaviour
    //           at a point of the integration range
    //
    if (fmax(fabs(a1), fabs(b2)) <= (0.1e+01 +
        0.1e+03 * epmach) * (fabs(a2) + 0.1e+04 * uflow)) {
      ier = 4;   // Bad integrand behavior
    }
    final_check=true;
    if (errsum <= errbnd) {
      // ***jump out of do-loop
 //     std::cout << "dqagpe: Success errsum <= errbnd" << std::endl;
      need_sum = true;
      final_check = false;
      break;
    }
    if (ier != 0) {
      // ***jump out of do-loop
//      std::cout << "dqagpe: Aborting with raw error code = " << ier << std::endl;
      break;
    }
    erlarg = erlarg - erlast;
    if (extrap)
      erlarg += erro12;
    if (subs.front().level == level_max) {
      if (noext) {
        continue;
      }
      //
      //           perform extrapolation.
      //
      numrl2++;
      rlist2.at(numrl2-1) = area;
      if (numrl2 > 2) {
        double reseps, abseps;
        dqelg_v2(numrl2, rlist2, reseps, abseps, res3la, nres);
        ktmin++;
        if (ktmin > 10 && abserr < 0.1e-02 * errsum) {
          ier = 5;    // Roundoff error in extrapolation table
        }
        if (abseps < abserr) {
//          std::cout << "Using reseps = " << reseps
//                    << ", vs result = " << result << std::endl;
          ktmin = 0;
          abserr = abseps;
          result = reseps;
          correc = erlarg;
          ertest = fmax(epsabs, epsrel * fabs(reseps));
          // ***jump out of do-loop
          if (abserr < ertest) {
 //           std::cout << "dqagpe: Success.  abserr from dqelg < errtest." << std::endl;
            break;
          }
        }
        //
        //           prepare bisection of the smallest interval.
        //
        if (numrl2 == 1) {
//          std::cout << "Performed extrapolation but numrl2 == 1" << std::endl;
          noext = true;
        }
        if (ier >= 5) {
//          std::cout << "dqagpe: Aborting because of roundoff error in extrapolation table" << std::endl;
          break;
        }
      }
      extrap = false;
      erlarg = errsum;
    } else {
      extrap = true;
    }
  } // main do loop
  //
  //           set the final result.
  //           ---------------------
  //
  if (final_check){
    if (abserr == oflow) {
      test_for_divergence=false;
      need_sum=true;
    } else if ((ier + ierro) == 0) {
      test_for_divergence=true;
    } else {
      if (ierro == 3) {
//        std::cout << "dqagpe: ierro = 3 adding " << correc << " to abserr." << std::endl;
        abserr += correc;
      }
      if (ier == 0) {
        ier = 3;        // Roundoff error detected
      }
      if (result != 0.0e+00 && area != 0.0e+00) {
        if (abserr / fabs(result) > errsum / fabs(area)) {
          need_sum=true;
          test_for_divergence=false;
        } else {
          test_for_divergence=true;
        }
      } else if (abserr > errsum) {
         need_sum = true;
         test_for_divergence=false;
      } else if (area == 0.0e+00) {
        need_sum=false;
        test_for_divergence=false;
      } else {
        test_for_divergence=true;
      }
    }
    //
    //           test on divergence.
    //
    if (test_for_divergence){
      need_sum=false;
      if (ksgn != (-1) && fmax(fabs(result), fabs(
          area)) > resabs * 0.1e-01) {
        if (0.1e-01 > (result / area) || (result / area) > 0.1e+03 ||
            errsum > fabs(area)) {
          ier = 6;      // Divergent integral
        }
      }
    } // test for divergence
  } // final_check
  //
  //           compute global integral sum.
  //
  if (need_sum){
//    std::cout << "Calculating sum" << std::endl;
    result = 0.0e+00;
    for (int k=1; k<=last; k++) {
      result += subs.at(k-1).r;
    }
    abserr = errsum;
  }
  if (ier > 2) {
    ier = ier - 1;
  }
  result = result * sign;
}




