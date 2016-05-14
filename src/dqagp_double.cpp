#include <limits>
#include <cmath>
#include <array>
#include <algorithm>
#include "dqagp_double.h"
#include "kronrod_double.h"
#include <stdio.h>
#include <iostream>
#include <iomanip>

#include <RcppEigen.h>
#define cout Rcpp::Rcout
#define cerr Rcpp::Rcerr
/*
using std::cout;
using std::cerr;
*/
using std::endl;
using std::setw;
using std::setprecision;
using std::right;
using std::array;
using std::sort;
using std::make_heap;
using std::pop_heap;
using std::push_heap;


template<typename T>
class sort_data {
public:
    T a;
    int index;
    sort_data(T a, int index) : a(a), index(index){}
};

template<typename T>
bool sd_lt(const sort_data<T> &lhs,const sort_data<T>& rhs) {
    return lhs.a < rhs.a;
}

integr_ctl_double::integr_ctl_double(const bool noext, const double wtg_level, const int N,
                       const double epsabs, const double epsrel, const int limit, const int verbose)
    : noext(noext), wtg_level(wtg_level), N(N), epsabs(epsabs), epsrel(epsrel),
      limit(limit), verbose(verbose) {
    if (epsabs==0 && epsrel < fmax(50*std::numeric_limits<double>::epsilon(),1.e-25)) {
        integr_ctl_double::epsrel = fmax(50*std::numeric_limits<double>::epsilon(),1.e-25);
        cerr << "integr_ctl_double: Resetting epsrel to minimum allowable: " << integr_ctl_double::epsrel << endl;
    }
    vector<double> x_gauss(N);  // We don't keep it because the nodes are in x_kronrod
    w_gauss.resize(N);
    x_kronrod.resize(2*N+1);
    w_kronrod.resize(2*N+1);
    int M = std::max(2*N,int(ceil(3.*N/2.)+1));
    vector<double> a0(M);
    vector<double> b0(M);
    r_jacobi01(M,double(0),double(0), a0, b0);
    toms726(N, a0, b0, x_gauss, w_gauss, verbose);
    vector<double> a(2*N+1);
    vector<double> b(2*N+1);
    r_kronrod(N, a0, b0, a, b);
    toms726(2*N+1, a, b, x_kronrod, w_kronrod, verbose);
    for (int i = 0; i<N; i++) {
        x_gauss.at(i)=2*x_gauss.at(i)-1;
        w_gauss.at(i)=2*w_gauss.at(i);
    }

    for (int i = 0; i<2*N+1; i++) {
        x_kronrod.at(i) = 2*x_kronrod.at(i)-1;
        w_kronrod.at(i) = 2*w_kronrod.at(i);
    }
}

ostream& operator<< (ostream& os, integr_ctl_double& ctl) {
    os << setw(15) << "noext = " << ctl.noext << endl
    << setw(15) << "wtg_level = " << ctl.wtg_level << endl
    << setw(15) << "N = " << ctl.N << endl
    << setw(15) << "epsabs = " << ctl.epsabs << endl
    << setw(15) << "epsrel = " << ctl.epsrel << endl
    << setw(15) << "limit = " << ctl.limit << endl;
    os << "Gauss nodes and weights: " << endl << endl;
    for (int i=0; i<ctl.N; i++) {
        os << setw(25) << setprecision(16) << ctl.x_kronrod.at(1+2*i)
        << setw(25) << setprecision(16) << ctl.w_gauss.at(i) << endl;
    }
    os << endl << "Kronrod nodes and weights:" << endl << endl;
    for (int i=0; i<2*ctl.N+1; i++)
        os << setw(25) << setprecision(16) << ctl.x_kronrod.at(i)
        << setw(25) << setprecision(16) << ctl.w_kronrod.at(i) << endl;
    return os;

}

bool subinterval_double::initialized = false;
double subinterval_double::epmach;
double subinterval_double::uflow;
integr_ctl_double* subinterval_double::ctl;

void subinterval_double::initialize(integr_ctl_double* ctl_in){
    ctl=ctl_in;
    epmach = std::numeric_limits<double>::epsilon();
    uflow = std::numeric_limits<double>::min();
    initialized=true;
};

void print_subs(const vector<subinterval_double> subs, const int last, const vector<double> points) {
    if (last==0) return;
    if (subs.size()<last) {
        throw std::range_error("print_subs: last is greater than the length of subs.");
    }
    // Determine the geometric order of the subintervals

    vector<sort_data<double> > srt_a;
    for (int i=0; i<last; i++) {
        sort_data<double> elem(subs.at(i).a,i);
        srt_a.push_back(elem);
    }
    sort(srt_a.begin(), srt_a.end(),sd_lt<double>);

    // Determine the error ranks for each subinterval

    vector<sort_data<double> > srt_eord;
    for (int i=0; i<last; i++) {
        sort_data<double> elem(subs.at(i).e,i);
        srt_eord.push_back(elem);
    }
    sort(srt_eord.begin(), srt_eord.end(),sd_lt<double>);

    vector<sort_data<int> > srt_iord;
    for (int i=0; i<last; i++) {
        sort_data<int> elem(srt_eord.at(i).index,last-1-i);
        srt_iord.push_back(elem);
    }
    sort(srt_iord.begin(), srt_iord.end(),sd_lt<int>);

    cout << " "
         << setw(13) << right << "a"
         << setw(13) << right << "b"
         << setw(13) << right << "length"
         << setw(13) << right << "r"
         << setw(13) << right << "average"
         << setw(13) << right << "e"
         << setw(5) << right << "rank"
         << setw(6) << right << "level" << endl;
    for (int i=0; i<last;i++){
        int j = srt_a[i].index;

        // Determine whether the left endpoint is in points.
        bool ispt = i==0;
        for (int ipt = 0; ipt<points.size(); ipt++) {
            ispt = ispt || subs.at(j).a==points.at(ipt);
            if (ispt) break;
        }
        if (ispt)
            cout << "*";
        else
            cout << " ";
        cout << setw(13) << setprecision(5) << subs.at(j).a
        << setw(13) << setprecision(5) << subs.at(j).b
        << setw(13) << setprecision(5) << subs.at(j).b-subs.at(j).a
        << setw(13) << setprecision(5) << subs.at(j).r
        << setw(13) << setprecision(5) << subs.at(j).r/(subs.at(j).b-subs.at(j).a)
        << setw(13) << setprecision(5) << subs.at(j).e
        << setw(5) << srt_iord.at(j).index+1
        << setw(6)  << subs.at(j).level << endl;
    }
}

/*
double eps_mat[50][50];
*/

void
dqelg(
  int& n,
  vector<double>& epstab,
  double& result,
  double& abserr,
  vector<double>& res3la,
  int& nres)
{
  //***begin prologue  dqelg_double
  //***refer to  dqagie,dqagoe,dqagpe_double,dqagse
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
  //***end prologue  dqelg_double
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
  //***first executable statement  dqelg_double

  double epmach = std::numeric_limits<double>::epsilon();
  double oflow = std::numeric_limits<double>::max();
  nres++;
  abserr = oflow;
  result = epstab.at(n-1);
  if (n < 3) {
    abserr = std::max(abserr, 0.5e+01 * epmach * fabs(result));
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
    double tol2 = std::max(fabs(e2), e1abs) * epmach;
    double delta3 = e1 - e0;
    double err3 = fabs(delta3);
    double tol3 = std::max(e1abs, fabs(e0)) * epmach;
    if (err2 <= tol2 && err3 <= tol3) {
      //
      //           if e0, e1 and e2 are equal to within machine
      //           accuracy, convergence is assumed.
      //           result = e2
      //           abserr = abs(e1-e0)+abs(e2-e1)
      //
      result = res;
      abserr = std::max(err2 + err3, 50 * epmach * fabs(result));
      return;
    }
    double e3 = epstab.at(k1-1);
    epstab.at(k1-1) = e1;
    double delta1 = e1 - e3;
    double err1 = fabs(delta1);
    double tol1 = std::max(e1abs, fabs(e3)) * epmach;
    //
    //           if two elements are very close to each other, omit
    //           a part of the table by adjusting the value of n
    //
    if (err1 <= tol1 || err2 <= tol2 || err3 <= tol3) {
      n = i + i - 1;
      break;
    }
    double ss = 1 / delta1 + 1 / delta2 - 1 / delta3;
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
    res = e1 + 1 / ss;
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
  abserr = std::max(abserr, 50 * epmach * fabs(result));
/*
  cout << "dqelg_double: "
            << "result = " << setw(20) << setprecision(15) << result
            << ", abserr = " << setw(20) << setprecision(15) << abserr
            << endl;


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
      cout << setw(11) << setprecision(5) << result-eps_mat[row][col];
    }
    cout << endl;
  }
*/
}

void
subinterval_double::dqk(
                          integr_fn_double f,
                          void* ex
                          )
{
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
    //***end prologue  dqk_double
    //
    //           the abscissae and weights are given for the interval (-1,1).
    //
    //           xgk    - abscissae of the 2*N+1-point kronrod rule
    //                    xgk(2), xgk(4), ...  abscissae of the N-point
    //                    gauss rule
    //                    xgk(1), xgk(3), ...  abscissae which are optimally
    //                    added to the N-point gauss rule
    //
    //           wgk    - weights of the 2N+1-point kronrod rule
    //
    //           wg     - weights of the N-point gauss rule
    //
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
    //***first executable statement  dqk_double
    //"
    int N = subinterval_double::ctl->N;
    double centr =  (a + b) / 2;
    double hlgth = (b - a) / 2;
    double dhlgth = fabs(hlgth);
    //
    //           compute the 2N+1-point kronrod approximation to
    //           the integral, and estimate the absolute error.
    //
    double resg = 0.;
    double resk = 0.;
    rabs = 0.;
    double fval = 0.;
    vector<double> fv(2*N+1);

    for (int j=0; j<2*N+1; j++) {
        fv[j]=centr+hlgth * subinterval_double::ctl->x_kronrod[j];
    }
    f(&fv[0], 2*N+1, ex);
    for (int j=0; j<2*N+1; j++) {
        fval = fv[j];
        resk += subinterval_double::ctl->w_kronrod[j] * fval;
        rabs += subinterval_double::ctl->w_kronrod[j] * fabs(fval);
    }
    for (int j=0; j<N; j++) {
        fval = fv[2*j+1];
        resg += subinterval_double::ctl->w_gauss[j]*fval;
    }
    double reskh = resk / 2;
    defabs = 0;
    for (int j=0; j<2*N+1 ; j++) {
        defabs += subinterval_double::ctl->w_kronrod[j] * fabs(fv[j] - reskh);
    }
    r = resk * hlgth;
    rabs = rabs * dhlgth;
    defabs = defabs * dhlgth;
    e = fabs((resk - resg) * hlgth);
    if (defabs != double(0) && e != double(0)) {
        e = defabs * std::min(1., std::pow((200 * e / defabs), 1.5));
    }
    if (rabs > uflow / (50 * epmach)) {
        e = std::max((epmach * 50) * rabs, e);
    }
}

void
dqagpe(
  integr_fn_double f,
  void* ex,
  const vector<double>& points,
  double& result,
  double& abserr,
  int& neval,
  int& ier,
  vector<subinterval_double>& subs,
  int& last)
{
  //***begin prologue  dqagpe_double
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
  //            points - double precision
  //                     vector of dimension at least 2 containing the
  //                     endpoints of the initial subintervals, which are
  //                     assumed to be in ascending order.
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
  //                     in the partition of (a,b), limit.ge.npts
  //                     if limit.lt.npts, the routine will end with
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
  //                             npts.lt.2 or
  //                             break points are specified outside
  //                             the integration range or
  //                             (epsabs.le.0 and
  //                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28))
  //                             or limit.lt.npts.
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
  //            last   - integer
  //                     number of subintervals actually produced in the
  //                     subdivisions process
  //
  //***references  (none)
  //***routines called  dqelg_double,dqk21_double,dqpsrt
  //***end prologue  dqagpe_double
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
  //***first executable statement  dqagpe_double
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
  int npts = int(points.size());
  int nint = npts-1;
  if (npts < 2 || subinterval_double::ctl->limit < nint) {
    cerr << "dqagp: Failed first test on input." << endl;
    if (npts < 2) cerr << "npts < 2" << endl;
    if (subinterval_double::ctl->limit < nint) cerr << "limit < nint" << endl;
    ier = 6;    // Invalid Input
    result = NAN;
    return;
  }
  if (npts==2 && points.at(0)==points.at(1)) {
        result =0;
        return;
  }
  for (vector<double>::const_iterator ppoint=points.begin(); ppoint < points.end()-1; ppoint++)
    if (*ppoint >= *(ppoint+1)) {
        cerr << "dqagp: points are either not distinct or not in ascending order" << endl;
        ier = 6;  //Invalid Input
        result = NAN;
        return;
    }

  //
  //            compute first integral and error approximations.
  //            ------------------------------------------------
  //
  double resabs = 0.0e+00;
  int i;
  for (i=1; i<=nint; i++) {
    subs.at(i-1).a = points.at(i-1);
    subs.at(i-1).b = points.at(i);
    subs.at(i-1).dqk(f, ex);
    abserr += subs.at(i-1).e;
    result += subs.at(i-1).r;
    if (subs.at(i-1).e == subs.at(i-1).defabs && subs.at(i-1).e != 0.0e+00) {
      subs.at(i-1).ndin = 1;
    } else {
      subs.at(i-1).ndin = 0;
    }
    resabs += subs.at(i-1).rabs;
    subs.at(i-1).level =0;
    min_length = std::min(min_length, fabs(subs.at(i-1).b-subs.at(i-1).a));
  }
  double errsum = 0.0e+00;
  for (i=1; i<=nint; i++) {
    errsum += subs.at(i-1).e;
  }
  //
  //           test on accuracy.
  //
  last = nint;
  neval = (2*subinterval_double::ctl->N+1) * nint;
  double dres = fabs(result);
  double errbnd = std::max(subinterval_double::ctl->epsabs, subinterval_double::ctl->epsrel * dres);
  if (abserr <= 100 * epmach * resabs && abserr > errbnd) {
    ier = 2;       // Roundoff error detected
  }
  make_heap(subs.begin(),subs.begin()+nint);
  if (subinterval_double::ctl->limit < npts) {
    ier = 1;         // Maximum number of subdivisions reached
  }
  if (ier != 0 || abserr <= errbnd) {
    if (ier > 2) {
      ier = ier - 1;
    }
    return;
  }
  //
  //           initialization
  //           --------------
  //
  vector<double> rlist2(50,0);
  rlist2.at(0) = result;
  vector<double> res3la(3,0);
  double correc = 0.;
  double area = result;
  int nres = 0;
  int numrl2 = 1;
  int ktmin = 0;
  bool extrap = false;
  bool noext = subinterval_double::ctl->noext;
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
  if (dres >= (1 - 50 * epmach) * resabs) {
    ksgn = 1;
  }
  bool final_check, test_for_divergence, need_sum;
  //
  //           main do-loop
  //           ------------
  //
  for (last=npts; last<=subinterval_double::ctl->limit ; last++) {
    //           At this point last is the number of intervals after the biscection.
    //
    //           bisect the subinterval with the largest score
    //           estimate.
    //
    pop_heap(subs.begin(), subs.begin()+last-1);
    int maxerr = last - 1;
    if (!subs.at(maxerr-1).is_divisible()) {
        last--;
        break;
    }
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
    subs.at(maxerr-1).dqk(f, ex);
    subs.at(last-1).dqk(f, ex);
    //
    //           improve previous approximations to integral
    //           and error and test for accuracy.

    neval += 2*(2*subinterval_double::ctl->N+1);
    double area12 = subs.at(maxerr-1).r + subs.at(last-1).r;
    double erro12 = subs.at(maxerr-1).e + subs.at(last-1).e;
    errsum += erro12 - errmax;
    area += area12 - rold;
    if (subs.at(maxerr-1).defabs != subs.at(maxerr-1).e
                && subs.at(last-1).defabs != subs.at(last-1).e
                && erro12 > subinterval_double::ctl->epsrel*fabs(area12)) {
      if ( fabs(rold - area12) <= 0.1e-04* fabs(
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
    push_heap(subs.begin(),subs.begin()+last-1);
    push_heap(subs.begin(),subs.begin()+last);
    errbnd = std::max(subinterval_double::ctl->epsabs, subinterval_double::ctl->epsrel * fabs(area));
    min_length = std::min(min_length,fabs(b2-a1)/2);
    //
    //           test for roundoff error and eventually set error flag.
    //
/*
    cout << "a1 = " << a1 << ", b2 = " << b2
              << ", level = " << level_cur <<endl
              << "rold = " << rold << ", rnew = " << area12 << endl
              << "errorold = " << errmax << ", errornew = " << erro12 << endl
              << "iroff1 = " << iroff1
              << ", iroff2 = " << iroff2
              << ", iroff3 = " << iroff3 << endl;
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
    if (last == subinterval_double::ctl->limit) {
      ier = 1;   // Maximum number of subdivisions reached
    }
    //
    //           set error flag in the case of bad integrand behaviour
    //           at a point of the integration range
    //
    if (fmax(fabs(a1), fabs(b2)) <= (1 +
        100 * epmach) * (fabs(a2) + 1000 * uflow)) {
      ier = 4;   // Bad integrand behavior
    }
    final_check=true;
    if (errsum <= errbnd) {
      // ***jump out of do-loop
 //     cout << "dqagpe: Success errsum <= errbnd" << endl;
      need_sum = true;
      final_check = false;
      break;
    }
    if (ier != 0) {
      // ***jump out of do-loop
//      cout << "dqagpe: Aborting with raw error code = " << ier << endl;
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
        dqelg(numrl2, rlist2, reseps, abseps, res3la, nres);
        ktmin++;
        if (ktmin > 10 && abserr < 0.1e-02 * errsum) {
          ier = 5;    // Roundoff error in extrapolation table
        }
        if (abseps < abserr) {
//          cout << "Using reseps = " << reseps
//                    << ", vs result = " << result << endl;
          ktmin = 0;
          abserr = abseps;
          result = reseps;
          correc = erlarg;
          ertest = std::max(subinterval_double::ctl->epsabs, subinterval_double::ctl->epsrel * fabs(reseps));
          // ***jump out of do-loop
          if (abserr < ertest) {
 //           cout << "dqagpe: Success.  abserr from dqelg < errtest." << endl;
            break;
          }
        }
        //
        //           prepare bisection of the smallest interval.
        //
        if (numrl2 == 1) {
//          cout << "Performed extrapolation but numrl2 == 1" << endl;
          noext = true;
        }
        if (ier >= 5) {
//          cout << "dqagpe: Aborting because of roundoff error in extrapolation table" << endl;
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
//        cout << "dqagpe: ierro = 3 adding " << correc << " to abserr." << endl;
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
//    cout << "Calculating sum" << endl;
    result = 0.0e+00;
    for (int k=1; k<=last; k++) {
      result += subs.at(k-1).r;
    }
    abserr = errsum;
  }
  if (ier > 2) {
    ier = ier - 1;
  }
}




