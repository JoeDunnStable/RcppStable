
#include "kronrod_double.h"

#include <Accelerate/Accelerate.h>

#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <cmath>

using std::cout;
using std::endl;
using std::setw;
using std::setprecision;
using std::right;
using std::vector;

template<typename T>
class sort_data {
public:
    T a;
    int index;
    sort_data() : a(0), index(0){}
    sort_data(T a, int index) : a(a), index(index){}
};

template<typename T>
bool sd_lt(const sort_data<T> &lhs,const sort_data<T>& rhs) {
    return lhs.a < rhs.a;
}

// toms726 generates the nodes x and weights w for a quadrature rule
// given the recurrence factors a and b for the orthogonal polynomials
//
void toms726(const int N, const vector<double>& a, const vector<double>& b, vector<double>& x, vector<double>& w,
             const int verbose){
    int sumbgt0=0;
    for (int j=0; j<N; j++) sumbgt0+=(b.at(j)>0);
    if(sumbgt0!=N) {
        throw std::range_error("toms726: b is not all >0");
    }
    if (N == 1) {
        x.at(0)=a.at(0);
        w.at(0)=b.at(0);
        return;
    }
    vector<double> J(N*N, 0);
    for (int k=0; k<N-1; k++){
        J.at(k+k*N)=a.at(k);
        J.at(k+(k+1)*N)=sqrt(b.at(k+1));
        J.at(k+1+k*N)=J.at(k+(k+1)*N);
    }
    J.at(N-1+(N-1)*N)=a.at(N-1);
    int lwork, info;
    lwork = -1;
    vector<double> work(1);
    vector<double> d(N), e(N);
    int N_tmp = N;
    char jobz[] = "V";
    char uplo[] = "U";
    dsyev_(jobz, uplo, &N_tmp, &J[0], &N_tmp, &d[0], &work[0], &lwork, &info);
    lwork = (int) double (work[0]);
    work.resize(std::max(1, lwork));
    //inverse matrix
    dsyev_(jobz, uplo, &N_tmp, &J[0], &N_tmp, &d[0], &work[0], &lwork, &info);
    for (int k=0; k<N; k++) {
        e.at(k)=b.at(0)*(pow(J.at(k*(N)),2));
    }
    vector<sort_data<double> > srt_d(N);
    for (int i=0; i<N; i++) {
        sort_data<double> elem(d.at(i),i);
        srt_d.at(i) = elem;
    }
    std::sort(srt_d.begin(), srt_d.end(),sd_lt<double>);
    for (int i=0; i<N; i++){
        x.at(i)=srt_d.at(i).a;
        w.at(i)=e.at(srt_d.at(i).index);
    }
}

// R_JACOBI01 Recurrence coefficients for monic Jacobi polynomials
// on [0,1].
//
//    R_JACOBI01(n,a,b,c,d) generates the first n recurrence
//    coefficients for monic Jacobi polynomials on [0,1] with
//    parameters a and b. These are orthogonal on [0,1] relative
//    to the weight function w(t)=(1-t)^a t^b. The n alpha-
//    coefficients are stored in c, the n beta-
//    coefficients in d.
//
//    Translated to C++ by Joseph Dunn, Dec. 26, 2015
//
void r_jacobi01(const int N, const double a, const double b, vector<double>& c, vector<double>& d){
    if((N<=0)||(a<=-1)||(b<=-1)){
        throw std::range_error("r_jacobi01: parameter(s) out of range");
    }
    r_jacobi(N, a, b, c, d);
    for (int i=0; i<N; i++) {
        c.at(i)=(1+c.at(i))/2;
        if (i==0)
            d.at(i)=d.at(0)/pow(2,a+b+1);
        else
            d.at(i)=d.at(i)/4;
    }
}

// R_JACOBI Recurrence coefficients for monic Jacobi polynomials.
//
//    R_JACOBI(n,a,b,c,d) generates the first n recurrence
//    coefficients for monic Jacobi polynomials with parameters
//    a and b. These are orthogonal on [-1,1] relative to the
//    weight function w(t)=(1-t)^a(1+t)^b. The n alpha-coefficients
//    are stored in c, the n beta-coefficients in d
//
//    Translated to C++ by Joseph Dunn, Dec. 26, 2015
//    Supplied by Dirk Laurie, 6-22-1998; edited by Walter
//    Gautschi, 4-4-2002.
//
void r_jacobi(const int N, const double a, const double b, vector<double>& c, vector<double>& d){
    if((N<=0)||(a<=-1)||(b<=-1)) {
        throw std::range_error("r_jacobi: parameter(s) out of range");
    }
    double nu = (b-a)/(a+b+2);
    double mu = pow(2,a+b+1)*tgamma(a+1)*tgamma(b+1)/tgamma(a+b+2);
    if (N==1) {
      c.at(0)=nu;
      d.at(0)=mu;
      return;
    }
    c.at(0)=nu;
    for (int n=1; n<N; n++){
        double nab=2*n+a+b;
        c.at(n) = (pow(b,2)-pow(a,2))*1/(nab*(nab+2));
    }
    d.at(0) = mu;
    d.at(1) = 4*(a+1)*(b+1)/(pow(a+b+2,2)*(a+b+3));
    for (int n=2; n<N; n++) {
        double nab=2*n+a+b;
        d.at(n)=4*(n+a)*(n+b)*n*(n+a+b)/(pow(nab,2)*(nab+1)*(nab-1));
    }
}

// R_KRONROD Jacobi-Kronrod matrix.
//
//    ab=R_KRONROD(n,a0,b0,a,b) produces the alpha- and beta-elements in
//    the Jacobi-Kronrod matrix of order 2n+1 for the weight
//    function (or measure) w. The input data for the weight
//    function w are the recurrence coefficients of the associated
//    orthogonal polynomials, which are stored in the vectors a0 & b0 of
//    dimension [ceiling(3*n/2)+1]. The alpha-elements are stored in
//    a, the beta-elements in b, of
//    dimension (2*n+1)
//
//    Supplied by Dirk Laurie, 6-22.1998
//
//    Translated to C++ by Joseph Dunn, Dec. 26 2015
//    Edited by Pavel Holoborodko, November 7, 2011:
//	 Added arbitrary precision support using
//	 Multiprecision Computing Toolbox for MATLAB
//	 http://advanpix.com
//
void r_kronrod(const int N, const vector<double>& a0, const vector<double>& b0, vector<double>& a, vector<double>& b) {
    if (a0.size()<ceil(3.*N/2.)+1) {
        throw std::range_error("r_kronrod: a0 is too short.");
    }
    for (int k=0; k<=ceil(3.*N/2.); k++) {
        a.at(k) = a0.at(k);
        b.at(k) = b0.at(k);
    }
    vector<double> s((N/2)+2);
    vector<double> t((N/2)+2);
    for (int k=0; k<(N/2)+2; k++)
        t.at(k)=s.at(k)=double(0);
    t.at(1)=b[N+1];	// P.H.
    for (int m=0; m<=(N-2); m++) {
        double cumsum = 0;
        vector<double> tmp((N/2)+2);
        tmp.at(0)=s.at(0);
        for (int k=((m+1)/2); k>=0; k--){
            int l = m-k;
            cumsum=cumsum+(a.at(k+N+1)-a.at(l))*t.at(k+1)+b.at(k+N+1)*s.at(k)-b.at(l)*s.at(k+1);
            tmp.at(k+1)=cumsum;
        }
        for (int k=((m+1)/2)+2; k<((N/2)+2); k++)
            tmp.at(k)=s.at(k);
        for (int k=0; k<((N/2)+2); k++){
            s.at(k)=t.at(k);
            t.at(k)=tmp.at(k);
	}
    }
    for (int j=(N/2); j>=0; j--)
        s.at(j+1)=s.at(j);
   for (int m=N-1; m<=(2*N-3); m++){
        double cumsum = 0;
        vector<double> tmp((N/2)+2);
        for (int k=m+1-N; k<=((m-1)/2); k++){
            int l=m-k;
            int j=N-1-l;
            cumsum=cumsum+-(a.at(k+N+1)-a.at(l))*t.at(j+1)-b.at(k+N+1)*s.at(j+1)+b.at(l)*s.at(j+2);
            tmp.at(j+1)=cumsum;
        }
        int j;
        for (int k=m+1-N; k<=((m-1)/2); k++){
            j=N-1-m+k;
            s.at(j+1)=tmp.at(j+1);
        }
        int k=((m+1)/2);
        if ((m % 2)==0) {
            a.at(k+N+1)=a.at(k)+(s.at(j+1)-b.at(k+N+1)*s.at(j+2))/t.at(j+2);
        }
        else {
            b.at(k+N+1)=s.at(j+1)/s.at(j+2);
        }
        for (int j=0; j<((N/2)+2); j++) {
            double swap = s.at(j);
            s.at(j)=t.at(j);
            t.at(j)=swap;
        }
    }
    a.at(2*N)=a.at(N-1)-b.at(2*N)*s.at(1)/t.at(1);
    return;
}
