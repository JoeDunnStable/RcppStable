#include <RcppArmadillo.h>
#include <iostream>
using namespace Rcpp;


// [[Rcpp::depends(RcppArmadillo]]

class cubicspline {
private:
  arma::vec knots;
  arma::mat coefs;
public:
  cubicspline(arma::vec x, arma::vec y, bool endp2nd, arma::vec der)
  {
    unsigned int n = x.n_elem;
    if (n<2) {
      Rcpp::stop("cubic_spline: Number of knots must be at least 2.\n");
    }
    arma::vec h = arma::diff(x);
    if (any(h<=0))
      Rcpp::stop("cubic_spline; the knots must be distinct and in ascending order.\n");
    knots=x;
    if (n==2) {
      coefs.set_size(1,4);
      coefs.zeros();
      coefs(0,0)=y(0);
      coefs(0,1)=(x(1)!=x(0)) ? (y(1)-y(0))/(x(1)-x(0)) : 0;
      return;
    }
    arma::vec e(n);
    arma::mat A(n,n);
    A.zeros();
    arma::vec d(n-1);
    e(0)= 2 *h(0);
    e.subvec(1,n-2) = 2*(h.subvec(0,n-3) + h.subvec(1,n-2));
    e(n-1)= 2*h(n - 2);
    A.diag()=e;
    A.diag(-1)=h;
    A.diag(1)=h;
    d = arma::diff(y)/h;

    arma::vec rhs(n);
    rhs.subvec(1,n-2) = 3 * (d.subvec(1,n-2) - d.subvec(0,n - 3));
    double der0 = der(0);
    double dern = der(1);
    if (endp2nd) {
      A(0,0) = 2 * h[0];
      A(0,1) =  h(0);
      A(n-1, n-1) = 2 * h(n - 2);
      A(n-2,n-3) = h(n-2);
      rhs(0) = 3 * (d(0) - der0);
      rhs(n-1) = 3 * (dern - d(n-2));
    }
    else {
      A.row(0).zeros();
      A(0,0) = 1;
      A.row(n-1).zeros();
      A(n-1, n-1) = 1;
      rhs(0) = der0;
      rhs(n-1)= dern;
    }

    coefs.set_size(n, 4);
    coefs.zeros();
    coefs.col(2) = solve(A, rhs);
    unsigned int m;
    for (m=0; m<n-1; m++) {
      coefs(m, 3) = (coefs(m+1,2) - coefs(m, 2))/3/h(m);
      coefs(m,1) = d(m) - h(m)/3 * (coefs(m + 1, 2) + 2 * coefs(m, 2));
      coefs(m,0) = y(m);
    }
    coefs.resize(n-1,4);
  } //constructor

  arma::vec operator() (arma::vec x){
    arma::vec ret(x.n_elem);
    arma::uvec s_index = sort_index(x);
    arma::uvec indices(x.n_elem);
    indices.zeros();
    arma::ivec isgood(x.n_elem);
    isgood.ones();
    arma::uword i=0, j=0;
    for (j=0; j<x.n_elem && (x(s_index(j))< knots(0)) ; j++)
      isgood(s_index(j)) = 0;
    while (j<x.n_elem){
      if (x(s_index(j))==knots(i)){
        indices(s_index(j)) = std::min(i,knots.n_elem-2);
        j++;
      } else if (x(s_index(j)) > knots(i)){
        i++;
        if (i>knots.n_elem-1){
          for (;j<x.n_elem;j++)
            isgood(s_index(j))=0;
          break;
        }
      } else {
        indices(s_index(j)) = std::min(i-1,knots.n_elem-2);
        j++;
      }
    }
    for (i=0; i<x.n_elem; i++) {
      if (isgood(i)) {
        j=indices(i);
        double dx=x(i)-knots(j);
        ret(i) = ((coefs(j,3)*dx + coefs(j,2))*dx + coefs(j,1))*dx + coefs(j,0);
      } else
        ret(i) = R_NaReal;
    }
    return ret;
  }
  void get_knots() {
    Rcpp::Environment global = Rcpp::Environment::global_env();
    global["knots_cpp"]=wrap(knots);
  }
  void get_coefs() {
    Rcpp::Environment global = Rcpp::Environment::global_env();
    global["coefs_cpp"]=wrap(coefs);
  }
};

