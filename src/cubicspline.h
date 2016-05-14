#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include <algorithm>
using std::sort;

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Mat;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vec;
typedef Eigen::Matrix<unsigned int, Eigen::Dynamic, 1> uVec;
typedef Eigen::Matrix<int, Eigen::Dynamic, 1> iVec;
typedef Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, unsigned int> pMat;

class CompareIndicesByAnotherVectorValues {
    const Vec& _values;
    public: CompareIndicesByAnotherVectorValues(const Vec &values) : _values(values) {}
    public: bool operator() (const unsigned int& a, const unsigned int& b) const { return (_values)(a) < (_values)(b); }
};

inline uVec sort_indexes(const Vec &v) {

    CompareIndicesByAnotherVectorValues comp(v);

    // initialize original index locations
    uVec idx(v.size());
    for (unsigned int i = 0; i != idx.size(); ++i) idx[i] = i;

    // sort indexes based on comparing values in v
    sort(idx.data(), idx.data()+idx.size(), comp);

    return idx;
}

class cubicspline {
private:
  Vec knots;
  Mat coefs;
public:
    cubicspline() {};
  cubicspline(Vec x, Vec y, bool endp2nd, Vec der)
  {
    unsigned int n = static_cast<unsigned int>(x.size());
    if (n<2) {
        return;
    }
    Vec h = x.tail(n-1)-x.head(n-1);
    if ((h.array()<=0).any())
        throw std::range_error("cubic_spline; the knots must be distinct and in ascending order.\n");
    knots=x;
    if (n==2) {
      coefs.setZero(1,4);
      coefs(0,0)=y(0);
      coefs(0,1)=(x(1)!=x(0)) ? (y(1)-y(0))/(x(1)-x(0)) : 0;
      return;
    }
    Vec e(n);
    Mat A;
    A.setZero(n,n);
    Vec d(n-1);
    e(0)= 2 *h(0);
    e.segment(1,n-2) = 2*(h.segment(0,n-2) + h.segment(1,n-2));
    e(n-1)= 2*h(n - 2);
    A.diagonal()=e;
    A.diagonal(-1)=h;
    A.diagonal(1)=h;
    d = (y.tail(n-1)-y.head(n-1)).array()/h.array();

    Vec rhs(n);
    rhs.segment(1,n-2) = 3 * (d.segment(1,n-2) - d.segment(0,n - 2));
    double der0 = endp2nd ? der(0) : 0;
    double dern = endp2nd ? der(1) : 0;
    if (endp2nd) {
      A(0,0) = 2 * h[0];
      A(0,1) =  h(0);
      A(n-1, n-1) = 2 * h(n - 2);
      A(n-2,n-3) = h(n-2);
      rhs(0) = 3 * (d(0) - der0);
      rhs(n-1) = 3 * (dern - d(n-2));
    }
    else {
      A.topRows(1).setZero();
      A(0,0) = 1;
      A.bottomRows(1).setZero();
      A(n-1, n-1) = 1;
      rhs(0) = der0;
      rhs(n-1)= dern;
    }

    coefs.setZero(n, 4);
      coefs.col(2) = A.colPivHouseholderQr().solve(rhs);
    unsigned int m;
    for (m=0; m<n-1; m++) {
      coefs(m, 3) = (coefs(m+1,2) - coefs(m, 2))/3/h(m);
      coefs(m,1) = d(m) - h(m)/3 * (coefs(m + 1, 2) + 2 * coefs(m, 2));
      coefs(m,0) = y(m);
    }
    coefs.conservativeResize(n-1,4);
//      cout << "knots:" << endl << knots << endl;
//      cout << "y_knots" << endl << y << endl;
//      cout << "coefs: " << endl << coefs << endl;
  } //constructor

  Vec operator() (Vec x){
    Vec ret(x.size());
    uVec s_index = sort_indexes(x);
    uVec indices;
    indices.setZero(x.size());
    iVec isgood;
    isgood.setOnes(x.size());
    unsigned int i=0, j=0;
    for (j=0; j<x.size() && (x(s_index(j))< knots(0)) ; j++)
      isgood(s_index(j)) = 0;
    while (j<x.size()){
      if (x(s_index(j))==knots(i)){
        indices(s_index(j)) = std::min(i,static_cast<unsigned int>(knots.size())-2);
        j++;
      } else if (x(s_index(j)) > knots(i)){
        i++;
        if (i>knots.size()-1){
          for (;j<x.size();j++)
            isgood(s_index(j))=0;
          break;
        }
      } else {
        indices(s_index(j)) = std::min(i-1,static_cast<unsigned int>(knots.size())-2);
        j++;
      }
    }
    for (i=0; i<x.size(); i++) {
      if (isgood(i)) {
        j=indices(i);
        double dx=x(i)-knots(j);
        ret(i) = ((coefs(j,3)*dx + coefs(j,2))*dx + coefs(j,1))*dx + coefs(j,0);
      } else
        ret(i) = NAN;
    }
    return ret;
  }
  Vec get_knots() {
     return knots;
  }
  Mat get_coefs() {
      return coefs;
  }
  unsigned int get_n_knots() {
      return knots.size();
  }
};

