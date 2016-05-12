//
//  kronrod.h
//  
//
//  Created by Joseph Dunn on 12/29/15.
//
//

#ifndef kronrod_h
#define kronrod_h
#include <vector>
using std::vector;

void toms726(const int N, const vector<double>& a, const vector<double>& b, vector<double>& x, vector<double>& w,
             const int verbose);
void r_jacobi01(const int N, const double a, const double b, vector<double>& c, vector<double>& d);
void r_jacobi(const int N, const double a, const double b, vector<double>& c, vector<double>& d);
void r_kronrod(int N, const vector<double>& a0, const vector<double>& b0, vector<double>& a, vector<double>& b);


#endif /* kronrod_h */
