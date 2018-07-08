#include "mex.h"
#include "Polynomial/Polynomial.hpp"

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
    if ( nlhs != 1 || nrhs != 1 )
    {
        mexErrMsgTxt("usage: roots = roots_mex(coeffs)");
    }
    int m = mxGetM(prhs[0]);
    int n = mxGetN(prhs[0]);
    if ( m != 1 )
    {
        mexErrMsgTxt("input must be of size [1 n+1] for polynomial of degree n");
    }
    double lb = mxGetScalar(prhs[1]);
    double ub = mxGetScalar(prhs[2]);
    double *coeffs = (double*)mxGetData(prhs[0]);
    Eigen::VectorXd coef(n);
    for ( int i = 0; i < n; i++ ) coef[i] = coeffs[i];
    Polynomial::Polynomial<Eigen::Dynamic> poly(coef);
    std::vector<double> real_roots;
    poly.realRoots(real_roots);
    
    plhs[0] = mxCreateDoubleMatrix(1,real_roots.size(),mxREAL);
    double *roots = (double*)mxGetData(plhs[0]);
    for ( int i = 0; i < real_roots.size(); i++ ) roots[i] = real_roots[i];
}
