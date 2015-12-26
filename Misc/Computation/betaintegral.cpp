#include <iostream>
#include <ostream>
#include <vector>
#include <math.h>
#include <iomanip>
#include <algorithm>
#include <random>
#include <functional>

#ifdef MATLAB_MEX_FILE
#include "mex.h"
#include "matrix.h"
#endif // MATLAB_MEX_FILE

using namespace std;


// Evaluates a two-term hypergeometric series of the form
//   sum_n [(a)_n (b)_n] / [(c)_n (d)_n]
double hyp2(double a, double b, double c, double d, size_t maxi)
{
    // Tolerance - series truncated when terms are below tol
    const double tol = 1e-15;

    // Accumulative sum of series
    double acc = 1;

    // Current and previous term
    double bcc = 1, bcc_ = 0;

    // Iteration variable
    size_t i = 0;

    // Continue while term is greater than tol
    while (bcc > tol) {
        // Save last term
        bcc_ = bcc;

        // Update term
        bcc *= ((a+i)*(b+i)) / ((c+i)*(d+i));

        // Accumulate series sum
        acc += bcc;

        // Break if maximum number of terms exceeded
        if (++i >= maxi) {
            cout << "Series did not converge, bcc=" << bcc << endl;
            cout << "  a=" << a << ", b=" << b << ", c=" << c << ", d=" << d << endl;
            break;
        }
    }
    //cout << "Terms evaluated: " << i << endl;

    // Estimate residual by fitting a line to the last two log-terms
    double residual = bcc*bcc/(bcc_-bcc);

    // Return result
    return log(acc + residual);
}

// Compute the difference of of two log numbers
//   log(exp(a)-exp(b))
double lsub(double a, double b)
{
    double maxab = max(a,b);
    return maxab + log(exp(a-maxab)-exp(b-maxab));
}

// Computes the value of the incomplete beta-beta integral
double betaincint(double a, double b, double c, double d, size_t maxi) // No input can be zero
{
    //cout << "Evaluating beta integral" << endl;
    //cout << "  a=" << a << ", b=" << b << ", c=" << c << ", d=" << d << endl;

    // Quality score of different series expansions
    vector<double> seriesDecay {(a+c)*(c+d)*(a+1)*(b+1)*(d+1),
                                (b+d)*(b+a)*(a+1)*(c+1)*(d+1),
                                (a+c)*(a+b)*(b+1)*(c+1)*(d+1),
                                (b+d)*(d+c)*(a+1)*(b+1)*(c+1)
                               };
    // Index of best series expansion
    size_t i = min_element(seriesDecay.begin(), seriesDecay.end()) - seriesDecay.begin();
    //cout << "Evaluating series " << i << endl;

    double H, logT, logI, logP;
    // Series 0
    if (i==0) {
        H = hyp2(c+a, c+d, c+1, a+b+c+d, maxi) - log(c);
        logP = lgamma(a+c) + lgamma(b+d) - lgamma(a+b+c+d) + H;
    }
    // Series 1
    else if (i==1) {
        H = hyp2(b+a, b+d, b+1, a+b+c+d, maxi) - log(b);
        logP = lgamma(a+c) + lgamma(b+d) - lgamma(a+b+c+d) + H;
    }
    // Series 2 or 3
    else {
        // Series 2
        if (i==2) {
            H = hyp2(a+b, a+c, a+1, a+b+c+d, maxi) - log(a);
        }
        // Series 3
        else {
            H = hyp2(d+b, d+c, d+1, a+b+c+d, maxi) - log(d);
        }
        logI = lgamma(a+c) + lgamma(b+d) - lgamma(a+b+c+d) + H;
        logT = lgamma(a) + lgamma(b) + lgamma(c) + lgamma(d) - lgamma(a+b) - lgamma(c+d);
        logP = lsub(logT, logI);
    }
    return logP;
}

#ifdef MATLAB_MEX_FILE
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    // Check inputs
    if(nrhs<4) {
        mexErrMsgIdAndTxt("betaintegral:nrhs","Four inputs required.");
    }
    mwSize M, N;
    M = mxGetM(prhs[0]);
    N = mxGetN(prhs[0]);
    for (size_t i = 0; i != 4; ++i) {
        if (mxGetM(prhs[i]) != M || mxGetN(prhs[i]) != N) {
            mexErrMsgIdAndTxt("betaintagral:nrhs", "Dimensions of first four inputs must be identical");
        }
    }

    // Inputs
    double *a, *b, *c, *d;
    size_t imax = 10000;

    // Get inputs
    a = mxGetPr(prhs[0]);
    b = mxGetPr(prhs[1]);
    c = mxGetPr(prhs[2]);
    d = mxGetPr(prhs[3]);
    if (nrhs>=5)
    {
        imax = mxGetScalar(prhs[4]);
    }

    // Create output
    plhs[0] = mxCreateDoubleMatrix(M, N, mxREAL);
    double *out;
    out = mxGetPr(plhs[0]);

    for (size_t i = 0; i != M*N; ++i) {
        out[i] = betaincint(a[i], b[i], c[i], d[i], imax);
    }
}
#endif // MATLAB_MEX_FILE


// Compute a large number of random integrals and prints out the sum
int main()
{
    // Maximum number of terms in series expansion
    size_t imax = 10000;

    // Random number generator
    default_random_engine generator;
    uniform_int_distribution<int> distribution(1, 10000);
    auto u = bind(distribution, generator);

    // Accumulated sum
    double acc = 0;

    // Loop over random integrals
    for (size_t k = 0; k != 500000; ++k)
    {
        acc += betaincint(u(), u(), u(), u(), imax);
    }

    // Print results
    cout << acc << endl;

    return 0;
}
