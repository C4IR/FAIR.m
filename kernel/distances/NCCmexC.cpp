/*
 * (c) Jan Modersitzki and Fabian Gigengack 2010/12/27, see FAIR.2 and FAIRcopyright.m.
 * http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
 * http://www.uni-muenster.de/EIMI/
 *
 * CPP code for Normalized Cross Correlation based distance measure.
 * See NCC for details.
 */

#include <math.h>
#include <mex.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    bool doDerivative = nlhs>2;
    if (nrhs<2)
        mexErrMsgTxt("Number of arguments must be 2");
    
    const double* T  = static_cast<double*>(mxGetData(prhs[0]));
    const double* R  = static_cast<double*>(mxGetData(prhs[1]));
    
    const int N = mxGetM(prhs[1]);
    mwSize dims[2]; dims[0] = N; dims[1] = 1;
    
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    double* Dc = static_cast<double*>(mxGetData(plhs[0]));
    plhs[1] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
    double* rc = static_cast<double*>(mxGetData(plhs[1]));
    
    double* dD = 0;
    double* dr = 0;
    double* d2psi = 0;
    if (doDerivative) {
        dims[0] = 1; dims[1] = N;
        plhs[2] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
        dD = static_cast<double*>(mxGetData(plhs[2]));
        plhs[3] = mxCreateDoubleMatrix(1, 1, mxREAL);
        dr = static_cast<double*>(mxGetData(plhs[3]));
        dr[0] = 1;
        plhs[4] = mxCreateDoubleMatrix(1, 1, mxREAL);
        d2psi = static_cast<double*>(mxGetData(plhs[4]));
    }
    
    double nR  = 0; // norm of R
    double nT2 = 0; // squared norm of T
    Dc[0] = 0;
    for (int i=0;i<N;i++) {
        nR    += R[i] * R[i];
        nT2   += T[i] * T[i];
        Dc[0] += R[i] * T[i];
        rc[i]  = T[i];
    }
    nR = sqrt(nR);
    
    if (doDerivative) {
        double val = Dc[0] / nR / nT2;
        double val1 = 2 * val / nR;
        double val2 = 2 * val * val;
        for (int i=0;i<N;i++) {
            dD[i] = (-val1 * R[i] + val2 * T[i]) * dr[0];
        }
        d2psi[0] = 2 / nT2;
    }
    
    Dc[0] = 1 - Dc[0] * Dc[0] / nT2 / nR / nR;
}

/* =========================================================================== */
