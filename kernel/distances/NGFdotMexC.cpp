/*
 * (c) Jan Modersitzki and Fabian Gigengack 2010/12/27, see FAIR.2 and FAIRcopyright.m.
 * http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
 * http://www.uni-muenster.de/EIMI/
 *
 * CPP code for Normalized Gradient Field based distance measure.
 * See NGFdot for details.
 */

#include <math.h>
#include <mex.h>

void NGF2D(double* Dc, double* rc, double* dD, double* drcP,
        mwIndex* drcI, mwIndex* drcJ, double* d2psi, const double* T,
        const double* R, const double* o, const double* m, const double* h,
        const double* edge, double* fac, bool doDerivative);
void NGF3D(double* Dc, double* rc, double* dD, double* drcP,
        mwIndex* drcI, mwIndex* drcJ, double* d2psi, const double* T,
        const double* R, const double* o, const double* m, const double* h,
        const double* edge, double* fac, bool doDerivative);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    mwIndex *drcI, *drcJ, k;
    double *drcP, *Dc, *rc, *dD, *d2psi, *fac;
    int dim;
    bool doDerivative;
    
    doDerivative = nlhs>2;
    if (nrhs<6)
        mexErrMsgTxt("Number of arguments must be 6");
    
    const double* T = static_cast<double*>(mxGetData(prhs[0]));
    const double* R = static_cast<double*>(mxGetData(prhs[1]));
    const double* o = static_cast<double*>(mxGetData(prhs[2]));
    const double* m = static_cast<double*>(mxGetData(prhs[3]));
    const double* h = static_cast<double*>(mxGetData(prhs[4]));
    const double* edge = static_cast<double*>(mxGetData(prhs[5]));
    
    const int N = mxGetM(prhs[1]);
    mwSize dims[2]; dims[0] = N; dims[1] = 1;
    
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    Dc = static_cast<double*>(mxGetData(plhs[0]));
    plhs[1] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
    rc = static_cast<double*>(mxGetData(plhs[1]));
    
    dD = 0;
    d2psi = 0;
    if (doDerivative) {
        dims[0] = 1; dims[1] = N;
        plhs[2] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
        dD = static_cast<double*>(mxGetData(plhs[2]));
        plhs[3] = mxCreateSparse(N, N, 6*N, mxREAL);
        drcP = mxGetPr(plhs[3]);
        drcI = mxGetIr(plhs[3]);
        drcJ = mxGetJc(plhs[3]);
        plhs[4] = mxCreateDoubleMatrix(1, 1, mxREAL);
        d2psi = static_cast<double*>(mxGetData(plhs[4]));
    }
    
    dims[0] = 3*N; dims[1] = 1;
    fac = 0;
    if (doDerivative) {
        fac = static_cast<double*>(mxGetData(mxCreateNumericArray(2, dims,
                mxDOUBLE_CLASS, mxREAL)));
    }
    
    dim = mxGetN(prhs[4]);
    switch (dim) {
        case 2:
            NGF2D(Dc, rc, dD, drcP, drcI, drcJ, d2psi, T, R, o, m, h, edge,
                    fac, doDerivative);
            break;
        case 3:
            NGF3D(Dc, rc, dD, drcP, drcI, drcJ, d2psi, T, R, o, m, h, edge,
                    fac, doDerivative);
            break;
        default:
            break;
    }
}

void NGF2D(double* Dc, double* rc, double* dD, double* drcP,
        mwIndex* drcI, mwIndex* drcJ, double* d2psi, const double* T,
        const double* R, const double* o, const double* m, const double* h,
        const double* edge, double* fac, bool doDerivative) {
    double hd, gradT1, gradT2, gradR1, gradR2;
    double lengthGT, lengthGR, r1, r2;
    int idx, m0, m1;
    mwIndex k;
    
    m0 = (int) m[0];
    m1 = (int) m[1];
    const int N = m0 * m1;
    hd = h[0] * h[1];
    Dc[0] = 0;
    
    for (int y=0;y<m1;y++) {
        for (int x=0;x<m0;x++) {
            idx = x+y*m0;
            if (x==0) {
                gradT1 = (T[idx+1] - T[idx]) / 2 / h[0];
                gradR1 = (R[idx+1] - R[idx]) / 2 / h[0];
            } else if (x==m0-1) {
                gradT1 = (T[idx] - T[idx-1]) / 2 / h[0];
                gradR1 = (R[idx] - R[idx-1]) / 2 / h[0];
            } else {
                gradT1 = (T[idx+1] - T[idx-1]) / 2 / h[0];
                gradR1 = (R[idx+1] - R[idx-1]) / 2 / h[0];
            }
            if (y==0) {
                gradT2 = (T[idx+m0] - T[idx]) / 2 / h[1];
                gradR2 = (R[idx+m0] - R[idx]) / 2 / h[1];
            } else if (y==m1-1) {
                gradT2 = (T[idx] - T[idx-m0]) / 2 / h[1];
                gradR2 = (R[idx] - R[idx-m0]) / 2 / h[1];
            } else {
                gradT2 = (T[idx+m0] - T[idx-m0]) / 2 / h[1];
                gradR2 = (R[idx+m0] - R[idx-m0]) / 2 / h[1];
            }
            lengthGT = sqrt(gradT1 * gradT1 + gradT2 * gradT2
                          + edge[0] * edge[0]);
            lengthGR = sqrt(gradR1 * gradR1 + gradR2 * gradR2
                          + edge[0] * edge[0]);
            r1 = gradR1 * gradT1 + gradR2 * gradT2;
            r2 = 1 / (lengthGT * lengthGR);
            if (doDerivative) {
                fac[idx]   = (r2 * gradR1 + r1 * (-1)
                        / (lengthGR * lengthGT * lengthGT *  lengthGT)
                        * gradT1) / 2 / h[0];
                fac[idx+N] = (r2 * gradR2 + r1 * (-1)
                        / (lengthGR * lengthGT * lengthGT *  lengthGT)
                        * gradT2) / 2 / h[1];
            }
            rc[idx] = r1 * r2;
            Dc[0]  += rc[idx] * rc[idx];
        }
    }
    Dc[0] = o[0] - hd * Dc[0];
    
    
    if (!doDerivative) {
        return;
    }
    k = 0;
    for (int y=0;y<m1;y++) {
        for (int x=0;x<m0;x++) {
            idx = x+y*m0;
            drcJ[idx] = k;
            dD[idx] = 0;
            if (y==0) {
                if (x==0) {
                    drcI[k] = idx;
                    drcP[k] = - fac[idx] - fac[idx+N];
                    dD[idx] += rc[drcI[k]] * drcP[k];
                    k++;
                    drcI[k] = idx+1;
                    drcP[k] = - fac[idx+1];
                    dD[idx] += rc[drcI[k]] * drcP[k];
                    k++;
                } else if (x==m0-1) {
                    drcI[k] = idx-1;
                    drcP[k] = fac[idx-1];
                    dD[idx] += rc[drcI[k]] * drcP[k];
                    k++;
                    drcI[k] = idx;
                    drcP[k] = fac[idx] - fac[idx+N];
                    dD[idx] += rc[drcI[k]] * drcP[k];
                    k++;
                } else {
                    drcI[k] = idx-1;
                    drcP[k] = fac[idx-1];
                    dD[idx] += rc[drcI[k]] * drcP[k];
                    k++;
                    drcI[k] = idx;
                    drcP[k] = - fac[idx+N];
                    dD[idx] += rc[drcI[k]] * drcP[k];
                    k++;
                    drcI[k] = idx+1;
                    drcP[k] = - fac[idx+1];
                    dD[idx] += rc[drcI[k]] * drcP[k];
                    k++;
                }
                drcI[k] = idx+m0;
                drcP[k] = - fac[idx+N+m0];
                dD[idx] += rc[drcI[k]] * drcP[k];
                k++;
            } else if (y==m1-1) {
                drcI[k] = idx-m0;
                drcP[k] = fac[idx+N-m0];
                dD[idx] += rc[drcI[k]] * drcP[k];
                k++;
                if (x==0) {
                    drcI[k] = idx;
                    drcP[k] = - fac[idx] + fac[idx+N];
                    dD[idx] += rc[drcI[k]] * drcP[k];
                    k++;
                    drcI[k] = idx+1;
                    drcP[k] = - fac[idx+1];
                    dD[idx] += rc[drcI[k]] * drcP[k];
                    k++;
                } else if (x==m0-1) {
                    drcI[k] = idx-1;
                    drcP[k] = fac[idx-1];
                    dD[idx] += rc[drcI[k]] * drcP[k];
                    k++;
                    drcI[k] = idx;
                    drcP[k] = fac[idx] + fac[idx+N];
                    dD[idx] += rc[drcI[k]] * drcP[k];
                    k++;
                } else {
                    drcI[k] = idx-1;
                    drcP[k] = fac[idx-1];
                    dD[idx] += rc[drcI[k]] * drcP[k];
                    k++;
                    drcI[k] = idx;
                    drcP[k] = fac[idx+N];
                    dD[idx] += rc[drcI[k]] * drcP[k];
                    k++;
                    drcI[k] = idx+1;
                    drcP[k] = - fac[idx+1];
                    dD[idx] += rc[drcI[k]] * drcP[k];
                    k++;
                }
            } else {
                drcI[k] = idx-m0;
                drcP[k] = fac[idx+N-m0];
                dD[idx] += rc[drcI[k]] * drcP[k];
                k++;
                if (x==0) {
                    drcI[k] = idx;
                    drcP[k] = - fac[idx];
                    dD[idx] += rc[drcI[k]] * drcP[k];
                    k++;
                    drcI[k] = idx+1;
                    drcP[k] = - fac[idx+1];
                    dD[idx] += rc[drcI[k]] * drcP[k];
                    k++;
                } else if (x==m0-1) {
                    drcI[k] = idx-1;
                    drcP[k] = fac[idx-1];
                    dD[idx] += rc[drcI[k]] * drcP[k];
                    k++;
                    drcI[k] = idx;
                    drcP[k] = fac[idx];
                    dD[idx] += rc[drcI[k]] * drcP[k];
                    k++;
                } else {
                    drcI[k] = idx-1;
                    drcP[k] = fac[idx-1];
                    dD[idx] += rc[drcI[k]] * drcP[k];
                    k++;
                    drcI[k] = idx+1;
                    drcP[k] = - fac[idx+1];
                    dD[idx] += rc[drcI[k]] * drcP[k];
                    k++;
                }
                drcI[k] = idx+m0;
                drcP[k] = - fac[idx+N+m0];
                dD[idx] += rc[drcI[k]] * drcP[k];
                k++;
            }
            dD[idx] = -2 * hd * dD[idx];
        }
    }
    drcJ[m0*m1] = k;
    d2psi[0]    = 2 * hd;
}

void NGF3D(double* Dc, double* rc, double* dD, double* drcP,
        mwIndex* drcI, mwIndex* drcJ, double* d2psi, const double* T,
        const double* R, const double* o, const double* m, const double* h,
        const double* edge, double* fac, bool doDerivative) {
    double hd, gradT1, gradT2, gradT3, gradR1, gradR2, gradR3;
    double lengthGT, lengthGR, r1, r2;
    int idx, m0, m1, m2;
    mwIndex k;
    
    m0 = (int) m[0];
    m1 = (int) m[1];
    m2 = (int) m[2];
    const int N = m0 * m1 * m2;
    hd = h[0] * h[1] * h[2];
    Dc[0] = 0;
    
    for (int z=0;z<m2;z++) {
        for (int y=0;y<m1;y++) {
            for (int x=0;x<m0;x++) {
                idx = x+y*m0+z*m0*m1;
                if (x==0) {
                    gradT1 = (T[idx+1] - T[idx]) / 2 / h[0];
                    gradR1 = (R[idx+1] - R[idx]) / 2 / h[0];
                } else if (x==m0-1) {
                    gradT1 = (T[idx] - T[idx-1]) / 2 / h[0];
                    gradR1 = (R[idx] - R[idx-1]) / 2 / h[0];
                } else {
                    gradT1 = (T[idx+1] - T[idx-1]) / 2 / h[0];
                    gradR1 = (R[idx+1] - R[idx-1]) / 2 / h[0];
                }
                if (y==0) {
                    gradT2 = (T[idx+m0] - T[idx]) / 2 / h[1];
                    gradR2 = (R[idx+m0] - R[idx]) / 2 / h[1];
                } else if (y==m1-1) {
                    gradT2 = (T[idx] - T[idx-m0]) / 2 / h[1];
                    gradR2 = (R[idx] - R[idx-m0]) / 2 / h[1];
                } else {
                    gradT2 = (T[idx+m0] - T[idx-m0]) / 2 / h[1];
                    gradR2 = (R[idx+m0] - R[idx-m0]) / 2 / h[1];
                }
                if (z==0) {
                    gradT3 = (T[idx+m0*m1] - T[idx]) / 2 / h[2];
                    gradR3 = (R[idx+m0*m1] - R[idx]) / 2 / h[2];
                } else if (z==m2-1) {
                    gradT3 = (T[idx] - T[idx-m0*m1]) / 2 / h[2];
                    gradR3 = (R[idx] - R[idx-m0*m1]) / 2 / h[2];
                } else {
                    gradT3 = (T[idx+m0*m1] - T[idx-m0*m1]) / 2 / h[2];
                    gradR3 = (R[idx+m0*m1] - R[idx-m0*m1]) / 2 / h[2];
                }
                
                lengthGT = sqrt(gradT1 * gradT1 + gradT2  * gradT2
                              + gradT3 * gradT3 + edge[0] * edge[0]);
                lengthGR = sqrt(gradR1 * gradR1 + gradR2  * gradR2
                              + gradR3 * gradR3 + edge[0] * edge[0]);
                r1 = gradR1 * gradT1 + gradR2 * gradT2 + gradR3 * gradT3;
                r2 = 1 / (lengthGT * lengthGR);
                if (doDerivative) {
                    fac[idx]     = (r2 * gradR1     + r1
                            * (-1) / (lengthGR * lengthGT
                            * lengthGT * lengthGT) * gradT1)
                            / 2 / h[0];
                    fac[idx+N]   = (r2 * gradR2   + r1
                            * (-1) / (lengthGR * lengthGT
                            * lengthGT * lengthGT) * gradT2)
                            / 2 / h[1];
                    fac[idx+2*N] = (r2 * gradR3 + r1
                            * (-1) / (lengthGR * lengthGT
                            * lengthGT * lengthGT) * gradT3)
                            / 2 / h[2];
                }
                rc[idx] = r1 * r2;
                Dc[0]  += rc[idx]*rc[idx];
            }
        }
    }
    Dc[0] = o[0] - hd * Dc[0];
    
    
    if (!doDerivative) {
        return;
    }
    k = 0;
    for (int z=0;z<m2;z++) {
        for (int y=0;y<m1;y++) {
            for (int x=0;x<m0;x++) {
                idx = x+y*m0+z*m0*m1;
                drcJ[idx] = k;
                dD[idx] = 0;
                if (z==0) {
                    if (y==0) {
                        if (x==0) {
                            drcI[k] = idx;
                            drcP[k] = - fac[idx] - fac[idx+N] - fac[idx+2*N];
                            dD[idx] += rc[drcI[k]] * drcP[k];
                            k++;
                            drcI[k] = idx+1;
                            drcP[k] = - fac[idx+1];
                            dD[idx] += rc[drcI[k]] * drcP[k];
                            k++;
                        } else if (x==m0-1) {
                            drcI[k] = idx-1;
                            drcP[k] = fac[idx-1];
                            dD[idx] += rc[drcI[k]] * drcP[k];
                            k++;
                            drcI[k] = idx;
                            drcP[k] = fac[idx] - fac[idx+N] - fac[idx+2*N];
                            dD[idx] += rc[drcI[k]] * drcP[k];
                            k++;
                        } else {
                            drcI[k] = idx-1;
                            drcP[k] = fac[idx-1];
                            dD[idx] += rc[drcI[k]] * drcP[k];
                            k++;
                            drcI[k] = idx;
                            drcP[k] = - fac[idx+N] - fac[idx+2*N];
                            dD[idx] += rc[drcI[k]] * drcP[k];
                            k++;
                            drcI[k] = idx+1;
                            drcP[k] = - fac[idx+1];
                            dD[idx] += rc[drcI[k]] * drcP[k];
                            k++;
                        }
                        drcI[k] = idx+m0;
                        drcP[k] = - fac[idx+N+m0];
                        dD[idx] += rc[drcI[k]] * drcP[k];
                        k++;
                    } else if (y==m1-1) {
                        drcI[k] = idx-m0;
                        drcP[k] = fac[idx+N-m0];
                        dD[idx] += rc[drcI[k]] * drcP[k];
                        k++;
                        if (x==0) {
                            drcI[k] = idx;
                            drcP[k] = - fac[idx] + fac[idx+N] - fac[idx+2*N];
                            dD[idx] += rc[drcI[k]] * drcP[k];
                            k++;
                            drcI[k] = idx+1;
                            drcP[k] = - fac[idx+1];
                            dD[idx] += rc[drcI[k]] * drcP[k];
                            k++;
                        } else if (x==m0-1) {
                            drcI[k] = idx-1;
                            drcP[k] = fac[idx-1];
                            dD[idx] += rc[drcI[k]] * drcP[k];
                            k++;
                            drcI[k] = idx;
                            drcP[k] = fac[idx] + fac[idx+N] - fac[idx+2*N];
                            dD[idx] += rc[drcI[k]] * drcP[k];
                            k++;
                        } else {
                            drcI[k] = idx-1;
                            drcP[k] = fac[idx-1];
                            dD[idx] += rc[drcI[k]] * drcP[k];
                            k++;
                            drcI[k] = idx;
                            drcP[k] = fac[idx+N] - fac[idx+2*N];
                            dD[idx] += rc[drcI[k]] * drcP[k];
                            k++;
                            drcI[k] = idx+1;
                            drcP[k] = - fac[idx+1];
                            dD[idx] += rc[drcI[k]] * drcP[k];
                            k++;
                        }
                    } else {
                        drcI[k] = idx-m0;
                        drcP[k] = fac[idx+N-m0];
                        dD[idx] += rc[drcI[k]] * drcP[k];
                        k++;
                        if (x==0) {
                            drcI[k] = idx;
                            drcP[k] = - fac[idx] - fac[idx+2*N];
                            dD[idx] += rc[drcI[k]] * drcP[k];
                            k++;
                            drcI[k] = idx+1;
                            drcP[k] = - fac[idx+1];
                            dD[idx] += rc[drcI[k]] * drcP[k];
                            k++;
                        } else if (x==m0-1) {
                            drcI[k] = idx-1;
                            drcP[k] = fac[idx-1];
                            dD[idx] += rc[drcI[k]] * drcP[k];
                            k++;
                            drcI[k] = idx;
                            drcP[k] = fac[idx] - fac[idx+2*N];
                            dD[idx] += rc[drcI[k]] * drcP[k];
                            k++;
                        } else {
                            drcI[k] = idx-1;
                            drcP[k] = fac[idx-1];
                            dD[idx] += rc[drcI[k]] * drcP[k];
                            k++;
                            drcI[k] = idx;
                            drcP[k] = - fac[idx+2*N];
                            dD[idx] += rc[drcI[k]] * drcP[k];
                            k++;
                            drcI[k] = idx+1;
                            drcP[k] = - fac[idx+1];
                            dD[idx] += rc[drcI[k]] * drcP[k];
                            k++;
                        }
                        drcI[k] = idx+m0;
                        drcP[k] = - fac[idx+N+m0];
                        dD[idx] += rc[drcI[k]] * drcP[k];
                        k++;
                    }
                    drcI[k] = idx+m0*m1;
                    drcP[k] = - fac[idx+2*N+m0*m1];
                    dD[idx] += rc[drcI[k]] * drcP[k];
                    k++;
                } else if (z==m2-1) {
                    drcI[k] = idx-m0*m1;
                    drcP[k] = fac[idx+2*N-m0*m1];
                    dD[idx] += rc[drcI[k]] * drcP[k];
                    k++;
                    if (y==0) {
                        if (x==0) {
                            drcI[k] = idx;
                            drcP[k] = - fac[idx] - fac[idx+N] + fac[idx+2*N];
                            dD[idx] += rc[drcI[k]] * drcP[k];
                            k++;
                            drcI[k] = idx+1;
                            drcP[k] = - fac[idx+1];
                            dD[idx] += rc[drcI[k]] * drcP[k];
                            k++;
                        } else if (x==m0-1) {
                            drcI[k] = idx-1;
                            drcP[k] = fac[idx-1];
                            dD[idx] += rc[drcI[k]] * drcP[k];
                            k++;
                            drcI[k] = idx;
                            drcP[k] = fac[idx] - fac[idx+N] + fac[idx+2*N];
                            dD[idx] += rc[drcI[k]] * drcP[k];
                            k++;
                        } else {
                            drcI[k] = idx-1;
                            drcP[k] = fac[idx-1];
                            dD[idx] += rc[drcI[k]] * drcP[k];
                            k++;
                            drcI[k] = idx;
                            drcP[k] = - fac[idx+N] + fac[idx+2*N];
                            dD[idx] += rc[drcI[k]] * drcP[k];
                            k++;
                            drcI[k] = idx+1;
                            drcP[k] = - fac[idx+1];
                            dD[idx] += rc[drcI[k]] * drcP[k];
                            k++;
                        }
                        drcI[k] = idx+m0;
                        drcP[k] = - fac[idx+N+m0];
                        dD[idx] += rc[drcI[k]] * drcP[k];
                        k++;
                    } else if (y==m1-1) {
                        drcI[k] = idx-m0;
                        drcP[k] = fac[idx+N-m0];
                        dD[idx] += rc[drcI[k]] * drcP[k];
                        k++;
                        if (x==0) {
                            drcI[k] = idx;
                            drcP[k] = - fac[idx] + fac[idx+N] + fac[idx+2*N];
                            dD[idx] += rc[drcI[k]] * drcP[k];
                            k++;
                            drcI[k] = idx+1;
                            drcP[k] = - fac[idx+1];
                            dD[idx] += rc[drcI[k]] * drcP[k];
                            k++;
                        } else if (x==m0-1) {
                            drcI[k] = idx-1;
                            drcP[k] = fac[idx-1];
                            dD[idx] += rc[drcI[k]] * drcP[k];
                            k++;
                            drcI[k] = idx;
                            drcP[k] = fac[idx] + fac[idx+N] + fac[idx+2*N];
                            dD[idx] += rc[drcI[k]] * drcP[k];
                            k++;
                        } else {
                            drcI[k] = idx-1;
                            drcP[k] = fac[idx-1];
                            dD[idx] += rc[drcI[k]] * drcP[k];
                            k++;
                            drcI[k] = idx;
                            drcP[k] = fac[idx+N] + fac[idx+2*N];
                            dD[idx] += rc[drcI[k]] * drcP[k];
                            k++;
                            drcI[k] = idx+1;
                            drcP[k] = - fac[idx+1];
                            dD[idx] += rc[drcI[k]] * drcP[k];
                            k++;
                        }
                    } else {
                        drcI[k] = idx-m0;
                        drcP[k] = fac[idx+N-m0];
                        dD[idx] += rc[drcI[k]] * drcP[k];
                        k++;
                        if (x==0) {
                            drcI[k] = idx;
                            drcP[k] = - fac[idx] + fac[idx+2*N];
                            dD[idx] += rc[drcI[k]] * drcP[k];
                            k++;
                            drcI[k] = idx+1;
                            drcP[k] = - fac[idx+1];
                            dD[idx] += rc[drcI[k]] * drcP[k];
                            k++;
                        } else if (x==m0-1) {
                            drcI[k] = idx-1;
                            drcP[k] = fac[idx-1];
                            dD[idx] += rc[drcI[k]] * drcP[k];
                            k++;
                            drcI[k] = idx;
                            drcP[k] = fac[idx] + fac[idx+2*N];
                            dD[idx] += rc[drcI[k]] * drcP[k];
                            k++;
                        } else {
                            drcI[k] = idx-1;
                            drcP[k] = fac[idx-1];
                            dD[idx] += rc[drcI[k]] * drcP[k];
                            k++;
                            drcI[k] = idx;
                            drcP[k] = fac[idx+2*N];
                            dD[idx] += rc[drcI[k]] * drcP[k];
                            k++;
                            drcI[k] = idx+1;
                            drcP[k] = - fac[idx+1];
                            dD[idx] += rc[drcI[k]] * drcP[k];
                            k++;
                        }
                        drcI[k] = idx+m0;
                        drcP[k] = - fac[idx+N+m0];
                        dD[idx] += rc[drcI[k]] * drcP[k];
                        k++;
                    }
                } else {
                    drcI[k] = idx-m0*m1;
                    drcP[k] = fac[idx+2*N-m0*m1];
                    dD[idx] += rc[drcI[k]] * drcP[k];
                    k++;
                    if (y==0) {
                        if (x==0) {
                            drcI[k] = idx;
                            drcP[k] = - fac[idx] - fac[idx+N];
                            dD[idx] += rc[drcI[k]] * drcP[k];
                            k++;
                            drcI[k] = idx+1;
                            drcP[k] = - fac[idx+1];
                            dD[idx] += rc[drcI[k]] * drcP[k];
                            k++;
                        } else if (x==m0-1) {
                            drcI[k] = idx-1;
                            drcP[k] = fac[idx-1];
                            dD[idx] += rc[drcI[k]] * drcP[k];
                            k++;
                            drcI[k] = idx;
                            drcP[k] = fac[idx] - fac[idx+N];
                            dD[idx] += rc[drcI[k]] * drcP[k];
                            k++;
                        } else {
                            drcI[k] = idx-1;
                            drcP[k] = fac[idx-1];
                            dD[idx] += rc[drcI[k]] * drcP[k];
                            k++;
                            drcI[k] = idx;
                            drcP[k] = - fac[idx+N];
                            dD[idx] += rc[drcI[k]] * drcP[k];
                            k++;
                            drcI[k] = idx+1;
                            drcP[k] = - fac[idx+1];
                            dD[idx] += rc[drcI[k]] * drcP[k];
                            k++;
                        }
                        drcI[k] = idx+m0;
                        drcP[k] = - fac[idx+N+m0];
                        dD[idx] += rc[drcI[k]] * drcP[k];
                        k++;
                    } else if (y==m1-1) {
                        drcI[k] = idx-m0;
                        drcP[k] = fac[idx+N-m0];
                        dD[idx] += rc[drcI[k]] * drcP[k];
                        k++;
                        if (x==0) {
                            drcI[k] = idx;
                            drcP[k] = - fac[idx] + fac[idx+N];
                            dD[idx] += rc[drcI[k]] * drcP[k];
                            k++;
                            drcI[k] = idx+1;
                            drcP[k] = - fac[idx+1];
                            dD[idx] += rc[drcI[k]] * drcP[k];
                            k++;
                        } else if (x==m0-1) {
                            drcI[k] = idx-1;
                            drcP[k] = fac[idx-1];
                            dD[idx] += rc[drcI[k]] * drcP[k];
                            k++;
                            drcI[k] = idx;
                            drcP[k] = fac[idx] + fac[idx+N];
                            dD[idx] += rc[drcI[k]] * drcP[k];
                            k++;
                        } else {
                            drcI[k] = idx-1;
                            drcP[k] = fac[idx-1];
                            dD[idx] += rc[drcI[k]] * drcP[k];
                            k++;
                            drcI[k] = idx;
                            drcP[k] = fac[idx+N];
                            dD[idx] += rc[drcI[k]] * drcP[k];
                            k++;
                            drcI[k] = idx+1;
                            drcP[k] = - fac[idx+1];
                            dD[idx] += rc[drcI[k]] * drcP[k];
                            k++;
                        }
                    } else {
                        drcI[k] = idx-m0;
                        drcP[k] = fac[idx+N-m0];
                        dD[idx] += rc[drcI[k]] * drcP[k];
                        k++;
                        if (x==0) {
                            drcI[k] = idx;
                            drcP[k] = - fac[idx];
                            dD[idx] += rc[drcI[k]] * drcP[k];
                            k++;
                            drcI[k] = idx+1;
                            drcP[k] = - fac[idx+1];
                            dD[idx] += rc[drcI[k]] * drcP[k];
                            k++;
                        } else if (x==m0-1) {
                            drcI[k] = idx-1;
                            drcP[k] = fac[idx-1];
                            dD[idx] += rc[drcI[k]] * drcP[k];
                            k++;
                            drcI[k] = idx;
                            drcP[k] = fac[idx];
                            dD[idx] += rc[drcI[k]] * drcP[k];
                            k++;
                        } else {
                            drcI[k] = idx-1;
                            drcP[k] = fac[idx-1];
                            dD[idx] += rc[drcI[k]] * drcP[k];
                            k++;
                            drcI[k] = idx+1;
                            drcP[k] = - fac[idx+1];
                            dD[idx] += rc[drcI[k]] * drcP[k];
                            k++;
                        }
                        drcI[k] = idx+m0;
                        drcP[k] = - fac[idx+N+m0];
                        dD[idx] += rc[drcI[k]] * drcP[k];
                        k++;
                    }
                    drcI[k] = idx+m0*m1;
                    drcP[k] = - fac[idx+2*N+m0*m1];
                    dD[idx] += rc[drcI[k]] * drcP[k];
                    k++;
                }
                dD[idx] = -2 * hd * dD[idx];
            }
        }
    }
    drcJ[m0*m1*m2] = k;
    d2psi[0]       = 2 * hd;
}

/* ========================================================================= */
    