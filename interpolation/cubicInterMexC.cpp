/*
 * (c) Fabian Gigengack 2011/04/13 see FAIR.2 and FAIRcopyright.m.
 * http://www.uni-muenster.de/EIMI/
 *
 * CPP code for cubic interpolation. See cubicInter for details.
 */

#include <math.h>
#include <mex.h>
#ifdef _OPENMP
	#include <omp.h>
#endif

inline double  cint(double p0, double p1, double p2, double p3, double v);
inline double dcint(double p0, double p1, double p2, double p3, double v);
void cubicInter1D(double *Tc, double *dT, const double *T, const double *omega, const double *m, int N, const double *X, double *h, bool doDerivative);
void cubicInter2D(double *Tc, double *dT, const double *T, const double *omega, const double *m, int N, const double *X, double *h, bool doDerivative);
void cubicInter3D(double *Tc, double *dT, const double *T, const double *omega, const double *m, int N, const double *X, double *h, bool doDerivative);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    int dim, i;
    if (nrhs<5)
        mexErrMsgTxt("Number of arguments must be 5!");
    
    // Get input
    const double* T     = static_cast<double*>(mxGetData(prhs[0]));
    const double* omega = static_cast<double*>(mxGetData(prhs[1]));
    const double* m     = static_cast<double*>(mxGetData(prhs[2]));
    const double* X     = static_cast<double*>(mxGetData(prhs[3]));
    bool doDerivative   = (bool) *mxGetLogicals(prhs[4]);
    
    // get dimension, h and number of elements N
    dim = mxGetN(prhs[1]) / 2;
    double h[3]; //h[dim];
    for (i=0; i<dim; i++) {
        h[i] = (omega[2*i+1]-omega[2*i]) / m[i];
    }
    const int N = mxGetM(prhs[3])/dim;
    
    //Allocate Tc
    mwSize dims[2]; dims[0] = N; dims[1] = 1;
    plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
    double* Tc = static_cast<double*>(mxGetData(plhs[0]));
    
    double* dT = 0;
    //Allocate dT
    if (doDerivative) {
        dims[1] = dim;
        plhs[1] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
        dT = static_cast<double*>(mxGetData(plhs[1]));
    } else {
        dims[0] = 1; dims[1] = 1;
        plhs[1] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
        dT = static_cast<double*>(mxGetData(plhs[1]));
    }
    
    switch (dim) {
        case 1:
            cubicInter1D(Tc, dT, T, omega, m, N, X, h, doDerivative);
            break;
        case 2:
            cubicInter2D(Tc, dT, T, omega, m, N, X, h, doDerivative);
            break;
        case 3:
            cubicInter3D(Tc, dT, T, omega, m, N, X, h, doDerivative);
            break;
        default:
            break;
    }
    
}

void cubicInter1D(double *Tc, double *dT, const double *T,
        const double *omega, const double *m, int N, const double *X,
        double *h, bool doDerivative) {
    // Increments in X and Y direction
    double x, p[4];
    int i, l, xf;
    
    #pragma omp parallel for default(shared) private(i, l, x, xf, p)
    for (i=0; i<N; i++) {
        // map x from [h/2,omega-h/2] -> [0,m-1],
        x = (X[i] - omega[0]) / h[0] + .5 - 1;
        
        // Valid?
        if (x<-1 || x>=m[0])
            continue;
        
        if (x>=1 && x<m[0]-2) {
            // No boundary treatment
            xf = floor(x);
            // Extract remainder
            x = x - xf;
            
            Tc[i] = cint(T[xf-1], T[xf], T[xf+1], T[xf+2], x);

            if (mxIsNaN(Tc[i]))
                mexPrintf("Is NAN");

            if (doDerivative) {
                dT[i] = dcint(T[xf-1], T[xf], T[xf+1], T[xf+2], x) / h[0];
            }
        } else {
            // Boundary treatment
            xf = floor(x);
            // Extract remainder
            x = x - xf;

            for (l=-1;l<3;l++) {
                p[l+1] = (xf+l<0 || xf+l>m[0]-1)? 0: T[xf+l];
            }
            Tc[i] = cint(p[0], p[1], p[2], p[3], x);

            if (mxIsNaN(Tc[i]))
                mexPrintf("Is NAN");

            if (doDerivative) {
                dT[i] = dcint(p[0], p[1], p[2], p[3], x) / h[0];
            }
        }
    }
}

void cubicInter2D(double *Tc, double *dT, const double *T,
        const double *omega, const double *m, int N, const double *X,
        double *h, bool doDerivative) {
    // Increments in X and Y direction
    double x, y, u[4], p[4];
    int i, j, l, xf, yf, i2 = (int)m[0], idx;
    
    #pragma omp parallel for default(shared) private(i, j, l, x, y, xf, yf, p, u, idx)
    for (i=0; i<N; i++) {
        // map x from [h/2,omega-h/2] -> [0,m-1],
        x = (X[i]   - omega[0]) / h[0] + .5 - 1;
        y = (X[i+N] - omega[2]) / h[1] + .5 - 1;
        
        // Valid?
        if (x<-1 || y<-1 || x>=m[0] || y>=m[1])
            continue;
        
        if (x>=1 && y>=1 && x<m[0]-2 && y<m[1]-2) {
            // No boundary treatment
            xf = floor(x);
            yf = floor(y);
            // Extract remainder
            x = x - xf;
            y = y - yf;
            for (j=-1; j<3; j++) {
                idx = xf+j+i2*yf;
                u[j+1] = cint(T[idx-i2], T[idx], T[idx+i2], T[idx+2*i2], y);
            }
            Tc[i] = cint(u[0], u[1], u[2], u[3], x);
            
            if (mxIsNaN(Tc[i]))
                mexPrintf("Is NAN");
            
            if (doDerivative) {
                dT[i] = dcint(u[0], u[1], u[2], u[3], x) / h[0];
                for (j=-1; j<3; j++) {
                    idx = xf+j+i2*yf;
                    u[j+1] = dcint(T[idx-i2], T[idx], T[idx+i2], T[idx+2*i2], y);
                }
                dT[i+N] = cint(u[0], u[1], u[2], u[3], x) / h[1];
            }
        } else {
            // Boundary treatment
            xf = floor(x);
            yf = floor(y);
            // Extract remainder
            x = x - xf;
            y = y - yf;
            for (j=-1; j<3; j++) {
                idx = xf+j+i2*yf;
                for (l=-1;l<3;l++) {
                    p[l+1] = (xf+j<0 || xf+j>m[0]-1 || yf+l<0 || yf+l>m[1]-1)? 0: T[idx+i2*l];
                }
                u[j+1] = cint(p[0], p[1], p[2], p[3], y);
            }
            Tc[i] = cint(u[0], u[1], u[2], u[3], x);
            
            if (mxIsNaN(Tc[i]))
                mexPrintf("Is NAN");
            
            if (doDerivative) {
                dT[i] = dcint(u[0], u[1], u[2], u[3], x) / h[0];
                for (j=-1; j<3; j++) {
                    idx = xf+j+i2*yf;
                    for (l=-1;l<3;l++) {
                        p[l+1] = (xf+j<0 || xf+j>m[0]-1 || yf+l<0 || yf+l>m[1]-1)? 0: T[idx+i2*l];
                    }
                    u[j+1] = dcint(p[0], p[1], p[2], p[3], y);
                }
                dT[i+N] = cint(u[0], u[1], u[2], u[3], x) / h[1];
            }
        }
    }
}

void cubicInter3D(double *Tc, double *dT, const double *T,
        const double *omega, const double *m, int N, const double *X,
        double *h, bool doDerivative) {
    // Increments in X,Y and Z direction
    double x, y, z, u[4], t[16], p[4];
    int i, j, k, l, xf, yf, zf, i2 = (int)m[0], i3 = i2*(int)m[1], idx;
    
    #pragma omp parallel for default(shared) private(i, j, k, l, x, y, z, xf, yf, zf, p, u, t, idx)
    for (i=0;i<N;i++) {
        // map x from [h/2,omega-h/2] -> [0,m-1],
        x = (X[i]     - omega[0]) / h[0] + .5 - 1;
        y = (X[i+N]   - omega[2]) / h[1] + .5 - 1;
        z = (X[i+2*N] - omega[4]) / h[2] + .5 - 1;
        
        // Valid?
        if (x<-1 || y<-1 || z<-1 || x>=m[0] || y>=m[1] || z>=m[2])
            continue;
        
        if (x>=1 && y>=1 && z>=1 && x<m[0]-2 && y<m[1]-2 && z<m[2]-2) {
            // No boundary treatment
            xf = floor(x);
            yf = floor(y);
            zf = floor(z);

            // Extract remainder
            x = x - xf;
            y = y - yf;
            z = z - zf;

            // cubic interpolation in z-dimension
            for (j=-1; j<3; j++) {
                for (k=-1; k<3; k++) {
                    idx = xf+j+i2*(yf+k)+i3*zf;
                    t[j+1+(k+1)*4] = cint(T[idx-i3], T[idx], T[idx+i3], T[idx+2*i3], z);
                }
            }
            // cubic interpolation in y-dimension
            for (j=0; j<4; j++) {
                u[j] = cint(t[j], t[j+4], t[j+8], t[j+12], y);
            }
            // cubic interpolation in x-dimension
            Tc[i] = cint(u[0], u[1], u[2], u[3], x);


            if (mxIsNaN(Tc[i]))
                mexPrintf("Is NAN");

            if (doDerivative) {
                // derivative of cubic interpolation in x-direction
                dT[i] = dcint(u[0], u[1], u[2], u[3], x) / h[0];

                // derivative of cubic interpolation in y-direction
                for (j=0; j<4; j++) {
                    u[j] = dcint(t[j], t[j+4], t[j+8], t[j+12], y);
                }
                dT[i+N] = cint(u[0], u[1], u[2], u[3], x) / h[1];

                // derivative of cubic interpolation in z-direction
                for (j=-1; j<3; j++) {
                    for (k=-1; k<3; k++) {
                        idx = xf+j+i2*(yf+k)+i3*zf;
                        t[j+1+(k+1)*4] = dcint(T[idx-i3], T[idx], T[idx+i3], T[idx+2*i3], z);
                    }
                }
                for (j=0; j<4; j++) {
                    u[j] = cint(t[j], t[j+4], t[j+8], t[j+12], y);
                }
                dT[i+2*N] = cint(u[0], u[1], u[2], u[3], x) / h[2];
            }
        } else {
            // Boundary treatment
            xf = floor(x);
            yf = floor(y);
            zf = floor(z);

            // Extract remainder
            x = x - xf;
            y = y - yf;
            z = z - zf;

            // cubic interpolation in z-dimension
            for (j=-1; j<3; j++) {
                for (k=-1; k<3; k++) {
                    idx = xf+j+i2*(yf+k)+i3*zf;
                    for (l=-1;l<3;l++) {
                        p[l+1] = (xf+j<0 || xf+j>m[0]-1 || yf+k<0 || yf+k>m[1]-1 || zf+l<0 || zf+l>m[2]-1)? 0: T[idx+i3*l];
                    }
                    t[j+1+(k+1)*4] = cint(p[0], p[1], p[2], p[3], z);
                }
            }
            // cubic interpolation in y-dimension
            for (j=0; j<4; j++) {
                u[j] = cint(t[j], t[j+4], t[j+8], t[j+12], y);
            }
            // cubic interpolation in x-dimension
            Tc[i] = cint(u[0], u[1], u[2], u[3], x);


            if (mxIsNaN(Tc[i]))
                mexPrintf("Is NAN");

            if (doDerivative) {
                // derivative of cubic interpolation in x-direction
                dT[i] = dcint(u[0], u[1], u[2], u[3], x) / h[0];

                // derivative of cubic interpolation in y-direction
                for (j=0; j<4; j++) {
                    u[j] = dcint(t[j], t[j+4], t[j+8], t[j+12], y);
                }
                dT[i+N] = cint(u[0], u[1], u[2], u[3], x) / h[1];

                // derivative of cubic interpolation in z-direction
                for (j=-1; j<3; j++) {
                    for (k=-1; k<3; k++) {
                        idx = xf+j+i2*(yf+k)+i3*zf;
                        for (l=-1;l<3;l++) {
                            p[l+1] = (xf+j<0 || xf+j>m[0]-1 || yf+k<0 || yf+k>m[1]-1 || zf+l<0 || zf+l>m[2]-1)? 0: T[idx+i3*l];
                        }
                        t[j+1+(k+1)*4] = dcint(p[0], p[1], p[2], p[3], z);
                    }
                }
                for (j=0; j<4; j++) {
                    u[j] = cint(t[j], t[j+4], t[j+8], t[j+12], y);
                }
                dT[i+2*N] = cint(u[0], u[1], u[2], u[3], x) / h[2];
            }
        }
    }
}

inline double cint(double p0, double p1, double p2, double p3, double v) {
    double v2 = v * v;
    double v3 = v * v * v;
    // Cubic interpolation function
    return    (- .5 * v3 +       v2 - .5 * v    ) * p0
            + ( 1.5 * v3 - 2.5 * v2          + 1) * p1
            + (-1.5 * v3 + 2   * v2 + .5 * v    ) * p2
            + (  .5 * v3 -  .5 * v2             ) * p3;
}

inline double dcint(double p0, double p1, double p2, double p3, double v) {
    double v2 = v * v;
    // Derivative of cubic interpolation function
    return    (-1.5 * v2 + 2 * v - .5) * p0
            + ( 4.5 * v2 - 5 * v     ) * p1
            + (-4.5 * v2 + 4 * v + .5) * p2
            + ( 1.5 * v2 -     v     ) * p3;
}

/*
 * =======================================================================================
 * FAIR: Flexible Algorithms for Image Registration, Version 2011
 * Copyright (c): Jan Modersitzki
 * Maria-Goeppert-Str. 1a, D-23562 Luebeck, Germany
 * Email: jan.modersitzki@mic.uni-luebeck.de
 * URL:   http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
 * =======================================================================================
 * No part of this code may be reproduced, stored in a retrieval system,
 * translated, transcribed, transmitted, or distributed in any form
 * or by any means, means, manual, electric, electronic, electro-magnetic,
 * mechanical, chemical, optical, photocopying, recording, or otherwise,
 * without the prior explicit written permission of the authors or their
 * designated proxies. In no event shall the above copyright notice be
 * removed or altered in any way.
 * 
 * This code is provided "as is", without any warranty of any kind, either
 * expressed or implied, including but not limited to, any implied warranty
 * of merchantibility or fitness for any purpose. In no event will any party
 * who distributed the code be liable for damages or for any claim(s) by
 * any other party, including but not limited to, any lost profits, lost
 * monies, lost data or data rendered inaccurate, losses sustained by
 * third parties, or any other special, incidental or consequential damages
 * arrising out of the use or inability to use the program, even if the
 * possibility of such damages has been advised against. The entire risk
 * as to the quality, the performace, and the fitness of the program for any
 * particular purpose lies with the party using the code.
 * =======================================================================================
 * Any use of this code constitutes acceptance of the terms of the above statements
 * =======================================================================================
 */