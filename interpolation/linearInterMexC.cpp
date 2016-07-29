/*
 * (c) Jan Modersitzki and Fabian Gigengack 2011/04/20, see FAIR.2 and FAIRcopyright.m.
 * http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
 * http://www.uni-muenster.de/EIMI/
 *
 * CPP code for linear interpolation. See linearInter for details.
 */

#include <math.h>
#include <mex.h>
#ifdef _OPENMP
	#include <omp.h>
#endif

void linearInter1D(double *Tc, double *dT, const double *T, const double *omega, const double *m, int N, const double *X, double *h, bool doDerivative, bool boundary);
void linearInter2D(double *Tc, double *dT, const double *T, const double *omega, const double *m, int N, const double *X, double *h, bool doDerivative, bool boundary);
void linearInter3D(double *Tc, double *dT, const double *T, const double *omega, const double *m, int N, const double *X, double *h, bool doDerivative, bool boundary);
inline int min(int a, int b) {return (a<b)?a:b;}
inline int max(int a, int b) {return (a>b)?a:b;}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    int dim, i;
    if (nrhs<6)
        mexErrMsgTxt("Number of arguments must be 6!");
    
    // Get input
    const double* T     = static_cast<double*>(mxGetData(prhs[0]));
    const double* omega = static_cast<double*>(mxGetData(prhs[1]));
    const double* m     = static_cast<double*>(mxGetData(prhs[2]));
    const double* X     = static_cast<double*>(mxGetData(prhs[3]));
    bool doDerivative   = (bool) *mxGetLogicals(prhs[4]);
    bool boundary       = (bool) *mxGetLogicals(prhs[5]);
    
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
            linearInter1D(Tc, dT, T, omega, m, N, X, h, doDerivative, boundary);
            break;
        case 2:
            linearInter2D(Tc, dT, T, omega, m, N, X, h, doDerivative, boundary);
            break;
        case 3:
            linearInter3D(Tc, dT, T, omega, m, N, X, h, doDerivative, boundary);
            break;
        default:
            break;
    }
    
}

void linearInter1D(double *Tc, double *dT, const double *T,
        const double *omega, const double *m, int N, const double *X,
        double *h, bool doDerivative, bool boundary) {
    // Increments in X and Y direction
    double x, p[2];
    int i, xf;
    
    #pragma omp parallel for default(shared) private(i, x, xf, p)
    for (i=0; i<N; i++) {
        // map x from [h/2,omega-h/2] -> [0,m-1],
        x = (X[i] - omega[0]) / h[0] + .5 - 1;
        
        xf = floor(x);
        
        if (boundary) {
            // Replicate boundary
            p[0] = T[min(m[0]-1,max(0,xf))];
            p[1] = T[min(m[0]-1,max(0,xf+1))];
        } else {
            // Zero padding
            // Valid?
            if (x<-1 || x>=m[0])
                continue;
            
            p[0] = (xf<0)?        0: T[xf];
            p[1] = (xf+1>m[0]-1)? 0: T[xf+1];
        }
        
        // Extract remainder
        x = x - xf;
        
        Tc[i] = p[0] * (1-x) + p[1] * x;
        
        if (mxIsNaN(Tc[i]))
            mexPrintf("Is NAN");
        
        if (doDerivative) {
            dT[i] = (p[1] - p[0]) / h[0];
        }
    }
}

void linearInter2D(double *Tc, double *dT, const double *T,
        const double *omega, const double *m, int N, const double *X,
        double *h, bool doDerivative, bool boundary) {
    // Increments in X and Y direction
    double x, y, p[4];
    int i, xf, yf, i2 = m[0];
    int x1, x2, y1, y2;
    
    #pragma omp parallel for default(shared) private(i, x, y, xf, yf, p, x1, x2, y1, y2)
    for (i=0; i<N; i++) {
        // map x from [h/2,omega-h/2] -> [0,m-1],
        x = (X[i]   - omega[0]) / h[0] + .5 - 1;
        y = (X[i+N] - omega[2]) / h[1] + .5 - 1;
        
        xf = floor(x);
        yf = floor(y);
        
        if (boundary) {
            // Replicate boundary
            x1 = min(m[0]-1,max(0,xf));
            x2 = min(m[0]-1,max(0,xf+1));
            y1 = min(m[1]-1,max(0,yf));
            y2 = min(m[1]-1,max(0,yf+1));

            p[0] = T[x1 + i2 * y1];
            p[1] = T[x2 + i2 * y1];
            p[2] = T[x1 + i2 * y2];
            p[3] = T[x2 + i2 * y2];
        } else {
            // Zero padding
            // Valid?
            if (x<-1 || y<-1 || x>=m[0] || y>=m[1])
                continue;
            
            p[0] = (xf<0        || yf<0)?        0: T[xf  +i2*yf];
            p[1] = (xf+1>m[0]-1 || yf<0)?        0: T[xf+1+i2*yf];
            p[2] = (xf<0        || yf+1>m[1]-1)? 0: T[xf  +i2*(yf+1)];
            p[3] = (xf+1>m[0]-1 || yf+1>m[1]-1)? 0: T[xf+1+i2*(yf+1)];
        }
        
        // Extract remainder
        x = x - xf;
        y = y - yf;
        
        Tc[i] = (p[0] * (1-x) + p[1] * x) * (1-y) + (p[2] * (1-x) + p[3] * x) * y;
        
        if (mxIsNaN(Tc[i]))
            mexPrintf("Is NAN");
        
        if (doDerivative) {
            dT[i]   = ((p[1] - p[0]) * (1-y) + (p[3] - p[2]) * y) / h[0];
            dT[i+N] = ((p[2] - p[0]) * (1-x) + (p[3] - p[1]) * x) / h[1];
        }
    }
}

void linearInter3D(double *Tc, double *dT, const double *T,
        const double *omega, const double *m, int N, const double *X,
        double *h, bool doDerivative, bool boundary) {
    // Increments in X,Y and Z direction
    double x, y, z, p[8];
    int i, xf, yf, zf;
    int i2 = m[0], i3 = i2*m[1];
    int x1, x2, y1, y2, z1, z2;
    
    #pragma omp parallel for default(shared) private(i, x, y, z, xf, yf, zf, p, x1, x2, y1, y2, z1, z2)
    for (i=0;i<N;i++) {
        // map x from [h/2,omega-h/2] -> [0,m-1],
        x = (X[i]     - omega[0]) / h[0] + .5 - 1;
        y = (X[i+N]   - omega[2]) / h[1] + .5 - 1;
        z = (X[i+2*N] - omega[4]) / h[2] + .5 - 1;
        
        xf = floor(x);
        yf = floor(y);
        zf = floor(z);
        
        if (boundary) {
            // Replicate boundary
            x1 = min(m[0]-1,max(0,xf));
            x2 = min(m[0]-1,max(0,xf+1));
            y1 = min(m[1]-1,max(0,yf));
            y2 = min(m[1]-1,max(0,yf+1));
            z1 = min(m[2]-1,max(0,zf));
            z2 = min(m[2]-1,max(0,zf+1));

            p[0] = T[x1 + i2 * y1 + i3 * z1];
            p[1] = T[x2 + i2 * y1 + i3 * z1];
            p[2] = T[x1 + i2 * y2 + i3 * z1];
            p[3] = T[x2 + i2 * y2 + i3 * z1];
            p[4] = T[x1 + i2 * y1 + i3 * z2];
            p[5] = T[x2 + i2 * y1 + i3 * z2];
            p[6] = T[x1 + i2 * y2 + i3 * z2];
            p[7] = T[x2 + i2 * y2 + i3 * z2];
        } else {
            // Zero padding
            // Valid?
            if (x<-1 || y<-1 || z<-1 || x>=m[0] || y>=m[1] || z>=m[2])
                continue;
            
            p[0] = (xf<0        || yf<0        || zf<0)?        0: T[xf  +i2*yf    +i3*zf];
            p[1] = (xf+1>m[0]-1 || yf<0        || zf<0)?        0: T[xf+1+i2*yf    +i3*zf];
            p[2] = (xf<0        || yf+1>m[1]-1 || zf<0)?        0: T[xf  +i2*(yf+1)+i3*zf];
            p[3] = (xf+1>m[0]-1 || yf+1>m[1]-1 || zf<0)?        0: T[xf+1+i2*(yf+1)+i3*zf];
            p[4] = (xf<0        || yf<0        || zf+1>m[2]-1)? 0: T[xf  +i2*yf    +i3*(zf+1)];
            p[5] = (xf+1>m[0]-1 || yf<0        || zf+1>m[2]-1)? 0: T[xf+1+i2*yf    +i3*(zf+1)];
            p[6] = (xf<0        || yf+1>m[1]-1 || zf+1>m[2]-1)? 0: T[xf  +i2*(yf+1)+i3*(zf+1)];
            p[7] = (xf+1>m[0]-1 || yf+1>m[1]-1 || zf+1>m[2]-1)? 0: T[xf+1+i2*(yf+1)+i3*(zf+1)];
        }

        // Extract remainder
        x = x - xf;
        y = y - yf;
        z = z - zf;
        
        Tc[i] = ((p[0] * (1-x) + p[1] * x) * (1-y)
              +  (p[2] * (1-x) + p[3] * x) *    y) * (1-z)
              + ((p[4] * (1-x) + p[5] * x) * (1-y)
              +  (p[6] * (1-x) + p[7] * x) *    y) *    z;
        
        if (mxIsNaN(Tc[i]))
            mexPrintf("Is NAN");
        
        if (doDerivative) {
            dT[i]     = (((p[1] - p[0]) * (1-y) + (p[3] - p[2]) * y) * (1-z)
                      +  ((p[5] - p[4]) * (1-y) + (p[7] - p[6]) * y) *    z) / h[0];
            
            dT[i+N]   = (((p[2] - p[0]) * (1-x) + (p[3] - p[1]) * x) * (1-z)
                      +  ((p[6] - p[4]) * (1-x) + (p[7] - p[5]) * x) *    z) / h[1];
            
            dT[i+2*N] = (((p[4] * (1-x) + p[5] * x) * (1-y)
                      +   (p[6] * (1-x) + p[7] * x) *    y)
                      -  ((p[0] * (1-x) + p[1] * x) * (1-y)
                      +   (p[2] * (1-x) + p[3] * x) *    y)) / h[2];
        }
    }
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