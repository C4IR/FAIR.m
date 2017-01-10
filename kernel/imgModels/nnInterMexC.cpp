/*
 * (c) Fabian Gigengack 2011/08/12, see FAIR.2 and FAIRcopyright.m.
 * http://www.uni-muenster.de/EIMI/
 *
 * CPP code for next neighbor interpolation. See nnInter for details.
 */

// #include <math.h>
#include <mex.h>
#ifdef _OPENMP
    #include <omp.h>
#endif

void nnInter1D(double *Tc, const double *T, const double *omega, const double *m, int N, const double *X, double *h, bool boundary);
void nnInter2D(double *Tc, const double *T, const double *omega, const double *m, int N, const double *X, double *h, bool boundary);
void nnInter3D(double *Tc, const double *T, const double *omega, const double *m, int N, const double *X, double *h, bool boundary);
inline int min(int a, int b) {return (a<b)?a:b;}
inline int max(int a, int b) {return (a>b)?a:b;}

//needed for windows with ms compiler
inline int round(double a) {return a<0?a-.5:a+.5;}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    int dim, i;
    if (nrhs<5)
        mexErrMsgTxt("Number of arguments must be 5!");
    
    // Get input
    const double* T     = static_cast<double*>(mxGetData(prhs[0]));
    const double* omega = static_cast<double*>(mxGetData(prhs[1]));
    
    const double* m     = static_cast<double*>(mxGetData(prhs[2]));
    const double* X     = static_cast<double*>(mxGetData(prhs[3]));
    bool boundary       = (bool) *mxGetLogicals(prhs[4]);
    
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
    
    switch (dim) {
        case 1:
            nnInter1D(Tc, T, omega, m, N, X, h, boundary);
            break;
        case 2:
            nnInter2D(Tc, T, omega, m, N, X, h, boundary);
            break;
        case 3:
            nnInter3D(Tc, T, omega, m, N, X, h, boundary);
            break;
        default:
            break;
    }
    
}

void nnInter1D(double *Tc, const double *T,
        const double *omega, const double *m, int N, const double *X,
        double *h, bool boundary) {
    double x;
    int i, xf;
    
    #pragma omp parallel for default(shared) private(i, x, xf)
    for (i=0; i<N; i++) {
        // map x from [h/2,omega-h/2] -> [0,m-1],
        x = (X[i] - omega[0]) / h[0] + .5 - 1;
        
        xf = round(x);
        
        if (boundary) {
            // Replicate boundary
            Tc[i] = T[min(m[0]-1,max(0,xf))];
        } else {
            // Zero padding
            Tc[i] = (xf<0 || xf>m[0]-1)? 0: T[xf];
        }

        if (mxIsNaN(Tc[i]))
            mexPrintf("Is NAN");
    }
}

void nnInter2D(double *Tc, const double *T,
        const double *omega, const double *m, int N, const double *X,
        double *h, bool boundary) {
    // Increments in X and Y direction
    double x, y;
    int i, xf, yf, i2 = m[0];
    
    #pragma omp parallel for default(shared) private(i, x, y, xf, yf)
    for (i=0; i<N; i++) {
        // map x from [h/2,omega-h/2] -> [0,m-1],
        x = (X[i]   - omega[0]) / h[0] + .5 - 1;
        y = (X[i+N] - omega[2]) / h[1] + .5 - 1;
        
        xf = round(x);
        yf = round(y);
        
        if (boundary) {
            // Replicate boundary
            Tc[i] = T[min(m[0]-1,max(0,xf)) + i2 * min(m[1]-1,max(0,yf))];
        } else {
            // Zero padding
            Tc[i] = (xf<0 || yf<0 || xf>m[0]-1 || yf>m[1]-1)? 0: T[xf+i2*yf];
        }

        if (mxIsNaN(Tc[i]))
            mexPrintf("Is NAN");
    }
}

void nnInter3D(double *Tc, const double *T,
        const double *omega, const double *m, int N, const double *X,
        double *h, bool boundary) {
    // Increments in X,Y and Z direction
    double x, y, z;
    int i, xf, yf, zf;
    int i2 = m[0], i3 = i2*m[1];
    
    #pragma omp parallel for default(shared) private(i, x, y, z, xf, yf, zf)
    for (i=0;i<N;i++) {
        // map x from [h/2,omega-h/2] -> [0,m-1],
        x = (X[i]     - omega[0]) / h[0] + .5 - 1;
        y = (X[i+N]   - omega[2]) / h[1] + .5 - 1;
        z = (X[i+2*N] - omega[4]) / h[2] + .5 - 1;
        
        xf = round(x);
        yf = round(y);
        zf = round(z);
        
        if (boundary) {
            // Replicate boundary
            Tc[i] = T[       min(m[0]-1,max(0,xf))
            + i2 * min(m[1]-1,max(0,yf))
            + i3 * min(m[2]-1,max(0,zf))];
        } else {
            // Zero padding
            Tc[i] = (   xf<0      || yf<0      || zf<0
                     || xf>m[0]-1 || yf>m[1]-1 || zf>m[2]-1)?
                        0: T[xf+i2*yf+i3*zf];
        }

        if (mxIsNaN(Tc[i]))
            mexPrintf("Is NAN");
    }
}

/*===================================================================================== */

