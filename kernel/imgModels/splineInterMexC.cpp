/*
 * (c) Jan Modersitzki and Fabian Gigengack 2011/02/02, see FAIR.2 and FAIRcopyright.m.
 * http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
 * http://www.uni-muenster.de/EIMI/
 *
 * CPP code for spline interpolation. See splineInter for details.
 */

#include <mex.h>
#include <math.h>
#ifdef _OPENMP
    #include <omp.h>
#endif

void splineInter1D(double *Tc,double *dT, const double *T, const double *omega, const double *m, int N, const double *X, double *h, bool doDerivative, bool boundary);
void splineInter2D(double *Tc,double *dT, const double *T, const double *omega, const double *m, int N, const double *X, double *h, bool doDerivative, bool boundary);
void splineInter3D(double *Tc,double *dT, const double *T, const double *omega, const double *m, int N, const double *X, double *h, bool doDerivative, bool boundary);
inline double  b0(int j, double xi);
inline double db0(int j, double xi);
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
    plhs[0] = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS, mxREAL);
    double* Tc = static_cast<double*>(mxGetData(plhs[0]));
    
    double* dT = 0;
    //Allocate dT
    if (doDerivative) {
        dims[1] = dim; 
        plhs[1] = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS, mxREAL);
        dT = static_cast<double*>(mxGetData(plhs[1]));
    } else {
        dims[0] = 1; dims[1] = 1;
        plhs[1] = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS, mxREAL);
        dT = static_cast<double*>(mxGetData(plhs[1]));
    }
    
    switch (dim) {
        case 1:
            splineInter1D(Tc,dT,T,omega,m,N,X,h,doDerivative, boundary);
            break;
        case 2:
            splineInter2D(Tc,dT,T,omega,m,N,X,h,doDerivative, boundary);
            break;
        case 3:
            splineInter3D(Tc,dT,T,omega,m,N,X,h,doDerivative, boundary);
            break;
        default:
            break;
    }
}

void splineInter1D(double *Tc, double *dT, const double *T,
        const double *omega, const double *m, int N, const double *X,
        double *h, bool doDerivative, bool boundary) {
    int xf, i, j, m0 = m[0];
    double x, p;
    
    #pragma omp parallel for default(shared) private(x, xf, p, i, j)
    for (i=0;i<N;i++) {
        Tc[i] = 0;
        
        x = (X[i] - omega[0]) / h[0] + .5 - 1; //subtract 1 for indexing purposes in C
        
        //check if it is a valid x
        if (!boundary && (x<=-2 || x>=m0+1))
            continue;
        
        xf = floor(x);
        x  = x - xf;
        
        if (doDerivative) {
            dT[i] = 0;
        }
        
        for (j=-1;j<3;j++) {
            if (boundary) {
                p = T[min(m0-1,max(0,xf+j))];
            } else {
                p = (xf+j<0 || xf+j>m0-1)? 0: T[xf+j];
            }
            Tc[i] += p*b0(3-j,x-j);
            if (doDerivative) {
                dT[i] += p*db0(3-j,x-j);
            }
        }
        
        if (mxIsNaN(Tc[i]))
            mexPrintf("Is NAN. ");
        
        if (doDerivative) {
            dT[i] = dT[i]/h[0];
        }
    }
}

void splineInter2D(double *Tc, double *dT, const double *T,
        const double *omega, const double *m, int N, const double *X,
        double *h, bool doDerivative, bool boundary) {
    // Increments in X and Y direction
    int xf, yf, i, j, k, m0 = m[0], m1 = m[1];
    double x, y, p;
    
    //for each value X compute the value in  the spline
    #pragma omp parallel for default(shared) private(x, y, xf, yf, p, i, j, k)
    for (i=0;i<N;i++) {
        Tc[i]   = 0;
        
        x = (X[i]   - omega[0]) / h[0] + .5 - 1;
        y = (X[i+N] - omega[2]) / h[1] + .5 - 1;

        //check if it is a valid x
        if (!boundary && (x<=-2 || y<=-2 || x>=m0+1 || y>=m1+1))
            continue;
        
        xf = floor(x);
        yf = floor(y);
        x  = x - xf;
        y  = y - yf;
        
        if (doDerivative) {
            dT[i]   = 0;
            dT[i+N] = 0;
        }

        for (j=-1;j<3;j++) {
            for (k=-1;k<3;k++) {
                if (boundary) {
                    p = T[min(m0-1,max(0,xf+j))+m0*(min(m1-1,max(0,yf+k)))];
                } else {
                    p = (xf+j<0 || xf+j>m0-1 ||
                         yf+k<0 || yf+k>m1-1)? 0: T[xf+j+m0*(yf+k)];
                }
                Tc[i] += p*b0(3-j,x-j)*b0(3-k,y-k);
                if (doDerivative) {
                    dT[i]   += p*db0(3-j,x-j)* b0(3-k,y-k);
                    dT[i+N] += p* b0(3-j,x-j)*db0(3-k,y-k);
                }
            }
        }
        
        if (mxIsNaN(Tc[i]))
            mexPrintf("Is NAN. ");
        
        if (doDerivative) {
            dT[i]   = dT[i]/h[0];
            dT[i+N] = dT[i+N]/h[1];
        }
    }    
}

void splineInter3D(double *Tc, double *dT, const double *T,
        const double *omega, const double *m, int N, const double *X,
        double *h, bool doDerivative, bool boundary) {
    int xf, yf, zf, i, j, k, l;
    int i1 = 1, m0 = m[0], m1 = m[1], m2 = m[2], m01 = m0*m1;
    double x, y, z, p;
    
    // Increments in X,Y and Z direction
    #pragma omp parallel for default(shared) private(x, y, z, xf, yf, zf, p, i, j, k, l)
    for (i=0;i<N;i++) {
        Tc[i] = 0;
        
        x = (X[i]     - omega[0]) / h[0] + .5 - 1;
        y = (X[i+N]   - omega[2]) / h[1] + .5 - 1;
        z = (X[i+2*N] - omega[4]) / h[2] + .5 - 1;
        
        xf = floor(x);
        yf = floor(y);
        zf = floor(z);
        
        //check if it is a valid x
        if (!boundary && (x<=-2 || y<=-2 || z<=-2 || x>=m[0]+1 || y>=m1+1|| z>=m2+1))
            continue;
        
        // Extract remainder
        x = x - xf;
        y = y - yf;
        z = z - zf;
        
        if (doDerivative) {
            dT[i]     = 0;
            dT[i+N]   = 0;
            dT[i+2*N] = 0;
        }
        
        for (j=-1;j<3;j++) {
            for (k=-1;k<3;k++) {
                for (l=-1;l<3;l++) {
                    if (boundary) {
                        p = T[min(m0-1,max(0,xf+j))+m0*(min(m1-1,max(0,yf+k)))+m01*(min(m2-1,max(0,zf+l)))];
                    } else {
                        p = (xf+j<0 || xf+j>m0-1 ||
                             yf+k<0 || yf+k>m1-1 ||
                             zf+l<0 || zf+l>m2-1)? 0: T[xf+j+m0*(yf+k)+m01*(zf+l)];
                    }
                    Tc[i] += p*b0(3-j,x-j)*b0(3-k,y-k)*b0(3-l,z-l);
                    if (doDerivative) {
                        dT[i]     += p*db0(3-j,x-j)* b0(3-k,y-k)* b0(3-l,z-l);
                        dT[i+N]   += p* b0(3-j,x-j)*db0(3-k,y-k)* b0(3-l,z-l);
                        dT[i+2*N] += p* b0(3-j,x-j)* b0(3-k,y-k)*db0(3-l,z-l);
                    }
                }
            }
        }
        
        if (mxIsNaN(Tc[i]))
            mexPrintf("Is NAN. ");
        
        if (doDerivative) {
            dT[i] = dT[i]/h[0]; dT[i+N] = dT[i+N]/h[1]; dT[i+2*N] = dT[i+2*N]/h[2];
        }
    }
}


inline double b0(int j, double xi) {
    switch (j) {
        case 1:
            return (2+xi) * (2+xi) * (2+xi); //  (2+xi).^3;
        case 2:
            return -(3*xi+6)*xi*xi+4;        // -(3*xi+6).*xi.^2+4;
        case 3:
            return  (3*xi-6)*xi*xi+4;        //  (3*xi-6).*xi.^2+4;
        case 4:
            return (2-xi) * (2-xi) * (2-xi); //  (2-xi).^3;
    }
    return 0;
}

inline double db0(int j, double xi) {
    switch (j) {
        case 1:
            return  3*(2+xi)*(2+xi);         //  3*(2+xi).^2;
        case 2:
            return -(9*xi+12)*xi;            // -(9*xi+12).*xi;
        case 3:
            return  (9*xi-12)*xi;            //  (9*xi-12).*xi;
        case 4:  
            return -3*(2-xi)*(2-xi);         // -3*(2-xi).^2;
    }
    return 0;
}

/* ==================================================================================== */
