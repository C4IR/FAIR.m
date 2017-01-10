/*
 * (c) Jan Modersitzki and Fabian Gigengack 2011/04/13, see FAIR.2 and FAIRcopyright.m.
 * http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
 * http://www.uni-muenster.de/EIMI/
 *
 * CPP code for smoothed linear interpolation. See linearInterSmooth for details.
 */

#include <math.h>
#include <mex.h>
#ifdef _OPENMP
    #include <omp.h>
#endif

void linearInterSmooth1D(double *Tc, double *dT, const double *T, const double *omega, const double *m, int N, const double *X, double *h, const double* eta, bool doDerivative);
void linearInterSmooth2D(double *Tc, double *dT, const double *T, const double *omega, const double *m, int N, const double *X, double *h, const double* eta, bool doDerivative);
void linearInterSmooth3D(double *Tc, double *dT, const double *T, const double *omega, const double *m, int N, const double *X, double *h, const double* eta, bool doDerivative);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    int dim;
    if (nrhs<6)
        mexErrMsgTxt("Number of arguments must be 6!");
    
    // Get input
    const double* T     = static_cast<double*>(mxGetData(prhs[0]));
    const double* omega = static_cast<double*>(mxGetData(prhs[1]));
    const double* m     = static_cast<double*>(mxGetData(prhs[2]));
    const double* X     = static_cast<double*>(mxGetData(prhs[3]));
    const double* eta   = static_cast<double*>(mxGetData(prhs[4]));
    bool doDerivative   = (bool) *mxGetLogicals(prhs[5]);
    
    // get dimension, h and number of elements N
    dim = mxGetN(prhs[1]) / 2;
    double h[3]; //h[dim];
    int i;
    for (i=0; i<dim; i++) {
        h[i] = (omega[2*i+1]-omega[2*i]) / m[i]; }
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
            linearInterSmooth1D(Tc,dT,T,omega,m,N,X,h,eta,doDerivative);
            break;
        case 2:
            linearInterSmooth2D(Tc,dT,T,omega,m,N,X,h,eta,doDerivative);
            break;
        case 3:
            linearInterSmooth3D(Tc,dT,T,omega,m,N,X,h,eta,doDerivative);
            break;
        default:
            break;
    }
}

void linearInterSmooth1D(double *Tc,double *dT,const double *T, const double *omega, const double *m,int N, const double *X,double *h,const double* eta,bool doDerivative){
    double x, xi, p[3];
    int xf;
    
    #pragma omp parallel for default(shared) private(x, xf, p, xi)
    for (int i=0;i<N;i++) {
        // map x from [h/2,omega-h/2] -> [0,m-1],
        x = (X[i] - omega[0]) / h[0] + .5 - 1;
        
        // Valid?
        if(x<=-1 || x>=m[0])
            continue;
        
        xf = floor(x);
        // Extract remainder
        xi = x - xf;
        
        if (xi>eta[0] && xi<(1-eta[0])) {
            // Zero padding
            p[0] = (xf  <0 || xf  >m[0]-1)? 0: T[xf];
            p[1] = (xf+1<0 || xf+1>m[0]-1)? 0: T[xf+1];
            
            Tc[i] = p[0] * (1-xi) + p[1] * xi;
            if (doDerivative) {
                dT[i] = (p[1] - p[0]) / h[0];
            }
            continue;
        }
        
        xf = floor(x + .5); // round
        // Extract remainder
        xi = x - xf;
        
        if (abs(xi)<=eta[0]) {
            // Zero padding
            p[0] = (xf-1<0 || xf-1>m[0]-1)? 0: T[xf-1];
            p[1] = (xf  <0 || xf  >m[0]-1)? 0: T[xf];
            p[2] = (xf+1<0 || xf+1>m[0]-1)? 0: T[xf+1];
            
            Tc[i] = p[1] + (p[1] - p[0]) * xi
                    + 1/(2*eta[0]) * (.5 * p[0] - p[1] + .5 * p[2])
                    * (xi+eta[0]) * (xi+eta[0]);
            if (doDerivative) {
                dT[i] = p[1] - p[0] + 1/eta[0] * (.5 * p[0] - p[1]
                        + .5 * p[2]) * (xi+eta[0]);
                dT[i] = dT[i]/h[0];
            }
        }
    }    
}

void linearInterSmooth2D(double *Tc,double *dT,const double *T, const double *omega, const double *m,int N, const double *X,double *h,const double* eta,bool doDerivative){
    double x, y, xif, yif, xir, yir, temp1, temp2, temp3, p[9];
    int j, k, xf, yf, xr, yr;
    // Increments in X and Y direction
    int i1 = 1;
    int i2 = m[0];
    
    #pragma omp parallel for default(shared) private(j, k, x, y, xf, yf, xr, yr, xif, yif, xir, yir, p, temp1, temp2, temp3)
    for (int i=0;i<N;i++) {
        // map x from [h/2,omega-h/2] -> [0,m-1],
        x = (X[i]   - omega[0]) / h[0] + .5 - 1;
        y = (X[i+N] - omega[2]) / h[1] + .5 - 1;
        
        // Valid?
        if (-eta[0]-1>=x || x>=m[0]+eta[0] || -eta[0]-1>=y || y>=m[1]+eta[0])
            continue;
        
        xf = floor(x);
        yf = floor(y);
        xr = floor(x + .5);
        yr = floor(y + .5);
        // Extract remainder
        xif = x - xf;
        yif = y - yf;
        xir = x - xr;
        yir = y - yr;
        
        // Case 1
        // xi1 and xi2 not in eta
        if ((xif>eta[0] && xif<(1-eta[0]) && yif>eta[0] && yif<(1-eta[0]))) {
            // no components in eta --> 2D linear interpolation
            // Zero padding
            for (j=0;j<2;j++) {
                for (k=0;k<2;k++) {
                    p[j+1+(k+1)*3] = (xf+j<0 || xf+j>m[0]-1 || yf+k<0 || yf+k>m[1]-1)? 0: T[xf+j+i2*(yf+k)];
                }
            }
            
            Tc[i] = (p[4] * (1-xif) + p[5] * xif) * (1-yif) + (p[7] * (1-xif) + p[8] * xif) * yif;
            
            if (doDerivative) {
                dT[i]   = ((p[5] - p[4]) * (1-yif) + (p[8] - p[7]) * yif) / h[0];
                dT[i+N] = ((p[7] - p[4]) * (1-xif) + (p[8] - p[5]) * xif) / h[1];
            }
            continue;
        }
        
        // Case 2
        // xi1 in eta and xi2 not in eta
        if ((xif<=eta[0] || xif>=(1-eta[0])) && yif>eta[0] && yif<(1-eta[0])) {
            // Zero padding
            for (j=-1;j<2;j++) {
                for (k=0;k<2;k++) {
                    p[j+1+(k+1)*3] = (xr+j<0 || xr+j>m[0]-1 || yf+k<0 || yf+k>m[1]-1)? 0: T[xr+j+i2*(yf+k)];
                }
            }
            
            temp1 = yif * (p[6] - p[3]) + p[3];
            temp2 = yif * (p[7] - p[4]) + p[4];
            temp3 = yif * (p[8] - p[5]) + p[5];
            
            Tc[i] = temp2 + (temp3 - temp2) * xir
                    + 1 / (4*eta[0]) * (temp3 - 2*temp2 + temp1) * (xir-eta[0]) * (xir-eta[0]);

            if (doDerivative) {
                dT[i]   = ((temp3-temp2) + 1/(2*eta[0])*(temp3-2*temp2+temp1) *(xir-eta[0])) / h[0];
                dT[i+N] = ((p[7]-p[4]) * (1-xir-1/(2*eta[0])*(xir-eta[0])*(xir-eta[0]))
                         + (p[8]-p[5]) * (  xir+1/(4*eta[0])*(xir-eta[0])*(xir-eta[0]))
                         + (p[6]-p[3]) * (      1/(4*eta[0])*(xir-eta[0])*(xir-eta[0]))) / h[1];
            }
            continue;
        }
        
        // Case 3
        // xi1 not in eta, xi2 in eta
        if ((yif<=eta[0] || yif>=(1-eta[0])) && xif>eta[0] && xif<(1-eta[0])) {
            // Zero padding
            for (j=0;j<2;j++) {
                for (k=-1;k<2;k++) {
                    p[j+1+(k+1)*3] = (xf+j<0 || xf+j>m[0]-1 || yr+k<0 || yr+k>m[1]-1)? 0: T[xf+j+i2*(yr+k)];
                }
            }
            
            temp1 = xif * (p[2] - p[1]) + p[1];
            temp2 = xif * (p[5] - p[4]) + p[4];
            temp3 = xif * (p[8] - p[7]) + p[7];
            
            Tc[i] = temp2 + (temp3 - temp2) * yir + 1 / (4*eta[0]) * (temp3 - 2*temp2 + temp1) * (yir-eta[0]) * (yir-eta[0]);
            
            if (doDerivative) {
                dT[i]   = ((p[2]-p[1]) * (      1/(4*eta[0])*(yir-eta[0])*(yir-eta[0]))
                         + (p[5]-p[4]) * (1-yir-1/(2*eta[0])*(yir-eta[0])*(yir-eta[0]))
                         + (p[8]-p[7]) * (  yir+1/(4*eta[0])*(yir-eta[0])*(yir-eta[0])))/h[0];
                dT[i+N] = ((temp3-temp2) + 1/(2*eta[0])*(temp3-2*temp2+temp1) * (yir-eta[0]))/h[1];
            }
            continue;
        }
        
        // Case 4
        // xi1 and xi2 in eta
        if ((xif<=eta[0] && yif<=eta[0]) || (xif>=(1-eta[0]) && yif<=eta[0])
            || (xif<=eta[0] && yif>=(1-eta[0]))
            || (xif>=(1-eta[0]) && yif>=(1-eta[0])) ) {
            // Zero padding
            for (j=-1;j<2;j++) {
                for (k=-1;k<2;k++) {
                    p[j+1+(k+1)*3] = (xr+j<0 || xr+j>m[0]-1 || yr+k<0 || yr+k>m[1]-1)? 0: T[xr+j+i2*(yr+k)];
                }
            }
            
            temp1 = p[3] + (p[6] - p[3]) * yir + 1/(4*eta[0])
                 * (p[6] - 2 * p[3] + p[0]) * (yir-eta[0]) * (yir-eta[0]);
            temp2 = p[4] + (p[7] - p[4]) * yir + 1/(4*eta[0])
                 * (p[7] - 2 * p[4] + p[1]) * (yir-eta[0]) * (yir-eta[0]);
            temp3 = p[5] + (p[8] - p[5]) * yir + 1/(4*eta[0])
                 * (p[8] - 2 * p[5] + p[2]) * (yir-eta[0]) * (yir-eta[0]);
            
            Tc[i] = temp2+(temp3-temp2) * xir + 1/(4*eta[0])
                 * (temp3-2*temp2+temp1) * (xir-eta[0]) * (xir-eta[0]);
            
            if (doDerivative) {
                dT[i] = ((temp3-temp2) + 1/(2*eta[0]) * (temp3-2*temp2+temp1) * (xir-eta[0]))/h[0];
                dT[i+N] = ((p[6]-p[3]+1/(2*eta[0])*(p[6]-2*p[3]+p[0])*(yir-eta[0]))
                    *(1/(4*eta[0])*(xir-eta[0])*(xir-eta[0]))
                    +(p[7]-p[4]+1/(2*eta[0])*(p[7]-2*p[4]+p[1])*(yir-eta[0]))
                    *(1-xir-1/(2*eta[0])*(xir-eta[0])*(xir-eta[0]))
                    +(p[8]-p[5]+1/(2*eta[0])*(p[8]-2*p[5]+p[2])*(yir-eta[0]))
                    *(xir+1/(4*eta[0])*(xir-eta[0])*(xir-eta[0])))/h[1];
            }
        }
    }
}

void linearInterSmooth3D(double *Tc,double *dT,const double *T, const double *omega, const double *m,int N, const double *X,double *h,const double* eta,bool doDerivative){
    double x, y, z, xif, yif, zif, xir, yir, zir, temp1, temp2, temp3;
    double t[9], p[27];
    int j, k, l, xf, yf, zf, xr, yr, zr;
    // Increments in X and Y direction
    int i1 = 1, i2 = m[0], i3 = i2*m[1];
    
    #pragma omp parallel for default(shared) private(j, k, l, x, y, z, xf, yf, zf, xr, yr, zr, xif, yif, zif, xir, yir, zir, temp1, temp2, temp3, t, p)
    for (int i=0;i<N;i++) {
        // map x from [h/2,omega-h/2] -> [0,m-1],
        x = (X[i]     - omega[0]) / h[0] + .5 - 1;
        y = (X[i+N]   - omega[2]) / h[1] + .5 - 1;
        z = (X[i+2*N] - omega[4]) / h[2] + .5 - 1;
        
        // Valid?
        if (-eta[0]-1>=x || x>=m[0]+eta[0] || -eta[0]-1>=y || y>=m[1]+eta[0] || -eta[0]-1>=z || z>=m[1]+eta[0])
            continue;
        
        xf = floor(x);
        yf = floor(y);
        zf = floor(z);
        xr = floor(x + .5);
        yr = floor(y + .5);
        zr = floor(z + .5);
        
        // Extract remainder
        xif = x - xf;
        yif = y - yf;
        zif = z - zf;
        xir = x - xr;
        yir = y - yr;
        zir = z - zr;
        
        // first case: no component in eta[0]
        if (xif>eta[0] && xif<(1-eta[0]) && yif>eta[0]
            && yif<(1-eta[0]) && zif>eta[0] && zif<(1-eta[0])) {
            // no component in eta[0] --> 3D linear interpolation
            // Zero padding
            for (j=-1;j<2;j++) {
                for (k=-1;k<2;k++) {
                    for (l=-1;l<2;l++) {
                        p[j+1+(k+1)*3+(l+1)*9] = (xf+j<0 || xf+j>m[0]-1 || yf+k<0 || yf+k>m[1]-1 || zf+l<0 || zf+l>m[2]-1)? 0: T[xf+j+i2*(yf+k)+i3*(zf+l)];
                    }
                }
            }
            
            Tc[i] = ((p[13]*(1-xif)+p[14]*xif)*(1-yif)
                  +  (p[16]*(1-xif)+p[17]*xif)*   yif)*(1-zif)
                  + ((p[22]*(1-xif)+p[23]*xif)*(1-yif)
                  +  (p[25]*(1-xif)+p[26]*xif)*   yif)*(zif);
            
            // compute derivative, if necessary
            if (doDerivative) {
                dT[i]     = ((p[14]-p[13])*(1-yif)+(p[17]-p[16])*yif)*(1-zif)
                           +((p[23]-p[22])*(1-yif)+(p[26]-p[25])*yif)*(zif);
                dT[i+N]   = ((p[16]-p[13])*(1-xif)+(p[17]-p[14])*xif)*(1-zif)
                           +((p[25]-p[22])*(1-xif)+(p[26]-p[23])*xif)*(zif);
                dT[i+2*N] = ((p[22]*(1-xif)+p[23]*xif)*(1-yif)+(p[25]*(1-xif)+p[26]*xif)*(yif))
                           -((p[13]*(1-xif)+p[14]*xif)*(1-yif)+(p[16]*(1-xif)+p[17]*xif)*(yif));
                
                dT[i] = dT[i]/h[0]; dT[i+N] = dT[i+N]/h[1]; dT[i+2*N] = dT[i+2*N]/h[2];
            }
            continue;
        }
        
        // --- quadratic part ------------------------------------------------------
        //
        //
        // second case: exactly one direction in eta[0]
        // there are three cases
        // case A: x3 in eta[0]  -  x1,x2 not in eta[0]
        // case B: x1 in eta[0]  -  x2,x3 not in eta[0]
        // case C: x2 in eta[0]  -  x1,x3 not in eta[0]
        // idea: compute 3 points via 2D linear interpolation on parallel planes
        // then quadratic interpolation in the third direction (component in eta[0])
        
        // case A: x3 in eta[0]  -  x1,x2 not
        if (xif>eta[0] && xif<(1-eta[0]) && yif>eta[0] && yif<(1-eta[0])
            && (zif<=eta[0] || zif>=(1-eta[0]))) {
            // Zero padding
            for (j=-1;j<2;j++) {
                for (k=-1;k<2;k++) {
                    for (l=-1;l<2;l++) {
                        p[j+1+(k+1)*3+(l+1)*9] = (xf+j<0 || xf+j>m[0]-1 || yf+k<0 || yf+k>m[1]-1 || zr+l<0 || zr+l>m[2]-1)? 0: T[xf+j+i2*(yf+k)+i3*(zr+l)];
                    }
                }
            }
            
            // linear interpolation on three parallel planes (x1,x2)
            temp1 = (p[4]  * (1-xif) + p[5]  * xif) * (1-yif)
                  + (p[7]  * (1-xif) + p[8]  * xif) * yif;
            temp2 = (p[13] * (1-xif) + p[14] * xif) * (1-yif)
                  + (p[16] * (1-xif) + p[17] * xif) * yif;
            temp3 = (p[22] * (1-xif) + p[23] * xif) * (1-yif)
                  + (p[25] * (1-xif) + p[26] * xif) * yif;
            
            // quadratic interpolation in x3-direction
            Tc[i] = temp2+(temp3-temp2)*zir
            +1/(4*eta[0])*(temp3-2*temp2+temp1)*(zir-eta[0])*(zir-eta[0]);
            
            // compute derivative, if necessary
            if (doDerivative) {
                dT[i+2*N] = (temp3-temp2)+1/(2*eta[0])*(temp3-2*temp2+temp1)*(zir-eta[0]);
                // dtemp/dx1
                temp1 = (p[5] -p[4]) *(1-yif)+(p[8] -p[7]) *yif;
                temp2 = (p[14]-p[13])*(1-yif)+(p[17]-p[16])*yif;
                temp3 = (p[23]-p[22])*(1-yif)+(p[26]-p[25])*yif;
                dT[i] = temp2+(temp3-temp2)*zir
                +1/(4*eta[0])*(temp3-2*temp2+temp1)*(zir-eta[0])*(zir-eta[0]);
                // dtemp/dx2
                temp1 = (p[7] -p[4]) *(1-xif)+(p[8] -p[5]) *xif;
                temp2 = (p[16]-p[13])*(1-xif)+(p[17]-p[14])*xif;
                temp3 = (p[25]-p[22])*(1-xif)+(p[26]-p[23])*xif;
                dT[i+N] = temp2+(temp3-temp2)*zir
                +1/(4*eta[0])*(temp3-2*temp2+temp1)*(zir-eta[0])*(zir-eta[0]);
                
                dT[i] = dT[i]/h[0]; dT[i+N] = dT[i+N]/h[1]; dT[i+2*N] = dT[i+2*N]/h[2];
            }
            continue;
        }
        
        // case B: x1 in eta[0]  -  x2, x3 not in eta[0]
        if (yif>eta[0] && yif<(1-eta[0]) && zif>eta[0] && zif<(1-eta[0])
            && (xif<=eta[0] || xif>=(1-eta[0]))) {
            // Zero padding
            for (j=-1;j<2;j++) {
                for (k=-1;k<2;k++) {
                    for (l=-1;l<2;l++) {
                        p[j+1+(k+1)*3+(l+1)*9] = (xr+j<0 || xr+j>m[0]-1 || yf+k<0 || yf+k>m[1]-1 || zf+l<0 || zf+l>m[2]-1)? 0: T[xr+j+i2*(yf+k)+i3*(zf+l)];
                    }
                }
            }
            
            // linear interpolation on three parallel planes (x2,x3)
            temp1 = (p[12] * (1-yif) + p[15] * yif) * (1-zif)
                  + (p[21] * (1-yif) + p[24] * yif) * zif;
            temp2 = (p[13] * (1-yif) + p[16] * yif) * (1-zif)
                  + (p[22] * (1-yif) + p[25] * yif) * zif;
            temp3 = (p[14] * (1-yif) + p[17] * yif) * (1-zif)
                  + (p[23] * (1-yif) + p[26] * yif) * zif;
            
            // quadratic interpolation in x1-direction
            Tc[i] = temp2+(temp3-temp2)*xir
            +1/(4*eta[0])*(temp3-2*temp2+temp1)*(xir-eta[0])*(xir-eta[0]);
            // compute derivative, if necessary
            if (doDerivative) {
                dT[i] = (temp3-temp2)+1/(2*eta[0])*(temp3-2*temp2+temp1)*(xir-eta[0]);
                // dtemp/dx2
                temp1 = (p[15]-p[12])*(1-zif)+(p[24]-p[21])*zif;
                temp2 = (p[16]-p[13])*(1-zif)+(p[25]-p[22])*zif;
                temp3 = (p[17]-p[14])*(1-zif)+(p[26]-p[23])*zif;
                dT[i+N] = temp2+(temp3-temp2)*xir
                +1/(4*eta[0])*(temp3-2*temp2+temp1)*(xir-eta[0])*(xir-eta[0]);
                // dtemp/dx3
                temp1 = (p[21]-p[12])*(1-yif)+(p[24]-p[15])*yif;
                temp2 = (p[22]-p[13])*(1-yif)+(p[25]-p[16])*yif;
                temp3 = (p[23]-p[14])*(1-yif)+(p[26]-p[17])*yif;
                dT[i+2*N] = temp2+(temp3-temp2)*xir
                +1/(4*eta[0])*(temp3-2*temp2+temp1)*(xir-eta[0])*(xir-eta[0]);
                
                dT[i] = dT[i]/h[0]; dT[i+N] = dT[i+N]/h[1]; dT[i+2*N] = dT[i+2*N]/h[2];
            }
            continue;
        }
        
        // case C: x2 in eta[0]  - x1, x3 not in eta[0]
        if (xif>eta[0] && xif<(1-eta[0]) && zif>eta[0] && zif<(1-eta[0])
            && (yif<=eta[0] || yif>=(1-eta[0]))) {
            // Zero padding
            for (j=-1;j<2;j++) {
                for (k=-1;k<2;k++) {
                    for (l=-1;l<2;l++) {
                        p[j+1+(k+1)*3+(l+1)*9] = (xf+j<0 || xf+j>m[0]-1 || yr+k<0 || yr+k>m[1]-1 || zf+l<0 || zf+l>m[2]-1)? 0: T[xf+j+i2*(yr+k)+i3*(zf+l)];
                    }
                }
            }
            
            // linear interpolation on three parallel planes (x1,x3)
            temp1 = (p[10]*(1-xif)+p[11]*xif)*(1-zif)
                   +(p[19]*(1-xif)+p[20]*xif)*zif;
            temp2 = (p[13]*(1-xif)+p[14]*xif)*(1-zif)
                   +(p[22]*(1-xif)+p[23]*xif)*zif;
            temp3 = (p[16]*(1-xif)+p[17]*xif)*(1-zif)
                   +(p[25]*(1-xif)+p[26]*xif)*zif;
            // quadratic interpolation in x2-direction
            Tc[i] = temp2+(temp3-temp2)*yir
            +1/(4*eta[0])*(temp3-2*temp2+temp1)*(yir-eta[0])*(yir-eta[0]);
            
            if (doDerivative) {
                dT[i+N] = temp3-temp2+1/(2*eta[0])*(temp3-2*temp2+temp1)*(yir-eta[0]);
                // dtemp/dx1
                temp1 = (p[11]-p[10])*(1-zif)+(p[20]-p[19])*zif;
                temp2 = (p[14]-p[13])*(1-zif)+(p[23]-p[22])*zif;
                temp3 = (p[17]-p[16])*(1-zif)+(p[26]-p[25])*zif;
                dT[i] = temp2+(temp3-temp2)*yir
                +1/(4*eta[0])*(temp3-2*temp2+temp1)*(yir-eta[0])*(yir-eta[0]);
                // dtemp/dx3
                temp1 = (p[19]-p[10])*(1-xif)+(p[20]-p[11])*xif;
                temp2 = (p[22]-p[13])*(1-xif)+(p[23]-p[14])*xif;
                temp3 = (p[25]-p[16])*(1-xif)+(p[26]-p[17])*xif;
                dT[i+2*N] = temp2+(temp3-temp2)*yir
                +1/(4*eta[0])*(temp3-2*temp2+temp1)*(yir-eta[0])*(yir-eta[0]);
                
                dT[i] = dT[i]/h[0]; dT[i+N] = dT[i+N]/h[1]; dT[i+2*N] = dT[i+2*N]/h[2];
            }
            continue;
        }
        
        // third case: exactly two compentents in eta[0]
        // again we have to consider 3 cases
        // case A: x1,x3 in eta[0] - x2 not in eta[0]
        // case B: x1,x2 in eta[0] - x3 not in eta[0]
        // case C: x2,x3 in eta[0] - x1 not in eta[0]
        // idea: compute 3 points on parallel planes, the planes are given
        //       by the condition - one component is in eta[0], the other is not in eta[0]
        //       This points can be handled as in 2D.
        //       At last a quadratic interpolation can be performed
        
        // case A: x1, x3 in eta[0], x2 not in eta[0]
        if ((xif<=eta[0] || xif>=(1-eta[0]))&& yif>eta[0] && yif<(1-eta[0])
            && (zif<=eta[0] || zif>=(1-eta[0]))) {
            // Zero padding
            for (j=-1;j<2;j++) {
                for (k=-1;k<2;k++) {
                    for (l=-1;l<2;l++) {
                        p[j+1+(k+1)*3+(l+1)*9] = (xr+j<0 || xr+j>m[0]-1 || yf+k<0 || yf+k>m[1]-1 || zr+l<0 || zr+l>m[2]-1)? 0: T[xr+j+i2*(yf+k)+i3*(zr+l)];
                    }
                }
            }
            
            // compute 3 points on 3 parallel planes (x1,x2)
            // for each of the 3 points on the planes, 3 more points are needed
            // (computed by linear interpolation in x2-direction)
            t[0] = yif *(p[6] -p[3]) +p[3];
            t[1] = yif *(p[7] -p[4]) +p[4];
            t[2] = yif *(p[8] -p[5]) +p[5];
            t[3] = yif *(p[15]-p[12])+p[12];
            t[4] = yif *(p[16]-p[13])+p[13];
            t[5] = yif *(p[17]-p[14])+p[14];
            t[6] = yif *(p[24]-p[21])+p[21];
            t[7] = yif *(p[25]-p[22])+p[22];
            t[8] = yif *(p[26]-p[23])+p[23];
            // compute points on planes via quadratic interpolation in x1-direction
            temp1 = t[1]+(t[2]-t[1])*xir
            +1/(4*eta[0])*(t[2]-2*t[1]+t[0])*(xir-eta[0])*(xir-eta[0]);
            temp2 = t[4]+(t[5]-t[4])*xir
            +1/(4*eta[0])*(t[5]-2*t[4]+t[3])*(xir-eta[0])*(xir-eta[0]);
            temp3 = t[7]+(t[8]-t[7])*xir
            +1/(4*eta[0])*(t[8]-2*t[7]+t[6])*(xir-eta[0])*(xir-eta[0]);
            // quadratic interpolation in x3-direction
            Tc[i] = temp2+(temp3-temp2)*zir
            +1/(4*eta[0])*(temp3-2*temp2+temp1)*(zir-eta[0])*(zir-eta[0]);
            
            if (doDerivative) {
                dT[i+2*N] = temp3-temp2+1/(2*eta[0])*(temp3-2*temp2+temp1)*(zir-eta[0]);
                // dtemp/dx1
                temp1 = t[2]-t[1]+1/(2*eta[0])*(t[2]-2*t[1]+t[0])*(xir-eta[0]);
                temp2 = t[5]-t[4]+1/(2*eta[0])*(t[5]-2*t[4]+t[3])*(xir-eta[0]);
                temp3 = t[8]-t[7]+1/(2*eta[0])*(t[8]-2*t[7]+t[6])*(xir-eta[0]);
                dT[i] = temp2+(temp3-temp2)*zir
                +1/(4*eta[0])*(temp3-2*temp2+temp1)*(zir-eta[0])*(zir-eta[0]);
                // dt/dx2
                t[0] = p[6] -p[3];
                t[1] = p[7] -p[4];
                t[2] = p[8] -p[5];
                t[3] = p[15]-p[12];
                t[4] = p[16]-p[13];
                t[5] = p[17]-p[14];
                t[6] = p[24]-p[21];
                t[7] = p[25]-p[22];
                t[8] = p[26]-p[23];
                // dtemp/dx2
                temp1 = t[1]+(t[2]-t[1])*xir+1/(4*eta[0])*(t[2]-2*t[1]+t[0])*(xir-eta[0])*(xir-eta[0]);
                temp2 = t[4]+(t[5]-t[4])*xir+1/(4*eta[0])*(t[5]-2*t[4]+t[3])*(xir-eta[0])*(xir-eta[0]);
                temp3 = t[7]+(t[8]-t[7])*xir+1/(4*eta[0])*(t[8]-2*t[7]+t[6])*(xir-eta[0])*(xir-eta[0]);
                dT[i+N] = temp2+(temp3-temp2)*zir
                +1/(4*eta[0])*(temp3-2*temp2+temp1)*(zir-eta[0])*(zir-eta[0]);
                
                dT[i] = dT[i]/h[0]; dT[i+N] = dT[i+N]/h[1]; dT[i+2*N] = dT[i+2*N]/h[2];
            }
            continue;
        }
        
        // case B: x1, x2 in eta[0], x3 not in eta[0]
        if ((xif<=eta[0] || xif >=(1-eta[0])) && (yif<=eta[0] || yif >=(1-eta[0]))
            && zif>eta[0] && zif<(1-eta[0])) {
            // Zero padding
            for (j=-1;j<2;j++) {
                for (k=-1;k<2;k++) {
                    for (l=-1;l<2;l++) {
                        p[j+1+(k+1)*3+(l+1)*9] = (xr+j<0 || xr+j>m[0]-1 || yr+k<0 || yr+k>m[1]-1 || zf+l<0 || zf+l>m[2]-1)? 0: T[xr+j+i2*(yr+k)+i3*(zf+l)];
                    }
                }
            }
            
            // compute 3 points on 3 parallel planes (x1,x3)
            // for each of the 3 points on the planes, 3 more points are needed
            // (computed by linear interpolation in x3-direction)
            t[0] = zif*(p[18]-p[9]) +p[9];
            t[1] = zif*(p[19]-p[10])+p[10];
            t[2] = zif*(p[20]-p[11])+p[11];
            t[3] = zif*(p[21]-p[12])+p[12];
            t[4] = zif*(p[22]-p[13])+p[13];
            t[5] = zif*(p[23]-p[14])+p[14];
            t[6] = zif*(p[24]-p[15])+p[15];
            t[7] = zif*(p[25]-p[16])+p[16];
            t[8] = zif*(p[26]-p[17])+p[17];
            // compute points on planes via quadratic interpolation in x1-direction
            temp1 = t[1]+(t[2]-t[1])*xir
            +1/(4*eta[0])*(t[2]-2*t[1]+t[0])*(xir-eta[0])*(xir-eta[0]);
            temp2 = t[4]+(t[5]-t[4])*xir
            +1/(4*eta[0])*(t[5]-2*t[4]+t[3])*(xir-eta[0])*(xir-eta[0]);
            temp3 = t[7]+(t[8]-t[7])*xir
            +1/(4*eta[0])*(t[8]-2*t[7]+t[6])*(xir-eta[0])*(xir-eta[0]);
            // quadratic interpolation in x2-direction
            Tc[i] = temp2+(temp3-temp2)*yir
            +1/(4*eta[0])*(temp3-2*temp2+temp1)*(yir-eta[0])*(yir-eta[0]);
            // compute derivative, if necessary
            if (doDerivative) {
                dT[i+N] = temp3-temp2+1/(2*eta[0])*(temp3-2*temp2+temp1)*(yir-eta[0]);
                // dtemp/dx1
                temp1 = t[2]-t[1]+1/(2*eta[0])*(t[2]-2*t[1]+t[0])*(xir-eta[0]);
                temp2 = t[5]-t[4]+1/(2*eta[0])*(t[5]-2*t[4]+t[3])*(xir-eta[0]);
                temp3 = t[8]-t[7]+1/(2*eta[0])*(t[8]-2*t[7]+t[6])*(xir-eta[0]);
                dT[i] = temp2+(temp3-temp2)*yir
                +1/(4*eta[0])*(temp3-2*temp2+temp1)*(yir-eta[0])*(yir-eta[0]);
                //dt/dx3
                t[0] = p[18]-p[9];
                t[1] = p[19]-p[10];
                t[2] = p[20]-p[11];
                t[3] = p[21]-p[12];
                t[4] = p[22]-p[13];
                t[5] = p[23]-p[14];
                t[6] = p[24]-p[15];
                t[7] = p[25]-p[16];
                t[8] = p[26]-p[17];
                // dtemp/dx3
                temp1 = t[1]+(t[2]-t[1])*xir+1/(4*eta[0])*(t[2]-2*t[1]+t[0])*(xir-eta[0])*(xir-eta[0]);
                temp2 = t[4]+(t[5]-t[4])*xir+1/(4*eta[0])*(t[5]-2*t[4]+t[3])*(xir-eta[0])*(xir-eta[0]);
                temp3 = t[7]+(t[8]-t[7])*xir+1/(4*eta[0])*(t[8]-2*t[7]+t[6])*(xir-eta[0])*(xir-eta[0]);
                dT[i+2*N] = temp2+(temp3-temp2)*yir
                +1/(4*eta[0])*(temp3-2*temp2+temp1)*(yir-eta[0])*(yir-eta[0]);
                
                dT[i] = dT[i]/h[0]; dT[i+N] = dT[i+N]/h[1]; dT[i+2*N] = dT[i+2*N]/h[2];
            }
            continue;
        }
        
        // case C: x2, x3 in eta[0], x1 not in eta[0]
        if (xif>eta[0] && xif<(1-eta[0]) && (yif<=eta[0] || yif >=(1-eta[0]))
            && (zif<=eta[0] || zif >=(1-eta[0]))) {
            // Zero padding
            for (j=-1;j<2;j++) {
                for (k=-1;k<2;k++) {
                    for (l=-1;l<2;l++) {
                        p[j+1+(k+1)*3+(l+1)*9] = (xf+j<0 || xf+j>m[0]-1 || yr+k<0 || yr+k>m[1]-1 || zr+l<0 || zr+l>m[2]-1)? 0: T[xf+j+i2*(yr+k)+i3*(zr+l)];
                    }
                }
            }
            
            // compute 3 points on 3 parallel planes (x1,x3)
            // for each of the 3 points on the planes, 3 more points are needed
            // (computed by linear interpolation in x1-direction)
            t[0] = xif*(p[2] -p[1]) +p[1];
            t[1] = xif*(p[11]-p[10])+p[10];
            t[2] = xif*(p[20]-p[19])+p[19];
            t[3] = xif*(p[5] -p[4]) +p[4];
            t[4] = xif*(p[14]-p[13])+p[13];
            t[5] = xif*(p[23]-p[22])+p[22];
            t[6] = xif*(p[8] -p[7]) +p[7];
            t[7] = xif*(p[17]-p[16])+p[16];
            t[8] = xif*(p[26]-p[25])+p[25];
            
            // compute points on planes via quadratic interpolation in x1-direction
            temp1 = t[1]+(t[2]-t[1])*zir+1/(4*eta[0])*(t[2]-2*t[1]+t[0])*(zir-eta[0])*(zir-eta[0]);
            temp2 = t[4]+(t[5]-t[4])*zir+1/(4*eta[0])*(t[5]-2*t[4]+t[3])*(zir-eta[0])*(zir-eta[0]);
            temp3 = t[7]+(t[8]-t[7])*zir+1/(4*eta[0])*(t[8]-2*t[7]+t[6])*(zir-eta[0])*(zir-eta[0]);
            
            // quadratic interpolation in x2-direction
            Tc[i] = temp2+(temp3-temp2)*yir
            +1/(4*eta[0])*(temp3-2*temp2+temp1)*(yir-eta[0])*(yir-eta[0]);
            
            // compute derivative, if necessary
            if (doDerivative) {
                dT[i+N] = temp3-temp2+1/(2*eta[0])*(temp3-2*temp2+temp1)*(yir-eta[0]);
                // dtemp/dx3
                temp1 = t[2]-t[1]+1/(2*eta[0])*(t[2]-2*t[1]+t[0])*(zir-eta[0]);
                temp2 = t[5]-t[4]+1/(2*eta[0])*(t[5]-2*t[4]+t[3])*(zir-eta[0]);
                temp3 = t[8]-t[7]+1/(2*eta[0])*(t[8]-2*t[7]+t[6])*(zir-eta[0]);
                dT[i+2*N] = temp2+(temp3-temp2)*yir
                +1/(4*eta[0])*(temp3-2*temp2+temp1)*(yir-eta[0])*(yir-eta[0]);
                // dt/dx1
                t[0] = p[2] -p[1];
                t[1] = p[11]-p[10];
                t[2] = p[20]-p[19];
                t[3] = p[5] -p[4];
                t[4] = p[14]-p[13];
                t[5] = p[23]-p[22];
                t[6] = p[8] -p[7];
                t[7] = p[17]-p[16];
                t[8] = p[26]-p[25];
                // dtemp/dx1
                temp1 = t[1]+(t[2]-t[1])*zir+1/(4*eta[0])*(t[2]-2*t[1]+t[0])*(zir-eta[0])*(zir-eta[0]);
                temp2 = t[4]+(t[5]-t[4])*zir+1/(4*eta[0])*(t[5]-2*t[4]+t[3])*(zir-eta[0])*(zir-eta[0]);
                temp3 = t[7]+(t[8]-t[7])*zir+1/(4*eta[0])*(t[8]-2*t[7]+t[6])*(zir-eta[0])*(zir-eta[0]);
                dT[i] = temp2+(temp3-temp2)*yir
                +1/(4*eta[0])*(temp3-2*temp2+temp1)*(yir-eta[0])*(yir-eta[0]);
                
                dT[i] = dT[i]/h[0]; dT[i+N] = dT[i+N]/h[1]; dT[i+2*N] = dT[i+2*N]/h[2];
            }
            continue;
        }
        
        // third case: all components in eta[0]
        // idea: compute 3 points on parallel planes the same way as in 2D
        //       (case: both components in eta[0]). Then again quadratic interpolation
        //       in the remaining direction
        
        // all components in eta[0]
        if ((xif<=eta[0] || xif>=(1-eta[0])) && (yif<=eta[0] || yif >=(1-eta[0]))
            && (zif<=eta[0] || zif >=(1-eta[0]))) {
            // Zero padding
            for (j=-1;j<2;j++) {
                for (k=-1;k<2;k++) {
                    for (l=-1;l<2;l++) {
                        p[j+1+(k+1)*3+(l+1)*9] = (xr+j<0 || xr+j>m[0]-1 || yr+k<0 || yr+k>m[1]-1 || zr+l<0 || zr+l>m[2]-1)? 0: T[xr+j+i2*(yr+k)+i3*(zr+l)];
                    }
                }
            }
            
            // quadratic interpolation in x2-direction
            t[0] = p[3]+(p[6]-p[3])*yir+1/(4*eta[0])*(p[6]-2*p[3]+p[0])
                *(yir-eta[0])*(yir-eta[0]);
            t[1] = p[4]+(p[7]-p[4])*yir+1/(4*eta[0])*(p[7]-2*p[4]+p[1])
                *(yir-eta[0])*(yir-eta[0]);
            t[2] = p[5]+(p[8]-p[5])*yir+1/(4*eta[0])*(p[8]-2*p[5]+p[2])
                *(yir-eta[0])*(yir-eta[0]);
            t[3] = p[12]+(p[15]-p[12])*yir+1/(4*eta[0])*(p[15]-2*p[12]+p[9])
                *(yir-eta[0])*(yir-eta[0]);
            t[4] = p[13]+(p[16]-p[13])*yir+1/(4*eta[0])*(p[16]-2*p[13]+p[10])
                *(yir-eta[0])*(yir-eta[0]);
            t[5] = p[14]+(p[17]-p[14])*yir+1/(4*eta[0])*(p[17]-2*p[14]+p[11])
                *(yir-eta[0])*(yir-eta[0]);
            t[6] = p[21]+(p[24]-p[21])*yir+1/(4*eta[0])*(p[24]-2*p[21]+p[18])
                *(yir-eta[0])*(yir-eta[0]);
            t[7] = p[22]+(p[25]-p[22])*yir+1/(4*eta[0])*(p[25]-2*p[22]+p[19])
                *(yir-eta[0])*(yir-eta[0]);
            t[8] = p[23]+(p[26]-p[23])*yir+1/(4*eta[0])*(p[26]-2*p[23]+p[20])
                *(yir-eta[0])*(yir-eta[0]);
            
            // quadratic interpolation in x1-direction
            temp1 = t[1]+(t[2]-t[1])*xir
            +1/(4*eta[0])*(t[2]-2*t[1]+t[0])*(xir-eta[0])*(xir-eta[0]);
            temp2 = t[4]+(t[5]-t[4])*xir
            +1/(4*eta[0])*(t[5]-2*t[4]+t[3])*(xir-eta[0])*(xir-eta[0]);
            temp3 = t[7]+(t[8]-t[7])*xir
            +1/(4*eta[0])*(t[8]-2*t[7]+t[6])*(xir-eta[0])*(xir-eta[0]);
            
            // quadratic interpolation in x3-direction
            Tc[i] = temp2+(temp3-temp2)*zir
            +1/(4*eta[0])*(temp3-2*temp2+temp1)*(zir-eta[0])*(zir-eta[0]);
            // compute derivative, if necessary
            if (doDerivative) {
                dT[i+2*N] = temp3-temp2+1/(2*eta[0])*(temp3-2*temp2+temp1)*(zir-eta[0]);
                
                // dtemp/dx1
                temp1 = t[2]-t[1]+1/(2*eta[0])*(t[2]-2*t[1]+t[0])*(xir-eta[0]);
                temp2 = t[5]-t[4]+1/(2*eta[0])*(t[5]-2*t[4]+t[3])*(xir-eta[0]);
                temp3 = t[8]-t[7]+1/(2*eta[0])*(t[8]-2*t[7]+t[6])*(xir-eta[0]);
                
                dT[i] = temp2+(temp3-temp2)*zir
                +1/(4*eta[0])*(temp3-2*temp2+temp1)*(zir-eta[0])*(zir-eta[0]);
                // dt/dx2
                t[0] = p[6] -p[3] +1/(2*eta[0])*(p[6] -2*p[3] +p[0]) *(yir-eta[0]);
                t[1] = p[7] -p[4] +1/(2*eta[0])*(p[7] -2*p[4] +p[1]) *(yir-eta[0]);
                t[2] = p[8] -p[5] +1/(2*eta[0])*(p[8] -2*p[5] +p[2]) *(yir-eta[0]);
                t[3] = p[15]-p[12]+1/(2*eta[0])*(p[15]-2*p[12]+p[9]) *(yir-eta[0]);
                t[4] = p[16]-p[13]+1/(2*eta[0])*(p[16]-2*p[13]+p[10])*(yir-eta[0]);
                t[5] = p[17]-p[14]+1/(2*eta[0])*(p[17]-2*p[14]+p[11])*(yir-eta[0]);
                t[6] = p[24]-p[21]+1/(2*eta[0])*(p[24]-2*p[21]+p[18])*(yir-eta[0]);
                t[7] = p[25]-p[22]+1/(2*eta[0])*(p[25]-2*p[22]+p[19])*(yir-eta[0]);
                t[8] = p[26]-p[23]+1/(2*eta[0])*(p[26]-2*p[23]+p[20])*(yir-eta[0]);
                // dtemp/dx2
                temp1 = t[1]+(t[2]-t[1])*xir+1/(4*eta[0])*(t[2]-2*t[1]+t[0])*(xir-eta[0])*(xir-eta[0]);
                temp2 = t[4]+(t[5]-t[4])*xir+1/(4*eta[0])*(t[5]-2*t[4]+t[3])*(xir-eta[0])*(xir-eta[0]);
                temp3 = t[7]+(t[8]-t[7])*xir+1/(4*eta[0])*(t[8]-2*t[7]+t[6])*(xir-eta[0])*(xir-eta[0]);
                
                dT[i+N] = temp2+(temp3-temp2)*zir
                +1/(4*eta[0])*(temp3-2*temp2+temp1)*(zir-eta[0])*(zir-eta[0]);
                
                dT[i] = dT[i]/h[0]; dT[i+N] = dT[i+N]/h[1]; dT[i+2*N] = dT[i+2*N]/h[2];
            }
            continue;
        }
    }
}

/*==================================================================================== */
