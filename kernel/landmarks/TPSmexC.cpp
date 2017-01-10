//authors: Thomas Polzin, Jan Ruehaak
//date: 15.08.2013

#include <math.h>
#include <mex.h>
#include <algorithm>
#include <stdio.h>
#include <iostream>
//#include "TimingUtils.h"

enum eRegistrationCodes {
    NoError = 0,
    UnknownError,
    OutOfMemoryError,
    InvalidInputImageError,
    InvalidInputImageWarning,
    InvalidInputTransformationError,
    InvalidOutputTransformationError,
    InvalidParameterError,
    NumericalError,
    InternalError,
    NullPointerError,
    FinalLevelReached
};


using namespace std;

void ComputeCoefficients(double* c, const double* t, double* r, const int dim, const int L, const double theta);
void EvaluateTPS(double* yTPS, double* yr, double* c, double* xc, double* r,  const int dim, const int L, const int N);
//Evaluate TPS on cell-centered grid --> compute x from omega and m
void EvaluateTPScc(double* yTPS, double* yr, double* c, double* omega, double* m, double* r,  const int dim, const int L);
int solve(double* LS, const double* rhs, double* x, const int m);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    //Initialize output
    double *yTPS, *yr, *c; //transformed grid and transformed reference landmarks and spline coefficients
    yTPS = NULL;
    c = NULL;
    yr = NULL;
    double t0, t1;
    
    //Read input
    const double *t = static_cast<double*>(mxGetData(prhs[0])); // template landmarks as Lxd matrix, L is number of landmarks, d is spatial dimension
    double *r = static_cast<double*>(mxGetData(prhs[1])); // reference landmarks as Lxd matrix, L is number of landmarks, d is spatial dimension
    double *xc = static_cast<double*>(mxGetData(prhs[2])); // cell centered grid which should be transformed, x is prod(m)*dim
    const double *ptheta = static_cast<double*>(mxGetData(prhs[3])); //regularization parameter enabling smoother solutions approximating the landmark constrains
    double *omega, *m;
    omega = NULL;
    m = NULL;
    
    const double theta = *ptheta;
    
    // Dimension equals number of columns of r
    const int dim = mxGetN(prhs[1]);
        
    // Number of landmarks equals number of rows of r
    const int L = mxGetM(prhs[1]);
        
    // length of xc
    const int N = mxGetM(prhs[2]);
    
    //if N = 1, we assume, that TPS should be evaluated on a regular grid (1 = cell-centered) and load omega and m
    if (N==1){
        omega = static_cast<double*>(mxGetData(prhs[4]));
        m     = static_cast<double*>(mxGetData(prhs[5]));
    }
    
    int temp = dim;
    mwSize dimsYTPS[2];
    if (N == 1){
        for (int i=0; i<dim; i++){
            temp *= m[i];
        }
        dimsYTPS[0] = (int) temp;
    }
    else{
        dimsYTPS[0] = N;
    }
    dimsYTPS[1] = 1;
    
    mwSize dimsYR[2];
    dimsYR[0] = L;
    dimsYR[1] = dim;
    
    mwSize dimsC[2];
    dimsC[0] = L + dim +1;
    dimsC[1] = dim;
    
    plhs[0] = mxCreateNumericArray(2, dimsYTPS, mxDOUBLE_CLASS, mxREAL);
    yTPS = static_cast<double*>(mxGetData(plhs[0]));
    
    plhs[1] = mxCreateNumericArray(2, dimsYR, mxDOUBLE_CLASS, mxREAL);
    yr = static_cast<double*>(mxGetData(plhs[1]));
    
    
    
    plhs[2] = mxCreateNumericArray(2, dimsC, mxDOUBLE_CLASS, mxREAL);
    c = static_cast<double*>(mxGetData(plhs[2]));
    //t0 = getCurrentTime();
    ComputeCoefficients(c, t, r, dim, L, theta);
    //t1 = getCurrentTime();
    //mexPrintf("runTimeSolve=%f \n",t1-t0);
    
    if (N==1){
        //t0 = getCurrentTime();
        EvaluateTPScc(yTPS, yr, c, omega, m, r, dim, L);
        //t1 = getCurrentTime();
        //mexPrintf("runTimeGrid with use of adhoc computed grid=%f \n",t1-t0);
    }
    else{
        //t0 = getCurrentTime();
        EvaluateTPS(yTPS, yr, c, xc, r, dim, L, N);
        //t1 = getCurrentTime();
        //mexPrintf("runTimeGrid with use of predefined Grid=%f \n",t1-t0);
    }
}

void ComputeCoefficients(double* c, const double* t, double* r, const int dim, const int L, const double theta){
    const int numCols  = L+dim+1;
    double* A = new double[numCols*numCols];
    double* rhs = new double[numCols];
    for (int i=L; i < numCols; i++){
        rhs[i] = 0;
    }
    //double t0 = getCurrentTime();
    switch (dim){
        
        case 2:
            //A-Block
            for (int i=0; i < L; i++){
                A[i+i*numCols] = theta;
                for (int j=0; j < i; j++){
                    double res1 = r[i]-r[j];
                    double res2 = r[i+L]-r[j+L];
                    double norm = sqrt(res1*res1 + res2*res2);
                    if (norm > 0){
                        double temp = norm*norm * log(norm);
                        A[i+j*numCols] = temp;
                        A[j+i*numCols] = temp;
                    }
                    else{
                        A[i+j*numCols]=0;
                        A[j+i*numCols]=0;
                    }
                }
            }
            break;
        case 3:
            //A-Block
            for (int i=0; i < L; i++){
                A[i+i*numCols] = theta;
                for (int j=0; j < i; j++){
                    double res1 = r[i]-r[j];
                    double res2 = r[i+L]-r[j+L];
                    double res3 = r[i+2*L]-r[j+2*L];
                    double norm = sqrt(res1*res1 + res2*res2 + res3*res3);
                    A[i+j*numCols] = norm;
                    A[j+i*numCols] = norm;
                }
            }
            break;
            
        default:
            mexPrintf("Only 2D and 3D data supported!");
            return;
            
    }
    //B-Block
    for (int i=0;i < L; i++){
        A[i+L*numCols] = 1;
        A[L+i*numCols] = 1;
        for (int j=0; j < dim; j++){
            A[i+(L+j+1)*numCols] = r[i+j*L];
            A[L+j+1+i*numCols] = r[i+j*L];
        }
    }
    
    //0-Block
    for (int i=L;i < numCols;i++){
        for (int j=L;j<numCols;j++){
            A[j+i*numCols] = 0;
        }
    }
    //double t1 = getCurrentTime();
    //mexPrintf("Time for Building A: %f \n",t1-t0);
    //Solve system
    for (int k=0; k<dim; k++){
        for (int i=0; i < L; i++){
            rhs[i] = t[i+k*L];
        }
        //double t0 = getCurrentTime();
        solve(A,rhs,c+k*numCols,numCols);
        //double t1 = getCurrentTime();
        //mexPrintf("Time Solving %d LES: %f \n",k,t1-t0);
    }
       
    delete[] rhs;
    delete[] A;
}

void EvaluateTPS(double* yTPS, double* yr, double* c, double* xc, double* r, const int dim, const int L, const int N){
    const int prodm = N/dim;
    
    const int cstride = dim+1+L;
    
    
    switch(dim){
        case 2:
            #pragma omp parallel for
            for(int k=0; k<prodm; k++){
                double x1 = xc[k];
                double x2 = xc[k+prodm];
                
                //Polynomial part
                yTPS[k]         = c[L]           + x1 * c[L+1]             + x2 * c[L+2];
                yTPS[k+prodm]   = c[L+cstride]   + x1 * c[cstride + L+1]   + x2 * c[cstride + L+2];
                
                for(int j=0; j < L; j++){
                    //adding weighted radial basis functions
                    double res1 = x1 - r[j];
                    double res2 = x2 - r[j+L];
                    
                    double norm = sqrt(res1*res1 + res2*res2);
                    if (norm > 0){
                        yTPS[k]           += c[j]           * norm * norm * log(norm);
                        yTPS[prodm + k]   += c[j+cstride]   * norm * norm * log(norm);
                    }
                    
                }
                    } // end of k-loop
            #pragma omp parallel for
            for(int j=0; j < L; j++){
                double r1 = r[j];
                double r2 = r[j+L];
                
                //Polynomial part
                yr[j]       = c[L]           + r1 * c[L+1]             + r2 * c[L+2];
                yr[j + L]   = c[L+cstride]   + r1 * c[cstride + L+1]   + r2 * c[cstride + L+2];
                
                for(int l=0; l<L; l++){
                    double res1 = r[l]     - r1;
                    double res2 = r[l+L]   - r2;
                    double norm = sqrt(res1*res1 + res2*res2);
                    if (norm > 0){
                        yr[j]       += c[l]           * norm * norm * log(norm);
                        yr[L + j]   += c[l+cstride]   * norm * norm * log(norm);
                    }
                    
                }
            }
            break;
        case 3:
            #pragma omp parallel for
            for(int k=0; k<prodm; k++){
                double x1 = xc[k];
                double x2 = xc[k+prodm];
                double x3 = xc[k+2*prodm];
    
                //Polynomial part
                yTPS[k]         = c[L]           + x1 * c[L+1]             + x2 * c[L+2]             + x3 * c[L+3];
                yTPS[k+prodm]   = c[L+cstride]   + x1 * c[cstride + L+1]   + x2 * c[cstride + L+2]   + x3 * c[cstride + L+3];
                yTPS[k+2*prodm] = c[L+2*cstride] + x1 * c[2*cstride + L+1] + x2 * c[2*cstride + L+2] + x3 * c[2*cstride + L+3];
    
                for(int j=0; j < L; j++){
                //adding weighted radial basis functions
                    double res1 = x1 - r[j];
                    double res2 = x2 - r[j+L];
                    double res3 = x3 - r[j+2*L];
        
                    double norm = sqrt(res1*res1 + res2*res2 + res3*res3);
        
                    yTPS[k]           += c[j]           * norm;
                    yTPS[prodm + k]   += c[j+cstride]   * norm;
                    yTPS[2*prodm + k] += c[j+2*cstride] * norm;
        
                }
            } // end of k-loop
            #pragma omp parallel for
            for(int j=0; j < L; j++){
                double r1 = r[j];
                double r2 = r[j+L];
                double r3 = r[j+2*L];
    
                //Polynomial part
                yr[j]       = c[L]           + r1 * c[L+1]             + r2 * c[L+2]             + r3 * c[L+3];
                yr[j + L]   = c[L+cstride]   + r1 * c[cstride + L+1]   + r2 * c[cstride + L+2]   + r3 * c[cstride + L+3];
                yr[j + 2*L] = c[L+2*cstride] + r1 * c[2*cstride + L+1] + r2 * c[2*cstride + L+2] + r3 * c[2*cstride + L+3];
    
                for(int l=0; l<L; l++){
                    double res1 = r[l]     - r1;
                    double res2 = r[l+L]   - r2;
                    double res3 = r[l+2*L] - r3;
                    double norm = sqrt(res1*res1 + res2*res2 + res3*res3);
                    yr[j]       += c[l]           * norm;
                    yr[L + j]   += c[l+cstride]   * norm;
                    yr[2*L + j] += c[l+2*cstride] * norm;
                }
            }
            break;
        default:
            mexPrintf("Only 2D and 3D data supported!");
            break;
    }
} 

void EvaluateTPScc(double* yTPS, double* yr, double* c, double* omega, double* m, double* r, const int dim, const int L){
    
    int prodm = 1;
    for (int i=0;i<dim;i++){
        prodm *= (int) m[i];
    }
    
    const int cstride = dim+1+L;
    const int m0 = (int) m[0];
    const int yStride = m0;
    const int m1 = (int) m[1];
    
    const double omega0 = omega[0];
    const double omega2 = omega[2];
    
    double omega4 = 0;
    
    const double h0 = (omega[1] - omega[0]) / m[0];
    const double h1 = (omega[3] - omega[2]) / m[1];
    int m2=0, zStride=0;
    
    double h2=0;
    
    if (dim > 2){
        m2 = (int) m[2];
        omega4 = omega[4];
        zStride = m0*m1;
        h2 = (omega[5] - omega[4]) / m[2];
    }
   
    switch(dim){
        case 2:
            #pragma omp parallel for
            for(int j=0; j<m1; j++){
                double x2 = omega2 + (j+0.5)*h1;
                for(int i=0; i<m0; i++){
                   int index = i + yStride*j;
                   double x1 = omega0 + (i+0.5)*h0;
                        
                   //Polynomial part
                   yTPS[index]         = c[L]           + x1 * c[L+1]             + x2 * c[L+2];
                   yTPS[index+prodm]   = c[L+cstride]   + x1 * c[cstride + L+1]   + x2 * c[cstride + L+2];
            
                   for(int l=0; l < L; l++){
                        //adding weighted radial basis functions
                            double res1 = x1 - r[l];
                            double res2 = x2 - r[l+L];
                                                
                            double norm = sqrt(res1*res1 + res2*res2);
        
                            yTPS[index]           += c[l]           * norm;
                            yTPS[prodm + index]   += c[l+cstride]   * norm;
                     
                        
                    }
                }
            }        // end of j-loop
            #pragma omp parallel for
            for(int j=0; j < L; j++){
                double r1 = r[j];
                double r2 = r[j+L];
                
                //Polynomial part
                yr[j]       = c[L]           + r1 * c[L+1]             + r2 * c[L+2];
                yr[j + L]   = c[L+cstride]   + r1 * c[cstride + L+1]   + r2 * c[cstride + L+2];
                
                for(int l=0; l<L; l++){
                    double res1 = r[l]     - r1;
                    double res2 = r[l+L]   - r2;
                    double norm = sqrt(res1*res1 + res2*res2);
                    if (norm > 0){
                        yr[j]       += c[l]           * norm * norm * log(norm);
                        yr[L + j]   += c[l+cstride]   * norm * norm * log(norm);
                    }
                    
                }
            } 
        break;
        
        case 3:
                                   
            #pragma omp parallel for
            for(int k=0; k<m2; k++){
                double x3 = omega4 + (k+0.5)*h2;
                for(int j=0; j<m1; j++){
                    double x2 = omega2 + (j+0.5)*h1;
                    for(int i=0; i<m0; i++){
                        int index = i + yStride*j + zStride*k;
                        double x1 = omega0 + (i+0.5)*h0;
                        
                        //Polynomial part
                        yTPS[index]         = c[L]           + x1 * c[L+1]             + x2 * c[L+2]             + x3 * c[L+3];
                        yTPS[index+prodm]   = c[L+cstride]   + x1 * c[cstride + L+1]   + x2 * c[cstride + L+2]   + x3 * c[cstride + L+3];
                        yTPS[index+2*prodm] = c[L+2*cstride] + x1 * c[2*cstride + L+1] + x2 * c[2*cstride + L+2] + x3 * c[2*cstride + L+3];
    
                        for(int l=0; l < L; l++){
                        //adding weighted radial basis functions
                            double res1 = x1 - r[l];
                            double res2 = x2 - r[l+L];
                            double res3 = x3 - r[l+2*L];
                                    
                            double norm = sqrt(res1*res1 + res2*res2 + res3*res3);
        
                            yTPS[index]           += c[l]           * norm;
                            yTPS[prodm + index]   += c[l+cstride]   * norm;
                            yTPS[2*prodm + index] += c[l+2*cstride] * norm;
                        }
                    }
                }
            }        // end of k-loop
            
            
            
             
            #pragma omp parallel for
            for(int j=0; j < L; j++){
                double r1 = r[j];
                double r2 = r[j+L];
                double r3 = r[j+2*L];
    
                //Polynomial part
                yr[j]       = c[L]           + r1 * c[L+1]             + r2 * c[L+2]             + r3 * c[L+3];
                yr[j + L]   = c[L+cstride]   + r1 * c[cstride + L+1]   + r2 * c[cstride + L+2]   + r3 * c[cstride + L+3];
                yr[j + 2*L] = c[L+2*cstride] + r1 * c[2*cstride + L+1] + r2 * c[2*cstride + L+2] + r3 * c[2*cstride + L+3];
    
                for(int l=0; l<L; l++){
                    double res1 = r[l]     - r1;
                    double res2 = r[l+L]   - r2;
                    double res3 = r[l+2*L] - r3;
                    double norm = sqrt(res1*res1 + res2*res2 + res3*res3);
                    yr[j]       += c[l]           * norm;
                    yr[L + j]   += c[l+cstride]   * norm;
                    yr[2*L + j] += c[l+2*cstride] * norm;
                }
            } 
        break;
        default:
            mexPrintf("Only 2D and 3D data supported!");
            break;
    }    
    
}

int solve(double* LS, const double* rhs, double* x, const int m) {
    
    // check pointers
    if ((LS == NULL) || (rhs == NULL) || (x == NULL)) {
        return NullPointerError;
    }
    
    // allocate mem and init
    double* A = new double[m*m];
    if (A == NULL) {
        return NullPointerError;
    }
    for (int i=0; i<m*m; ++i) {
        A[i]=LS[i];
    }
    
    double* b = new double[m];
    if (b == NULL) {
        delete[] A;
        return NullPointerError;
    }
    for (int i=0; i<m; ++i){
        b[i]=rhs[i];
    }
    
    const double tol = 1e-15;
    double c = 0;
    int i = 0, j = 0, k = 0, tmp = 0;
    
    int end = m-1;
    int ret = NoError;
    
    int* I = new int[m];
    
    if (I == NULL) {
        delete[] A;
        delete[] b;
        return NullPointerError;
    }
    for (i=0; i<=end; ++i) {
        I[i] = i;
    }
    // end of allocation
    
    for (i=0; i<end ; ++i) {
        // find pivot
        for (k=i+1; k<=end ; ++k) {
            if ( fabs(A[I[k]+m*i]) > fabs(A[I[i]+m*i])) {
                tmp  = I[i];
                I[i] = I[k];
                I[k] = tmp;
            }
        }
        // Absolute value of pivot smaller than tol. System nearly singular.
        if (fabs(A[I[i]+m*i]) < tol) {
            mexPrintf("Absolute value of pivot smaller than tol. System nearly singular.");
            ret = NumericalError;
            break;
        }
        
        for (k=i+1; k<=end ; ++k) {
            c = A[I[k]+m*i] / A[I[i]+m*i];
            for (j=i; j<=end ; ++j) {
                A[I[k]+m*j] -= c*A[I[i]+m*j];
            }
            b[I[k]]  -=  c*b[I[i]];
        }
    }
    
    if (ret == NoError) {
        // substitute
        x[end] = b[I[end]] / A[I[end]+m*end];
        for (i=end-1; i>=0; --i) {
            x[i] = b[I[i]];
            for (j=i+1; j<=end; ++j) {
                x[i] -=  A[I[i]+m*j]*x[j];
            }
            x[i] /= A[I[i]+m*i];
        }
    } else {
        for (int i = 0; i < end; ++i) {
            x[i] = 0;
        }
    }
    
    delete[] A;
    delete[] b;
    delete[] I;
    return ret;
    
}


