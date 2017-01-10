#include <math.h>
#include <mex.h>
#include <algorithm>
#include <stdio.h>
#include <iostream>

using namespace std;

void diagBB2D(double* D, const double* omega, const double* m, const double* alpha);
void diagBB3D(double* D, const double* omega, const double* m, const double* alpha);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  
  //Initialize output
  double *D;
  
  D = NULL;
  
  //Read input
  const double *omega = static_cast<double*>(mxGetData(prhs[0]));
  const double *m = static_cast<double*>(mxGetData(prhs[1]));
  const double *alpha = static_cast<double*>(mxGetData(prhs[2]));
  
  // Dimension equals number of columns of m
  int dim = mxGetN(prhs[1]);
  
  // N = length of vector uc
  int N = dim;
  for (int i=0; i<dim;i++){
    N *= m[i];
  }
   
  mwSize dimsD[2]; dimsD[0] = N; dimsD[1] = 1;
    
  if (nrhs != 3){
    mexErrMsgTxt("Number of arguments must be 3");
  }
  
  // D = alpha * hd * diag(B' *B)  
  plhs[0] = mxCreateNumericArray(2, dimsD, mxDOUBLE_CLASS, mxREAL);
  D = static_cast<double*>(mxGetData(plhs[0]));
    
  switch (dim) {
    case 2:
      diagBB2D(D, omega, m, alpha);
      break;
    case 3:
      diagBB3D(D, omega, m, alpha);
      break;
    default:
      break;
  }
}

void diagBB2D(double* D, const double* omega, const double* m, const double* alpha){
  
  const int m0 = (int) m[0];
  const int m1 = (int) m[1];
  const int prodm = m0 * m1;
   
  const double h0 = (omega[1] - omega[0]) / m[0];
  const double h1 = (omega[3] - omega[2]) / m[1];
  const double coeff = alpha[0] * h0 * h1;
  
  const double invH2_0 = 1.0/(h0*h0);
  const double invH2_1 = 1.0/(h1*h1);
  const double sumInvH2 = -2*(invH2_0 + invH2_1);
  
  
  const int yStride = m0;
    
  #pragma omp parallel for
  for (int j=0; j<m1; j++){
      const int j1 = max(0,j-1);
      const int j2 = min(m1-1,j+1);
      for (int i=0; i<yStride; i++){
        double buffer[5];
        const int index  = i                  + j*yStride;
        const int indexL = max(0,        i-1) + j*yStride;
        const int indexR = min(yStride-1,i+1) + j*yStride;
        
        const int indexO = i                  + j1*yStride;
        const int indexU = i                  + j2*yStride;
               
        buffer[0] = 0;
        buffer[1] = 0;
        buffer[2] = sumInvH2;
        buffer[3] = 0;
        buffer[4] = 0;
        
        
        buffer[2-(index-indexL)] += invH2_0;
        buffer[2+(indexR-index)] += invH2_0;
        buffer[2-2*(index-indexO)/yStride] += invH2_1;
        buffer[2+2*(indexU-index)/yStride] += invH2_1;
        
        
        D[index] = coeff * (buffer[0]*buffer[0] + buffer[1]*buffer[1] + buffer[2]*buffer[2] + buffer[3]*buffer[3] + buffer[4]*buffer[4]);
        D[index+prodm]   = D[index];
                  
     }  
  }
}

void diagBB3D(double* D, const double* omega, const double* m, const double* alpha){
  
  const int m0 = (int) m[0];
  const int m1 = (int) m[1];
  const int m2 = (int) m[2];
  const int prodm = m0 * m1 * m2;
   
  const double h0 = (omega[1] - omega[0]) / m[0];
  const double h1 = (omega[3] - omega[2]) / m[1];
  const double h2 = (omega[5] - omega[4]) / m[2];
  const double coeff = alpha[0] * h0 * h1 * h2;
  
  const double invH2_0 = 1.0/(h0*h0);
  const double invH2_1 = 1.0/(h1*h1);
  const double invH2_2 = 1.0/(h2*h2);
  const double sumInvH2 = -2*(invH2_0 + invH2_1 + invH2_2);
  
  
  const int yStride = m0;
  const int zStride = m0*m1;
  
  #pragma omp parallel for
  for (int k=0; k<m2; k++){
    const int k1 = max(0,k-1);
    const int k2 = min(m2-1,k+1);
    for (int j=0; j<m1; j++){
      const int j1 = max(0,j-1);
      const int j2 = min(m1-1,j+1);
      for (int i=0; i<yStride; i++){
        double buffer[7];
        const int index  = i                  + j*yStride             + k*zStride;
        const int indexL = max(0,        i-1) + j*yStride             + k*zStride;
        const int indexR = min(yStride-1,i+1) + j*yStride             + k*zStride;
        
        const int indexO = i                  + j1*yStride + k*zStride;
        const int indexU = i                  + j2*yStride + k*zStride;
        
        const int indexV = i                  + j*yStride             + k1*zStride;
        const int indexH = i                  + j*yStride             + k2*zStride;
        
        buffer[0] = 0;
        buffer[1] = 0;
        buffer[2] = 0;
        buffer[3] = sumInvH2;
        buffer[4] = 0;
        buffer[5] = 0;
        buffer[6] = 0;
        
        buffer[3-(index-indexL)] += invH2_0;
        buffer[3+(indexR-index)] += invH2_0;
        buffer[3-2*(index-indexO)/yStride] += invH2_1;
        buffer[3+2*(indexU-index)/yStride] += invH2_1;
        buffer[3-3*(index-indexV)/zStride] += invH2_2;
        buffer[3+3*(indexH-index)/zStride] += invH2_2;
        
        D[index] = coeff * (buffer[0]*buffer[0] + buffer[1]*buffer[1] + buffer[2]*buffer[2] + buffer[3]*buffer[3] + buffer[4]*buffer[4] + buffer[5]*buffer[5] + buffer[6]*buffer[6]);
        D[index+prodm]   = D[index];
        D[index+2*prodm] = D[index];
          
      }  
    }
  }
}