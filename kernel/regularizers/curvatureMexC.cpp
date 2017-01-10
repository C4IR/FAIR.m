//author: Thomas Polzin
//allows for matrixfree Computation of Curvature regularization and its 1st derivative.

#include <math.h>
#include <mex.h>
#include <algorithm>
#include <stdio.h>
#include <iostream>

using namespace std;

void Curvature3D(double* Sc, double* dS, double* xc, const double* omega, const double* m, const double* alpha);
void Curvature2D(double* Sc, double* dS, double* xc, const double* omega, const double* m, const double* alpha);
void BTimesx2D(double* Bx, double* xc, const double* omega, const double* m);
void BTimesx3D(double* Bx, double* xc, const double* omega, const double* m);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  
  //Ausgaben initialisieren
  double *Sc, *dS; 
  dS = NULL;
  

  
  //Eingaben einlesen
  double *uc = static_cast<double*>(mxGetData(prhs[0]));
  const double *omega = static_cast<double*>(mxGetData(prhs[1]));
  const double *m = static_cast<double*>(mxGetData(prhs[2]));
  const double *alpha = static_cast<double*>(mxGetData(prhs[3]));
  
  // Dimension gleich Anzahl Spalten von m 
  int dim = mxGetN(prhs[2]);
  
  // Laenge des betrachteten Vektors uc
  const int N = mxGetM(prhs[0]);
  
  mwSize dims[2]; 
  dims[0] = 1; 
  dims[1] = N;
    
  if (nrhs < 4){
    mexErrMsgTxt("Number of arguments must be 4");
  }
  
  // Sc als Funktionswert
  plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
  Sc = static_cast<double*>(mxGetData(plhs[0]));
  Sc[0] = 0; 
  
  plhs[1] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
  dS = static_cast<double*>(mxGetData(plhs[1]));
    
  switch (dim) {
    case 2:
      Curvature2D(Sc,dS,uc,omega,m,alpha);
      break;
    case 3:
      Curvature3D(Sc,dS,uc,omega,m,alpha);
      break;
    default:
      break;
  }
}

void Curvature2D(double* Sc, double* dS, double* xc, const double* omega, const double* m, const double* alpha){
  
  double result = 0;
  const int m0 = (int) m[0];
  const int m1 = (int) m[1];
  const int prodm = m0 * m1;
  
  const double h0 = (omega[1] - omega[0]) / m[0];
  const double h1 = (omega[3] - omega[2]) / m[1];
  const double coeff = alpha[0] * h0 * h1;
  
  double* Bx = new double[2*prodm];
  
  BTimesx2D(Bx, xc, omega, m);
  BTimesx2D(dS, Bx, omega, m);
  
#pragma omp parallel for reduction(+: result)
  for (int i=0; i < 2 * prodm; i++){
    dS[i] = dS[i] * coeff;
    result += 0.5 * dS[i] * xc[i];
  }
  
  Sc[0] = result;

// free dynamically allocated buffer
  delete[] Bx;
}


void Curvature3D(double* Sc, double* dS, double* xc, const double* omega, const double* m, const double* alpha){
  
  double result = 0;
  const int m0 = (int) m[0];
  const int m1 = (int) m[1];
  const int m2 = (int) m[2];
  const int prodm = m0 * m1 * m2;
  
  const double h0 = (omega[1] - omega[0]) / m[0];
  const double h1 = (omega[3] - omega[2]) / m[1];
  const double h2 = (omega[5] - omega[4]) / m[2];
  const double coeff = alpha[0] * h0 * h1 * h2;
  
  double* Bx = new double[3*prodm];
  
  BTimesx3D(Bx, xc, omega, m);
  BTimesx3D(dS, Bx, omega, m);
  
#pragma omp parallel for reduction(+: result)
  for (int i=0; i < 3 * prodm; i++){
    dS[i] = dS[i] * coeff;
    result += 0.5 * dS[i] * xc[i];
  }
  
  Sc[0] = result;

// free dynamically allocated buffer
  delete[] Bx;
}

void BTimesx2D(double* Bx, double* xc, const double* omega, const double* m){
// Code adaptiert nach Jan Ruehaak und Till Kipshagen
// xc als uc
  const int m0 = (int) m[0];
  const int m1 = (int) m[1];
  const int prodm = m0 * m1;
  
  const double h0 = (omega[1] - omega[0]) / m[0];
  const double h1 = (omega[3] - omega[2]) / m[1];
    
  const double invH2_0 = 1.0/(h0*h0);
  const double invH2_1 = 1.0/(h1*h1);
  const double sumInvH2 = -2*(invH2_0 + invH2_1);
  
  double * xc_y = xc + prodm; // access to y-component
  double * Bx_y = Bx + prodm;
    
  const int yStride = m0;
    
#pragma omp parallel 
{

#pragma omp for
    for (int j=0; j<m1; j++){
      for (int i=0; i<m0; i++){

        const int index = i + j*yStride;
        const int indexL = max(0,     i-1) + j*yStride;
        const int indexR = min(m0-1,i+1) + j*yStride ;

        const int indexO = i + max(0,     j-1)*yStride;
        const int indexU = i + min(m1-1,j+1)*yStride;

        Bx[index] = invH2_0*(xc[indexL] + xc[indexR])
        + invH2_1*(xc[indexO] + xc[indexU])
        +sumInvH2*xc[index];
        
      }
    }
  

#pragma omp for
    for (int j=0; j<m1; j++){
      for (int i=0; i<m0; i++){

        const int index = i + j*yStride;
        const int indexL = max(0,     i-1) + j*yStride;
        const int indexR = min(m0-1,i+1) + j*yStride;

        const int indexO = i + max(0,     j-1)*yStride;
        const int indexU = i + min(m1-1,j+1)*yStride;

        Bx_y[index] = invH2_0*(xc_y[indexL] + xc_y[indexR])
        + invH2_1*(xc_y[indexO] + xc_y[indexU])
        +sumInvH2*xc_y[index];
      }
    }
  

 }
  return;

}

void BTimesx3D(double* Bx, double* xc, const double* omega, const double* m){
// Code adaptiert nach Jan Ruehaak und Till Kipshagen
// xc als uc
  const int m0 = (int) m[0];
  const int m1 = (int) m[1];
  const int m2 = (int) m[2];
  const int prodm = m0 * m1 * m2;
  
  const double h0 = (omega[1] - omega[0]) / m[0];
  const double h1 = (omega[3] - omega[2]) / m[1];
  const double h2 = (omega[5] - omega[4]) / m[2];
  
  const double invH2_0 = 1.0/(h0*h0);
  const double invH2_1 = 1.0/(h1*h1);
  const double invH2_2 = 1.0/(h2*h2);
  const double sumInvH2 = -2*(invH2_0 + invH2_1 + invH2_2);
  
  double * xc_y = xc + prodm; // access to y-component
  double * Bx_y = Bx + prodm;
  double * xc_z = xc + 2*prodm; // access to z-component
  double * Bx_z = Bx + 2*prodm;
  
  const int yStride = m0;
  const int zStride = m0*m1;
  
#pragma omp parallel 
{

#pragma omp for
  for (int k=0; k<m2; k++){
    for (int j=0; j<m1; j++){
      for (int i=0; i<m0; i++){

        const int index = i + j*yStride + k*zStride;
        const int indexL = max(0,     i-1) + j*yStride + k*zStride;
        const int indexR = min(m0-1,i+1) + j*yStride + k*zStride;

        const int indexO = i + max(0,     j-1)*yStride + k*zStride;
        const int indexU = i + min(m1-1,j+1)*yStride + k*zStride;

        const int indexV = i + j*yStride + max(0,     k-1)*zStride;
        const int indexH = i + j*yStride + min(m2-1,k+1)*zStride;

        Bx[index] = invH2_0*(xc[indexL] + xc[indexR])
        + invH2_1*(xc[indexO] + xc[indexU])
        + invH2_2*(xc[indexV] + xc[indexH])
        +sumInvH2*xc[index];
        
      }
    }
  }

#pragma omp for
  for (int k=0; k<m2; k++){
    for (int j=0; j<m1; j++){
      for (int i=0; i<m0; i++){

        const int index = i + j*yStride + k*zStride;
        const int indexL = max(0,     i-1) + j*yStride + k*zStride;
        const int indexR = min(m0-1,i+1) + j*yStride + k*zStride;

        const int indexO = i + max(0,     j-1)*yStride + k*zStride;
        const int indexU = i + min(m1-1,j+1)*yStride + k*zStride;

        const int indexV = i + j*yStride + max(0,     k-1)*zStride;
        const int indexH = i + j*yStride + min(m2-1,k+1)*zStride;

        Bx_y[index] = invH2_0*(xc_y[indexL] + xc_y[indexR])
        + invH2_1*(xc_y[indexO] + xc_y[indexU])
        + invH2_2*(xc_y[indexV] + xc_y[indexH])
        +sumInvH2*xc_y[index];
      }
    }
  }

#pragma omp for
  for (int k=0; k<m2; k++){
    for (int j=0; j<m1; j++){
      for (int i=0; i<m0; i++){

        const int index = i + j*yStride + k*zStride;
        const int indexL = max(0,     i-1) + j*yStride + k*zStride;
        const int indexR = min(m0-1,i+1) + j*yStride + k*zStride;

        const int indexO = i + max(0,     j-1)*yStride + k*zStride;
        const int indexU = i + min(m1-1,j+1)*yStride + k*zStride;

        const int indexV = i + j*yStride + max(0,     k-1)*zStride;
        const int indexH = i + j*yStride + min(m2-1,k+1)*zStride;

        Bx_z[index] = invH2_0*(xc_z[indexL] + xc_z[indexR])
        + invH2_1*(xc_z[indexO] + xc_z[indexU])
        + invH2_2*(xc_z[indexV] + xc_z[indexH])
        +sumInvH2*xc_z[index];

      }
    }
  }
}
  return;
}


