/*
 * (c) Jan Modersitzki 2010/12/27, see FAIR.2 and FAIRcopyright.m
 * http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
 *
 * function [rho,drho] = rhoSplineC(T,R,minT,maxT,nT,minR,maxR,nR,doDerivative)
 * 
 * Parzen-Window Based density estimator for gray values of T and R, spline based.
 * 
 * 
 * Input: 
 *  T, R                template and reference
 *  minT, maxT, nT     discreted points in gray value range where splines are located
 *  minR, maxR, nR     discreted points in gray value range where splines are located
 *  doSerivative        FLAG for derivative compuatation
 * 
 * Output:
 *  rho          joint density estimator
 *  drho         derivative of joint density estimator
 * See thoSpline.m for more details
 */

#include "mex.h"
#include "matrix.h"
#include <math.h>

// sparse matrix output for the derivative; see MATLAB format
typedef struct {
  double* S;
  mwIndex* I;
  mwIndex* J;
  mwSize nzmax;
  mwSize rows;  
  mwSize cols;  
} sparseM;


// ========================================================================
// Here is where all the fun is.
// ========================================================================

void dpdfestimate(
  double*   TD,             // pointer to discretized image TD  
  double*   RD,             // pointer to discretized image RD  
  mwSize    m,              // length of TD and RD
  double    minT,           // first bin in TD  
  double    maxT,           // last bin in TD  
  mwSize    nT,             // number of splines for TD  
  double    minR,           // first bin in RD  
  double    maxR,           // last bin in RD  
  mwSize    nR,             // number of splines for RD  
  double*   rho,            // output: joint density
  sparseM   drho            // output: derivative of joint density
) 
{
  mwIndex   k;              // loop indices
  mwSize    i,j;            // loop indices
  int       p;              // a value is devide into floor p and remainder
  double    tk, rk;         // gray values of TD and RD
  int       It[4], Ir[4];   // indices for contributing splines
  double    st[4], sr[4];   // spline values
  double    dst[4];         // spline derivative values
  mwIndex   nzcnt = 0;      // counter for non-zero entires

  double widthT = (maxT-minT)/((double) nT);
  double widthR = (maxR-minR)/((double) nR);
  double facw     = 1.0/((double) m);
  double facdw    = 1.0/(((double) m)* widthT);

/*
  printf("bT(1)= %g,..,bT(%d) = %g (widthT = %g)\n",minT,nT,maxT,widthT);
  printf("bR(1)= %g,..,bT(%d) = %g (widthR = %g)\n",minR,nR,maxR,widthR);
/* */

  // initialize rho  
  for (i=0; i<nR*nT; ++i)
    rho[i] = 0.0;

  //nzcnt = 0;                    // counter for non-zero elements of drho

  for (k=0; k<m; ++k) {         // loop overt length of TD and RD

    // normalize widthT to 1
    tk  = (TD[k] - minT)/widthT + 0.5;
    p   = (int) floor(tk);
    tk  = tk - floor(tk);

    // indices of contibuting splines
    It[0] = p+2; 
    It[1] = p+1; 
    It[2] = p; 
    It[3] = p-1; 

    // spline values
    st[0] = 1.0/6.0 * tk*tk*tk;
    st[3] = 1.0/6.0 + 1.0/2.0*tk*(tk-1.0) - st[0];
    st[1] = tk - 2.0*st[0] + st[3];
    st[2] = 1.0 - st[0] - st[1] - st[3];

    // spline derivative values
    dst[0] = 0.5*tk*tk;
    dst[3] = -0.5 + tk - dst[0];
    dst[1] = 1 - 2*dst[0] + dst[3];
    dst[2] = -dst[0] - dst[1] - dst[3];

    // normalize widthR to 1
    rk = (RD[k] - minR)/widthR + 0.5;
    p  = (int) floor(rk);
    rk = rk -  floor(rk);

    // indices of contibuting splines
    Ir[0] = p+2; 
    Ir[1] = p+1; 
    Ir[2] = p; 
    Ir[3] = p-1; 
    
    // spline values
    sr[0] = 1.0/6.0 * rk*rk*rk;
    sr[3] = 1.0/6.0 + 1.0/2.0*rk*(rk-1.0) - sr[0];
    sr[1] = rk - 2.0*sr[0] + sr[3];
    sr[2] = 1.0 - sr[0] - sr[1] - sr[3];
    
    // check boundary points
    for (i=0; i<4; ++i) {
      
      if( (It[i] < 0) || (It[i] >= nT) ){
        It[i]  = 0;
        st[i]  = 0.0;
        dst[i] = 0.0;
      };

      if( (Ir[i] < 0) || (Ir[i] >= nR) ){
        Ir[i] = 0;
        sr[i] = 0.0;
      };
      
    };
    
    // blend the values
    drho.J[k] = nzcnt;
    for (i=0; i<4; ++i) {
      for (j=0; j<4; ++j) {
        rho[It[i]+nT*Ir[j]] += facw*st[i]*sr[j];
        drho.I[nzcnt]       = (mwIndex)(It[i]+nT*Ir[j]);
         drho.S[nzcnt]       = facdw*dst[i]*sr[j];
        nzcnt++;
      } 
    }
  }
  drho.J[m] = nzcnt;
}

// ========================================================================
// This is a copy of the prvious function but without the derivative part
// ========================================================================
void pdfestimate(
  double*   TD,             // pointer to discretized image TD  
  double*   RD,             // pointer to discretized image RD  
  mwSize    m,              // length of TD and RD
  double    minT,           // first bin in TD  
  double    maxT,           // last bin in TD  
  mwSize    nT,             // number of bin in TD  
  double    minR,           // first bin in RD  
  double    maxR,           // last bin in RD  
  mwSize    nR,             // number of bins in RD  
  double*   rho             // output: joint density
) 
{
  mwIndex   k,i,j;          // loop indices
  int       p;              // a gray value is devide into floor p and remainder
  double    tk, rk;         // gray values of TD and RD
  int       It[4], Ir[4];   // indices for contributing splines
  double    st[4], sr[4];   // spline values

  double widthT = (maxT-minT)/((double) nT);
  double widthR = (maxR-minR)/((double) nR);
  double facw     = 1.0/((double) m);
  double facdw    = 1.0/(((double) m)* widthT);

  // initialize rho  
  for( i=0; i<nR*nT; ++i ) rho[i] = 0.0;

  for( k=0; k<m; ++k ){// loop overt length of TD and RD

    // normalize widthT to 1
    tk  = (TD[k] - minT)/widthT + 0.5;
    p   = (int) floor(tk);
    tk  = tk - floor(tk);

    // indices of contibuting splines
    It[0] = p+2; 
    It[1] = p+1; 
    It[2] = p; 
    It[3] = p-1; 

    // spline values
    st[0] = 1.0/6.0 * tk*tk*tk;
    st[3] = 1.0/6.0 + 1.0/2.0*tk*(tk-1.0) - st[0];
    st[1] = tk - 2.0*st[0] + st[3];
    st[2] = 1.0 - st[0] - st[1] - st[3];

    // normalize widthR to 1
    rk = (RD[k] - minR)/widthR + 0.5;
    p  = (int) floor(rk);
    rk = rk -  floor(rk);

    // indices of contibuting splines
    Ir[0] = p+2; 
    Ir[1] = p+1; 
    Ir[2] = p; 
    Ir[3] = p-1; 
    
    // spline values
    sr[0] = 1.0/6.0 * rk*rk*rk;
    sr[3] = 1.0/6.0 + 1.0/2.0*rk*(rk-1.0) - sr[0];
    sr[1] = rk - 2.0*sr[0] + sr[3];
    sr[2] = 1.0 - sr[0] - sr[1] - sr[3];
    
    // check boundary points
    for (i=0; i<4; ++i) {
      
      if( (It[i] < 0) || (It[i] >= nT) ){
        It[i]  = 0;
        st[i]  = 0.0;
      };

      if( (Ir[i] < 0) || (Ir[i] >= nR) ){
        Ir[i] = 0;
        sr[i] = 0.0;
      };
      
    };

    // blend the values
    for( i=0; i<4; ++i ){
      for( j=0; j<4; ++j ){
        rho[It[i]+nT*Ir[j]] += facw*st[i]*sr[j];
      } 
    }
  };
}

// ========================================================================
// JM 2008/03/12 pdfestimate based on a Parzen Window estimator
// [rho,drho] = pdfestimateC(TD,RD,binT1,maxT,nT,binR1,binRN,nR);
// requires 8 input arguments:
// TD       discretized image T      RD       discretized image R
// minT     first bin in TD          minR     first bin in RD
// maxT     last  bin in TD          maxR     last  bin in RD
// nT       number of bins for TD    nR       number of bins for RD
// OUTPUT
// rho  
// drho
// ========================================================================

void mexFunction(
  int               nlhs,
  mxArray          *plhs[],
  int               nrhs,
  const mxArray    *prhs[]
)
{
                            // input variables
  double*   TD;             // pointer to discretized image TD  
  double*   RD;             // pointer to discretized image RD  
  double    minT;           // first bin in TD  
  double    maxT;           // last bin in TD  
  mwSize    nT;             // number of bin in TD  
  double    minR;           // first bin in RD  
  double    maxR;           // last bin in RD  
  mwSize    nR;             // number of bins in RD  

                            // output variables
  double*   rho;            // density estimator
  sparseM   drho;            // and its drivative    

                            // counter and such 
  mwSize m;


  if( nrhs != 8 ){
    mexErrMsgTxt("Need 8 input arguments.");
  };

  // the splines are located at the cell centered points of the bin discretiation
  // +--x--+--x--+--x--+--x--+
  // therefore, we have n contributors


  // transfer the input variables
  TD    = mxGetPr(prhs[0]);
  RD    = mxGetPr(prhs[1]);
  minT  = mxGetScalar(prhs[2]);
  maxT  = mxGetScalar(prhs[3]);
  nT    = (mwSize) mxGetScalar(prhs[4]); 
  minR  = mxGetScalar(prhs[5]);
  maxR  = mxGetScalar(prhs[6]);
  nR    = (mwSize) mxGetScalar(prhs[7]);

  // check whether RD and TD are of same length
  m = mxGetNumberOfElements(prhs[0]);
  if( mxGetNumberOfElements(prhs[1]) != m )
    mexErrMsgTxt("TD and RD must have same number of elements");
  
  // create the output variable rho
  plhs[0] = mxCreateDoubleMatrix(nT*nR,1,mxREAL);
  rho     = mxGetPr(plhs[0]);

  if( nlhs<2 ){
    //
    // compute without derivative
    //

    pdfestimate(TD,RD,m,minT,maxT,nT,minR,maxR,nR,rho);

  } else {  

    //
    // compute with derivative
    //
    // create the output variable rho
    plhs[1]    = mxCreateSparse(nT*nR, m,  16*m, mxREAL);
    drho.S     = mxGetPr(plhs[1]);
    drho.I     = mxGetIr(plhs[1]);
    drho.J     = mxGetJc(plhs[1]);
    drho.rows  = nT*nR;
    drho.cols  = m;
    drho.nzmax = 16*m;

    dpdfestimate(TD,RD,m,minT,maxT,nT,minR,maxR,nR,rho,drho);
  }
};

/* ========================================================================== */
