/*
 * (c) Fabian Gigengack 2011/04/13 see FAIR.2 and FAIRcopyright.m.
 * http://www.uni-muenster.de/EIMI/
 *
 * C code for kron(Q3,Q2,Q1) * w. See splineTransformation3Dsparse.
 */

#include "mex.h"
#ifdef _OPENMP
    #include<omp.h>
#endif

void mexFunction(int nlhs, mxArray *plhs[ ],int nrhs, const mxArray *prhs[ ])
{
    double *Qw,*w,*Q1,*Q2,*Q3,*p,*szW;
    int i,k,i1,i2,i3,j1,j2,j3;
    int sz,szQ1,szQ2,szQ3,szW1,szW2,szW3,szQ1Q2,szQ1Q2Q3,szW1W2,szW1W2W3;
    
    if (nrhs != 6){
        mexErrMsgTxt("Six input arguments required.");
    } else if (nlhs > 1) {
        mexErrMsgTxt("Too many output arguments.");
    }
    
    w   =  mxGetPr(prhs[0]);
    Q1  =  mxGetPr(prhs[1]);
    Q2  =  mxGetPr(prhs[2]);
    Q3  =  mxGetPr(prhs[3]);
    p   =  mxGetPr(prhs[4]);
    szW =  mxGetPr(prhs[5]);
    
    szQ1 = (int) p[0];
    szQ2 = (int) p[1];
    szQ3 = (int) p[2];
    szW1 = (int) szW[0];
    szW2 = (int) szW[1];
    szW3 = (int) szW[2];
    
    szQ1Q2   = szQ1 * szQ2;
    szQ1Q2Q3 = szQ1 * szQ2 * szQ3;
    szW1W2   = szW1 * szW2;
    szW1W2W3 = szW1 * szW2 * szW3;
    sz       = szQ1Q2Q3 * 3;
    
    plhs[0] = mxCreateDoubleMatrix(sz, 1, mxREAL);
    Qw = mxGetPr(plhs[0]);
    
    /* initialization */
    for (i=0;i<sz;i++){
        Qw[i] = 0;
    }
    /* compute tensor product */
    #pragma omp parallel for default(shared)
    for (k=0;k<3;k++){
        for (i3=0;i3<szQ3;i3++){
            for (i2=0;i2<szQ2;i2++){
                for (i1=0;i1<szQ1;i1++){
                    for (j3=0;j3<szW3;j3++){
                        for (j2=0;j2<szW2;j2++){
                            for (j1=0;j1<szW1;j1++){
                                Qw[i1+i2*szQ1+i3*szQ1Q2+k*szQ1Q2Q3] +=
                                        Q1[i1+j1*szQ1] *
                                        Q2[i2+j2*szQ2] *
                                        Q3[i3+j3*szQ3] *
                                        w[j1+j2*szW1+j3*szW1W2+k*szW1W2W3];
                            }
                        }
                    }
                }
            }
        }
    }
}

/* ==================================================================================== */
