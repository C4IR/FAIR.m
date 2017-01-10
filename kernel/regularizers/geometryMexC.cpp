/*
 *=======================================================================================
 * (c) Lars Ruthotto 2011/01/11, see FAIR.2 and FAIRcopyright.m.
 *
 * MATLAB wrapper for matrix free hyperelastic regularization as
 * described in
 *
 * @article{2011-BMR,
 *    Author = {Burger M., Modersitzki J., Ruthotto L. },
 *    Publisher = {University of Muenster},
 *    Title = {A hyperelastic regularization energy for image registration},
 *    Year = {2011}
 *  }
 * =======================================================================================
 */
#include "geometryC.h"

void mexFunction(int nlhs, mxArray *plhs[ ], int nrhs, const mxArray *prhs[ ]) {
    double  *result, *dSvol;
    const double *alpha, *h, *ARef, *VRef, *x, *yc, *m, *dT, *T, *Jac;
    const double *yc2, *Jac1, *Jac2, *dI1, *dI2, *I1, *I2, *v1, *v2;
    int dim, mrows, ncols, n, nn, nTetra, nTriangle, i;
    char *flag;
    bool doDerivative;
    
    if (nrhs < 3) {
        mexErrMsgTxt("At least three inputs required.");
    } else if (nlhs > 2) {
        mexErrMsgTxt("Too many output arguments.");
    }
    
    yc = mxGetPr(prhs[0]);
    m  = mxGetPr(prhs[1]);
    /* reading dim  */
    mrows = mxGetM(prhs[1]);
    ncols = mxGetN(prhs[1]);
    dim = (int) mrows*ncols;
    switch (dim) {
        case 2:
            nTetra    = 4;
            nTriangle = 1;
            break;
        case 3:
            nTetra    = 24;
            nTriangle = 24;
            break;
        default:
            break;
    }
    /* reading flag */
    mwSize buflen = mxGetNumberOfElements(prhs[2])+1;
    /* Copy the string data into buf. */
    flag = (char*) mxCalloc(buflen, sizeof(mxChar));
    if (flag == NULL)
        mexErrMsgTxt("Not enough heap space to hold converted string.");
    mxGetString(prhs[2], flag, buflen);
    n  = 1;
    nn = 1;
    for (i=0; i<dim; i++) {
        n  *= (int)  m[i];
        nn *= (int) (m[i]+1);
    }
    mrows = mxGetM(prhs[0]);
    ncols = mxGetN(prhs[0]);
    if (!mxIsDouble(prhs[0]) || !(mrows==dim*nn && ncols==1))
        mexErrMsgTxt("yc must be double column vector of length dim*prod(m+1).");
    /* ---------------------------------------------------------------
     *   CASE "V" : Compute volume
     * ---------------------------------------------------------------*/
    if (strcmp("V", flag)==0) {
        plhs[0] = mxCreateDoubleMatrix(nTetra*n, 1, mxREAL);
        result = mxGetPr(plhs[0]);
        for (i=0; i<nTetra*n; i++) {result[i]=0.0;}
        if (dim==2) {
            volume2D(result, yc, m, n);
        } else {
            volume3D(result, yc, m, n);
        }
        /* ---------------------------------------------------------------
         *   CASE "VRange" : Compute range of Vol
         * ---------------------------------------------------------------*/
    } else if (strcmp("VRange", flag)==0) {
        plhs[0] = mxCreateDoubleMatrix(2, 1, mxREAL);
        result = mxGetPr(plhs[0]);
        if (dim==2) {
            volumeRange2D(result, yc, m, n);
        } else {
            volumeRange3D(result, yc, m, n);
        }
        /* ---------------------------------------------------------------
         *   CASE "A" : Compute area
         * ---------------------------------------------------------------*/
    } else if (strcmp("A", flag)==0) {
        plhs[0] = mxCreateDoubleMatrix(nTriangle*n, 1, mxREAL);
        result = mxGetPr(plhs[0]);
        for (i=0; i<nTriangle*n; i++) {result[i]=0.0;}
        area3D(result, yc, m, n);
        /* ---------------------------------------------------------------
         *   CASE "Jac" : Compute Jacobian
         * ---------------------------------------------------------------*/
    } else if (strcmp("Jac", flag)==0) {
        plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
        result  = mxGetPr(plhs[0]);
        for (i=0; i<n; i++) {result[i]=0.0;}
        if (dim==2) {
            Jac2D(result, yc, m, n);
        } else {
            Jac3D(result, yc, m, n);
        }
        /* ---------------------------------------------------------------
         *   CASE "dJacx" : Compute dJac * x
         * ---------------------------------------------------------------*/
    } else if (strcmp("dJacx", flag)==0) {
        if (nrhs < 4)
            mexErrMsgTxt("Vector for dJacx missing.");
        mrows = mxGetM(prhs[3]);
        ncols = mxGetN(prhs[3]);
        if (!mxIsDouble(prhs[3]) || !(mrows==dim*nn && ncols==1))
            mexErrMsgTxt("Inner matrix dimensions must agree.");
        x       =  mxGetPr(prhs[3]);
        plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
        result  = mxGetPr(plhs[0]);
        for (i=0; i<n; i++) {result[i]=0.0;}
        if (dim==2) {
            dJacx2D(result, yc, x, m, n);
        } else {
            dJacx3D(result, yc, x, m, n);
        }
        /* ---------------------------------------------------------------
         *   CASE "dJacadjx" : Compute x'* dJac
         * ---------------------------------------------------------------*/
    } else if (strcmp("dJacadjx", flag)==0) {
        if (nrhs < 3)
            mexErrMsgTxt("Vector for dJacadjx missing.");
        mrows = mxGetM(prhs[3]);
        ncols = mxGetN(prhs[3]);
        if (!mxIsDouble(prhs[3]) || !(mrows== n && ncols==1))
            mexErrMsgTxt("Inner matrix dimensions must agree.");
        x       =  mxGetPr(prhs[3]);
        plhs[0] = mxCreateDoubleMatrix(dim*nn, 1, mxREAL);
        result  = mxGetPr(plhs[0]);
        for (i=0; i<dim*nn; i++) {result[i]=0.0;}
        if (dim==2) {
            dJacTx2D(result, yc, x, m, n);
        } else {
            dJacTx3D(result, yc, x, m, n);
        }
        /* ---------------------------------------------------------------
         *   CASE "S" : Compute S
         * ---------------------------------------------------------------*/
    } else if (strcmp("S", flag)==0) {
        mrows = mxGetM(prhs[4]);
        ncols = mxGetN(prhs[4]);
        if (!mxIsDouble(prhs[4]) || !(mrows==nTetra && ncols==1))
            mexErrMsgTxt("VRef has wrong size.");
        VRef  = mxGetPr(prhs[4]);
        mrows = mxGetM(prhs[5]);
        ncols = mxGetN(prhs[5]);
        if (!mxIsDouble(prhs[5]) || !(mrows==1 && ncols==3))
            mexErrMsgTxt("Alpha has wrong size.");
        alpha = mxGetPr(prhs[5]);
        mrows = mxGetM(prhs[6]);
        ncols = mxGetN(prhs[6]);
        if (!mxIsDouble(prhs[6]) || !(mrows==1 && ncols==dim))
            mexErrMsgTxt("h has wrong size.");
        h     =  mxGetPr(prhs[6]);
        doDerivative = *mxGetLogicals(prhs[7]);
        
        plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
        result = mxGetPr(plhs[0]);
        if (doDerivative) {
            plhs[1] = mxCreateDoubleMatrix(1, dim*nn, mxREAL);
            dSvol = mxGetPr(plhs[1]);
            for (i=0; i<dim*nn; i++) {dSvol[i]=0.0;}
        } else {
            plhs[1] =  mxCreateDoubleMatrix(1, 1, mxREAL);
            dSvol = mxGetPr(plhs[1]);
        }
        if (dim==2) {
            S2D(result, dSvol, yc, VRef,alpha,h, m, n, doDerivative);
        } else {
            mrows = mxGetM(prhs[3]);
            ncols = mxGetN(prhs[3]);
            if (!mxIsDouble(prhs[3]) || !(mrows==nTriangle && ncols==1))
                mexErrMsgTxt("ARef has wrong size.");
            ARef  = mxGetPr(prhs[3]);
            S3D(result, dSvol, yc, ARef,VRef,alpha,h, m, n, doDerivative);
        }
        /* ---------------------------------------------------------------
         *   CASE "d2S" : Compute d2S*x
         * ---------------------------------------------------------------*/
    } else if (strcmp("d2S", flag)==0) {
        mrows = mxGetM(prhs[3]);
        ncols = mxGetN(prhs[3]);
        if (!mxIsDouble(prhs[3]) || !(mrows==dim*nn && ncols==1))
            mexErrMsgTxt("x has wrong size.");
        x     =  mxGetPr(prhs[3]);
        mrows = mxGetM(prhs[5]);
        ncols = mxGetN(prhs[5]);
        if (!mxIsDouble(prhs[5]) || !(mrows==nTetra && ncols==1))
            mexErrMsgTxt("VRef has wrong size.");
        VRef     =  mxGetPr(prhs[5]);
        mrows = mxGetM(prhs[6]);
        ncols = mxGetN(prhs[6]);
        if (!mxIsDouble(prhs[6]) || !(mrows==1 && ncols==3))
            mexErrMsgTxt("Alpha has wrong size.");
        alpha     =  mxGetPr(prhs[6]);
        mrows = mxGetM(prhs[7]);
        ncols = mxGetN(prhs[7]);
        if (!mxIsDouble(prhs[7]) || !(mrows==1 && ncols==dim))
            mexErrMsgTxt("h has wrong size.");
        h     =  mxGetPr(prhs[7]);
        
        plhs[0] = mxCreateDoubleMatrix(dim*nn, 1, mxREAL);
        result = mxGetPr(plhs[0]);
        for (i=0; i<dim*nn; i++) {result[i]=0.0;}
        if (dim==2) {
            d2S2D(result, yc, x, VRef,alpha,h, m, n);
        } else {
            mrows = mxGetM(prhs[4]);
            ncols = mxGetN(prhs[4]);
            if (!mxIsDouble(prhs[4]) || !(mrows==nTriangle && ncols==1))
                mexErrMsgTxt("ARef has wrong size.");
            ARef  = mxGetPr(prhs[4]);
            d2S3D(result, yc, x, ARef, VRef, alpha, h, m, n);
        }
        /* ---------------------------------------------------------------
         *   CASE "d2Sdiag" : Compute diagonal of d2S
         * ---------------------------------------------------------------*/
    } else if (strcmp("d2Sdiag", flag)==0) {
        mrows = mxGetM(prhs[4]);
        ncols = mxGetN(prhs[4]);
        if (!mxIsDouble(prhs[4]) || !(mrows==nTetra && ncols==1))
            mexErrMsgTxt("VRef has wrong size.");
        VRef  = mxGetPr(prhs[4]);
        mrows = mxGetM(prhs[5]);
        ncols = mxGetN(prhs[5]);
        if (!mxIsDouble(prhs[5]) || !(mrows==1 && ncols==3))
            mexErrMsgTxt("Alpha has wrong size.");
        alpha = mxGetPr(prhs[5]);
        mrows = mxGetM(prhs[6]);
        ncols = mxGetN(prhs[6]);
        if (!mxIsDouble(prhs[6]) || !(mrows==1 && ncols==dim))
            mexErrMsgTxt("h has wrong size.");
        h     = mxGetPr(prhs[6]);
        
        plhs[0] = mxCreateDoubleMatrix(dim*nn, 1, mxREAL);
        result = mxGetPr(plhs[0]);
        for (i=0; i<dim*nn; i++) {result[i]=0.0;}
        if (dim==2) {
            d2Sdiag2D(result, yc, VRef,alpha,h, m, n);
        } else {
            mrows = mxGetM(prhs[3]);
            ncols = mxGetN(prhs[3]);
            if (!mxIsDouble(prhs[3]) || !(mrows==nTriangle && ncols==1))
                mexErrMsgTxt("ARef has wrong size.");
            ARef  = mxGetPr(prhs[3]);
            d2S3Ddiag(result, yc, ARef, VRef, alpha, h, m, n);
        }
        /* ---------------------------------------------------------------
         *   CASE "VAMPIREdiag" : Compute diag of d2J for VAMPIRE
         * ---------------------------------------------------------------*/
    } else if (strcmp("VAMPIREdiag", flag)==0) {
        //
        // Usage: diag = geometrymexC(yc,m,'VAMPIREdiag',Jac,dT,T,h)
        //
        // get Jac
        mrows = mxGetM(prhs[3]);
        ncols = mxGetN(prhs[3]);
        if (!mxIsDouble(prhs[3]) || !(mrows==n && ncols==1))
            mexErrMsgTxt("Jac has wrong size.");
        Jac   = mxGetPr(prhs[3]);
        // get dT
        mrows = mxGetM(prhs[4]);
        ncols = mxGetN(prhs[4]);
        if (!mxIsDouble(prhs[4]) || !(mrows==n*dim && ncols==1))
            mexErrMsgTxt("dT has wrong size.");
        dT    = mxGetPr(prhs[4]);
        // get T
        mrows = mxGetM(prhs[5]);
        ncols = mxGetN(prhs[5]);
        if (!mxIsDouble(prhs[5]) || !(mrows==n && ncols==1))
            mexErrMsgTxt("T has wrong size.");
        T     = mxGetPr(prhs[5]);
        
        mrows = mxGetM(prhs[6]);
        ncols = mxGetN(prhs[6]);
        if (!mxIsDouble(prhs[6]) || !(mrows==1 && ncols==dim))
            mexErrMsgTxt("h has wrong size.");
        h     = mxGetPr(prhs[6]);
        
        plhs[0] = mxCreateDoubleMatrix(dim*nn, 1, mxREAL);
        result  = mxGetPr(plhs[0]);
        
        for (i=0; i<dim*nn; i++) {result[i]=0.0;}
        if (dim==2) {
            VAMPIREdiag2D(result, yc,m,Jac,dT,T,h);
        } else {
            VAMPIREdiag3D(result, yc,m,Jac,dT,T,h);
        }
    }
    else {
        mexErrMsgTxt("Invalid Flag");
    }
    mxFree(flag);
}
/*
 *   =======================================================================================
 *   FAIR: Flexible Algorithms for Image Registration, Version 2011
 *   Copyright (c): Jan Modersitzki
 *   Maria-Goeppert-Str. 1a, D-23562 Luebeck, Germany
 *   Email: jan.modersitzki@mic.uni-luebeck.de
 *   URL:   http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
 *   =======================================================================================
 *   No part of this code may be reproduced, stored in a retrieval system,
 *   translated, transcribed, transmitted, or distributed in any form
 *   or by any means, means, manual, electric, electronic, electro-magnetic,
 *   mechanical, chemical, optical, photocopying, recording, or otherwise,
 *   without the prior explicit written permission of the authors or their
 *   designated proxies. In no event shall the above copyright notice be
 *   removed or altered in any way.
 *
 *   This code is provided "as is", without any warranty of any kind, either
 *   expressed or implied, including but not limited to, any implied warranty
 *   of merchantibility or fitness for any purpose. In no event will any party
 *   who distributed the code be liable for damages or for any claim(s) by
 *   any other party, including but not limited to, any lost profits, lost
 *   monies, lost data or data rendered inaccurate, losses sustained by
 *   third parties, or any other special, incidental or consequential damages
 *   arrising out of the use or inability to use the program, even if the
 *   possibility of such damages has been advised against. The entire risk
 *   as to the quality, the performace, and the fitness of the program for any
 *   particular purpose lies with the party using the code.
 *   =======================================================================================
 *   Any use of this code constitutes acceptance of the terms of the above statements
 *   =======================================================================================
 */