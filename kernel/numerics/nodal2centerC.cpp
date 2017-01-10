/*
 * (c) Fabian Gigengack 2012/05/15, see FAIRcopyright.m.
 * http://www.uni-muenster.de/EIMI/
 *
 * CPP code for the conversion of nodal to centered grids - or vice versa.
 */

#include <math.h>
#include <mex.h>
#ifdef _OPENMP
#include <omp.h>
#endif


void center2nodal(double *yOut, const double *yc, const int *m, const int dim);
void nodal2center(double *yOut, const double *yc, const int *m, const int dim);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    int i, dimLocal, n = 1, nn = 1, mLocal[3];
    if (nrhs<4)
        mexErrMsgTxt("Number of arguments must be 4!");
    
    // Get input
    const double *yc  = static_cast<double*>(mxGetData(prhs[0]));
    const double *m   = static_cast<double*>(mxGetData(prhs[1]));
    const double *dim = static_cast<double*>(mxGetData(prhs[2]));
    bool mode = (bool) *mxGetLogicals(prhs[3]);
    
    dimLocal = (int) dim[0];
    for (i=0; i<dimLocal; i++) {
        mLocal[i] = (int) m[i];
        n  *= mLocal[i];
        nn *= mLocal[i]+1;
    }
    
    mwSize dims[2];
    
    if (mode) {
        dims[0] = nn*dimLocal; dims[1] = 1;
        plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
        double *yOut = static_cast<double*>(mxGetData(plhs[0]));
        
        center2nodal(yOut, yc, mLocal, dimLocal);
    } else {
        dims[0] = n*dimLocal; dims[1] = 1;
        plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
        double *yOut = static_cast<double*>(mxGetData(plhs[0]));
        
        nodal2center(yOut, yc, mLocal, dimLocal);
    }
}

void center2nodal(double *yOut, const double *yc, const int *m, const int dim) {
    int idx,idxn,x,y,z,d;
    int offx,offy,offz;
    int m0 = m[0], m1 = m[1];
    int incx=1, incy=m0, incz, incd, incyn=m0+1, inczn, incdn;
    double tmp;
   
    if (dim==2) {
        incdn = (m0+1)*(m1+1);
        incd  = (m0)*(m1);
        for (d=0; d<dim; d++) {
            for (offx=0; offx<2; offx++) {
                for (offy=0; offy<2; offy++) {
                    #pragma omp parallel for private(idx, idxn,y,x,tmp)
                    for (y=offy; y<m1; y+=2) {
                        for (x=offx; x<m0; x+=2) {
                            idx  = x+y*incx +d*incd;
                            idxn = x+y*incx +d*incdn;
                            
                            tmp = 0.25 * yc[idx];
                            
                            yOut[idxn               ]  += tmp;
                            yOut[idxn + incx        ]  += tmp;
                            yOut[idxn        + incyn]  += tmp;
                            yOut[idxn + incx + incyn]  += tmp;
                        }
                    }
                    
                }
            }
        }
    } else {
        int m2 = m[2];
        
        inczn = (m0+1)*(m1+1);
        incz  = (m0)*(m1);
        incdn = (m0+1)*(m1+1)*(m2+1);
        incd  = (m0)*(m1)*m2;
        for (d=0; d<dim; d++) {
            for (offx=0; offx<2; offx++) {
                for (offy=0; offy<2; offy++) {
                    for (offz=0; offz<2; offz++) {
                        #pragma omp parallel for default(shared) private(idx, idxn,z,y,x,tmp)
                        for (z=offz; z<m2; z+=2) {
                            for (y=offy; y<m1; y+=2) {
                                for (x=offx; x<m0; x+=2) {
                                    idx  = x + y*incy  + z*incz  + d*incd;
                                    idxn = x + y*incyn + z*inczn + d*incdn;
                                    
                                    tmp  = 0.125 * yc[idx];
                                    
                                    // add to all eight surrounding nodes
                                    yOut[idxn                       ] += tmp;
                                    yOut[idxn + incx                ] += tmp;
                                    yOut[idxn + incx + incyn        ] += tmp;
                                    yOut[idxn        + incyn        ] += tmp;
                                    
                                    yOut[idxn                + inczn] += tmp;
                                    yOut[idxn + incx         + inczn] += tmp;
                                    yOut[idxn + incx + incyn + inczn] += tmp;
                                    yOut[idxn        + incyn + inczn] += tmp;
                                    
                                    
                                }
                            }
                        }
                    }
                }
            }
            
        }
    }
}


void nodal2center(double *yOut, const double *yc, const int *m, const int dim) {
    int idx,idxn,x,y,z,d;
    int m0 = m[0], m1 = m[1];
    
    if (dim==2) {
        #pragma omp parallel for private(idx, idxn,y,x)
        for (d=0; d<dim; d++) {
            for (y=0; y<m1; y++) {
                for (x=0; x<m0; x++) {
                    idx  = x+y*m0+d*m0*m1;
                    idxn = x+y*(m0+1)+d*(m0+1)*(m1+1);
                    yOut[idx] = 0.25 * (yc[idxn]
                                        + yc[idxn+1]
                                        + yc[idxn+(m0+1)]
                                        + yc[idxn+(m0+1)+1]);
                }
            }
        }
    } else {
        int m2 = m[2];
        //         #pragma omp parallel for private(idx, idxn)
        #pragma omp parallel for private(idx, idxn,z,y,x)
          for (d=0; d<dim; d++) {
            for (z=0; z<m2; z++) {
                for (y=0; y<m1; y++) {
                    for (x=0; x<m0; x++) {
                        idx  = x+y*m0+z*m0*m1+d*m0*m1*m2;
                        idxn = x+y*(m0+1)+z*(m0+1)*(m1+1)+d*(m0+1)*(m1+1)*(m2+1);
                        yOut[idx] = 0.125 * (yc[idxn]
                                             + yc[idxn+1]
                                             + yc[idxn+(m0+1)]
                                             + yc[idxn+(m0+1)+1]
                                             + yc[idxn+(m0+1)*(m1+1)]
                                             + yc[idxn+(m0+1)*(m1+1)+1]
                                             + yc[idxn+(m0+1)*(m1+1)+(m0+1)]
                                             + yc[idxn+(m0+1)*(m1+1)+(m0+1)+1]);
                    }
                }
            }
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