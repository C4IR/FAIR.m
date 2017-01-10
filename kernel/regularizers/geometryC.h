/*
 *=======================================================================================
 * (c) Lars Ruthotto 2011/01/11, see FAIR.2 and FAIRcopyright.m.
 *
 * contains all relevant functions for matrix free hyperelastic regularization as
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
// #pragma once

#include "mex.h"
#include <stdio.h>
#include <string.h>
#include "matrix.h"

#ifdef _OPENMP
#include<omp.h>
#endif

// 2D
void volume2D(double *volu, const double *yc, const double *m,const int n);
void volumeRange2D(double *range, const double *yc, const double *m,const int n);
void Jac2D(double *volu, const double *yc,const double *m,const int n);
double dVxTriangle2D(double *y1, double *y2, double *yM, double *x1, double *x2, double *xM);
void dJacx2D(double *result, const double *yc,const double *x,const double *m, int n);
// inline void dVTTriangle2D(double *y1, double *y2, double *yM, double x, double *r1, double *r2, double *rM);
void dJacTx2D(double *prod, const double *yc,const double *x,const double *m, int n);
// inline double SvolTriangle2D(double *y1, double *y2, double *yM,const double VRef, bool doDerivative, double *r1, double *r2, double *rM,double alphaVolume=1.0);

void S2D(double *S, double *dS, const double *yc,const double *VRef, const double *alpha ,const double *h,const double *m,const int n, const bool doDerivative);
void d2SvolTriangle2D(double *y1, double *y2, double *yM, double *x1, double *x2, double *xM, const double VRef, double *r1, double *r2, double *rM,double alphaVolume);
void d2S2D(double *prod, const double *yc,const double *x,const double *VRef,const double *alpha,const double *h,const double *m, const int n);
void d2SdiagTriangle(double *y1, double *y2, double *yM, const double VRef, double *r1, double *r2, double *r3, double *r4, double alphaVolume=1.0);
void d2Sdiag2D(double *prod, const double *yc,const double *VRef,const double *alpha,const double *h,const double *m, const int n);


// 3D
void volume3D(double *volu, const double *yc, const double *m, const int n);
void volumeRange3D(double *range, const double *yc, const double *m,const int n);
void Jac3D(double *volu, const double *yc,const double *m, int n);

double dVxTetra3D(double *cofDy, const double *y1, const double *y2,const double *y3, const double *yM, const double *x1,const double *x2, const double *x3,const double *xM);
void dJacx3D(double *result, const double *yc,const double *x,const double *m,const int n);
void dVTTetra3D(const double *yA,const double *yB,const double *yC,const double *yM,const double x, double *rA, double *rB, double *rC, double *rM);
void dJacTx3D(double *prod, const double *yc,const double *x,const double *m, int n);

double arTriangle3D(double *y1, double *y2, double *yM,const double ARef, double *r1, double *r2, double *rM, bool doDerivative,const double alphaArea);
double areaTriangle3D(double *y1, double *y2, double *yM);

double SvolTetra3D(double *yA, double *yB, double *yC, double *yM, const double VRef, bool doDerivative, double *rA, double *rB, double *rC, double *rM, double inv,const double alphaVolume);
// inline void d2SvolTetra3D(double *yA, double *yB, double *yC, double *yM, double *xA, double *xB, double *xC, double *xM,const double VRef, double *rA, double *rB, double *rC, double *rM, const double inv,const double alphaVolume);
// inline void d2SarTriangle3D( double *y1,double *y2, double *yM, double *x1, double *x2, double *xM,const double ARef, double *r1, double *r2, double *rM,const double alphaArea);
void S3D(double *S, double *dS, const double *yc,const double *ARef,const double *VRef, const double *alpha,const double *h,const double *m, int n, bool doDerivative);
void d2S3D(double *prod, const double *yc,const double *x,const double *ARef,const double *VRef,const  double *alpha,const double *h,const double *m, int n);
// inline void d2SdiagTetra(double *y1, double *y2, double *yF, double *yM,const double ARef, const double VRef, double *r1, double *r2, double *r3, double *r4, double *r5, double *r6, double *r7, double *r8,const double inv , double alphaArea, double alphaVolume);
void d2S3Ddiag(double *prod, const double *yc,const double *ARef,const double *VRef,const  double *alpha,const double *h,const double *m, int n);
void area3D(double *A,const double *yc, const double *m, int n);


// VAMPIRE
void VAMPIREdiagTriangle(double *y1, double *y2, double *yM, double *r1, double *r2, double *r3, double *r4, double Ti);
void VAMPIREdiag2D(double *diag, const double *yc, const double *m, const double *Jac, const double *dT, const double *Tc,const double *h);
// inline void VAMPIREdiagTetra(double *y1, double *y2, double *yF, double *yM, double *r1, double *r2, double *r3, double *r4, double *r5, double *r6, double *r7, double *r8,double Ti);
void VAMPIREdiag3D(double *diag, const double *yc, const double *m, const double *Jac, const double *dT, const double *Tc,const double *h);


// EPI
void EPIdiagTriangle(double *y1, double *y2, double *yM, double *r1, double *r2, double *r3, double *r4, const double *v, double Ti);
void EPIdiag2D(double *diag, const double *yc1, const double *yc2, const double *m, const double *Jac1, const double *Jac2, const double *dI1,const double *dI2, const double *I1, const double *I2, const double *v1, const double *v2,const double *h);
// inline void EPIdiagTetra(double *y1, double *y2, double *yF, double *yM, double *r1, double *r2, double *r3, double *r4, double *r5, double *r6, double *r7, double *r8, const double *v,double Ti);
void EPIdiag3D(double *diag, const double *yc1, const double *yc2, const double *m, const double *Jac1, const double *Jac2, const double *dI1,const double *dI2, const double *I1, const double *I2, const double *v1, const double *v2,const double *h);

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