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

// #include "geometryC.h"

#pragma mark -
#pragma mark helper
inline void vecminus(double *res,const double *yA,const double *yB,const int dim) {
    int i;
    for (i=0; i<dim; i++) {
        res[i] = yA[i] - yB[i];
    }
}

#pragma mark -
#pragma mark determinant, cofactor and penalty functions
inline double phicaley(const double x){
    // phi(x) = (x-1)^2 /x
    // return (x-1.0)*(x-1.0) /x;
    return (x-1.0)*(x-1.0)*(x-1.0)*(x-1.0) /(x*x);
}
inline double dphicaley(const double x){
    // dphi(x) = 1- 1/(x^2)
    // return 1.0- 1.0/(x*x);
    return 2*(x-1)*(x-1)*(x-1)*(x+1)/(x*x*x);
}
inline double d2phicaley(const double x){
    // d2phi(x) = 2 / (x^3)
    // return 2.0 / (x*x*x);
    double x4 = (x*x*x*x);
    return 2*(x4 - 4*x + 3)/x4;
}

inline double phiarea(const double x){
    // phi(x) = .5 (x-1)^2
    return .5*(x-1.0)*(x-1.0);
//     return (x-1.0)*(x-1.0) /x;
}
inline double dphiarea(const double x){
    // dphi(x) = (x-1)
    return (x-1.0);
//     return 1.0- 1.0/(x*x);
}
inline double d2phiarea(const double x){
    // d2phi(x) = 1
    return  1.0;
//     return 2.0 / (x*x*x);
}

//inline double phiarea(const double x){
//    // convex penalty (needed for existence)
//    //
//    //          ( .5 (x-1)^2, if x>1.0
//    // phi(x) = ( 0.0       , if x<1.0
//    if(x<1.0) {
//        return 0.0;
//    }
//    else {
//        return .5*(x-1.0)*(x-1.0);
//    }
//    //     return (x-1.0)*(x-1.0) /x;
//}
//inline double dphiarea(const double x){
//    // convex penalty (needed for existence)
//    //
//    //          ( x-1       , if x>1.0
//    // dphi(x) = ( 0.0       , if x<1.0
//    if(x<1.0) {
//        return 0.0;
//    }
//    else {
//        return (x-1.0);
//    }
//}
//inline double d2phiarea(const double x){
//    // convex penalty (needed for existence)
//    //
//    //            ( 1.0       , if x>1.0
//    // d2phi(x) = ( 0.0       , if x<1.0
//   if(x<1.0) {
//        return 0.0;
//    }
//        else {
//            return 1.0;
//    }
//}


inline void cofactor2D(double *cofDy,const double *y1,const double *y2,const double *yM, double x){
    x = x / 2.0;
    cofDy[0] =  x*(y2[1] - yM[1]); /* a_11 */
    cofDy[1] = -x*(y1[1] - yM[1]); /* a_12 */
    cofDy[2] = -x*(y2[0] - yM[0]); /* a_21 */
    cofDy[3] =  x*(y1[0] - yM[0]); /* a_22 */
}

inline double det3D(const double *col1, const double *col2, const double *col3, const double *colM) {
    return ((col1[0]-colM[0]) * (col2[1]-colM[1]) * (col3[2]-colM[2])
    + (col1[2]-colM[2]) * (col2[0]-colM[0]) * (col3[1]-colM[1])
    + (col1[1]-colM[1]) * (col2[2]-colM[2]) * (col3[0]-colM[0])
    - (col1[2]-colM[2]) * (col2[1]-colM[1]) * (col3[0]-colM[0])
    - (col1[0]-colM[0]) * (col2[2]-colM[2]) * (col3[1]-colM[1])
    - (col1[1]-colM[1]) * (col2[0]-colM[0]) * (col3[2]-colM[2])) / 6.0;
}

inline void cofactor3D(double *cofDy, const double* p1, const double *p2,const double *p3,const double *pM,double w){
    w = w/6.0;
    cofDy[0] = w*((p2[1]-pM[1]) * (p3[2]-pM[2]) - (p2[2]-pM[2]) * (p3[1]-pM[1])); // det(A_11)
    cofDy[1] = w*((p3[1]-pM[1]) * (p1[2]-pM[2]) - (p3[2]-pM[2]) * (p1[1]-pM[1])); //-det(A_12)
    cofDy[2] = w*((p1[1]-pM[1]) * (p2[2]-pM[2]) - (p1[2]-pM[2]) * (p2[1]-pM[1])); // det(A_13)
    cofDy[3] = w*((p3[0]-pM[0]) * (p2[2]-pM[2]) - (p3[2]-pM[2]) * (p2[0]-pM[0])); // -det(A_21)
    cofDy[4] = w*((p1[0]-pM[0]) * (p3[2]-pM[2]) - (p1[2]-pM[2]) * (p3[0]-pM[0])); //  det(A_22)
    cofDy[5] = w*((p2[0]-pM[0]) * (p1[2]-pM[2]) - (p2[2]-pM[2]) * (p1[0]-pM[0])); // -det(A_23)
    cofDy[6] = w*((p2[0]-pM[0]) * (p3[1]-pM[1]) - (p2[1]-pM[1]) * (p3[0]-pM[0])); // -det(A_21)
    cofDy[7] = w*((p3[0]-pM[0]) * (p1[1]-pM[1]) - (p3[1]-pM[1]) * (p1[0]-pM[0])); //  det(A_22)
    cofDy[8] = w*((p1[0]-pM[0]) * (p2[1]-pM[1]) - (p1[1]-pM[1]) * (p2[0]-pM[0])); // -det(A_23)
    
}

#pragma mark -
#pragma mark mfVolume

inline double volTriangle2D(double *y1, double *y2, double *yM){
    return  ((y1[0]-yM[0]) * (y2[1]-yM[1]) - (y1[1]-yM[1]) * (y2[0]-yM[0]))/2.0;
}

void volume2D( double *volu, const double *yc, const double *m,const int n){
    double yA[2], yB[2], yC[2], yD[2], yM[2];
    double vol[4];
    int i1, i2, r,id, idn,comp;
    
    const int m1 = m[0], m2 = m[1];
    const int nn = (m1+1)*(m2+1);
    // increments
    const int inc_x=1;
    const int inc_y= m1+1;
    //pointers
    const double *ptr;
#pragma omp parallel for  default(shared) private(ptr,id,idn, i1, i2,r, vol,comp, yA, yB, yC, yD, yM)
        for (i2=0; i2<m2; i2++){
            for (i1=0; i1<m1; i1++){
                // get edges and cell center (shared --> private)
                idn = i1+i2*inc_y;
                comp = 0;
                for (r=0; r<2; r++) {
                    ptr = yc + idn+comp;
                    yA[r] = *ptr; // yc[idn           +comp];
                    ptr = ptr+inc_x;
                    yB[r] = *ptr; //yc[idn+1         +comp];
                    ptr = ptr+inc_y;
                    yD[r] = *ptr; //yc[idn+1+(m1+1)  +comp];
                    ptr = ptr-inc_x;
                    yC[r] = *ptr; //* yc[idn  +(m1+1)  +comp];
                    yM[r] = 0.25*(yA[r]+yB[r]+yC[r]+yD[r]);
                    // init private result vectors (private)
                    comp+=nn;
                }
                
                // compute volumes of triangles (private)
                vol[0] = volTriangle2D(yA, yB, yM);
                vol[1] = volTriangle2D(yB, yD, yM);
                vol[2] = volTriangle2D(yD, yC, yM);
                vol[3] = volTriangle2D(yC, yA, yM);
                
                // write into result vector (private --> shared)
                id = i1+i2*m1;
                for(r=0;r<4;r++){
                    volu[id] = vol[r];
                    id+=n;
                }
                
            }
        }
        return;
}

void volumeRange2D( double *range, const double *yc, const double *m,const int n){
    double yA[2], yB[2], yC[2], yD[2], yM[2];
    double vol, mini=10000.0, maxi=0.0;
    int i1, i2, r,id, idn,comp;
    
    const int m1 = m[0], m2 = m[1];
    const int nn = (m1+1)*(m2+1);
    // increments
    const int inc_x=1;
    const int inc_y= m1+1;
    //pointers
    const double *ptr;
    for (i2=0; i2<m2; i2++){
        for (i1=0; i1<m1; i1++){
            // get edges and cell center (shared --> private)
            idn = i1+i2*inc_y;
            comp = 0;
            for (r=0; r<2; r++) {
                ptr = yc + idn+comp;
                yA[r] = *ptr; // yc[idn           +comp];
                ptr = ptr+inc_x;
                yB[r] = *ptr; //yc[idn+1         +comp];
                ptr = ptr+inc_y;
                yD[r] = *ptr; //yc[idn+1+(m1+1)  +comp];
                ptr = ptr-inc_x;
                yC[r] = *ptr; //* yc[idn  +(m1+1)  +comp];
                yM[r] = 0.25*(yA[r]+yB[r]+yC[r]+yD[r]);
                // init private result vectors (private)
                comp+=nn;
            }
            
            // compute volumes of triangles (private)
            vol = volTriangle2D(yA, yB, yM);
            if (mini>vol) {    mini=vol;} if (maxi<vol) {    maxi=vol;}
            vol = volTriangle2D(yB, yD, yM);
            if (mini>vol) {    mini=vol;}if (maxi<vol) {    maxi=vol;}
            vol = volTriangle2D(yD, yC, yM);
            if (mini>vol) {    mini=vol;}if (maxi<vol) {    maxi=vol;}
            vol = volTriangle2D(yC, yA, yM);
            if (mini>vol) {    mini=vol;}if (maxi<vol) {    maxi=vol;}
        }
    }
    *range = mini; range++;
    *range = maxi;
    return;
}

void Jac2D( double *volu, const double *yc,const double *m,const int n){
    double yA[2], yB[2], yC[2], yD[2], yM[2], vol;
    int i1, i2, r,id, idn,comp;
    
    const int m1 = m[0], m2 = m[1];
    const int nn = (m1+1)*(m2+1);
    // increments
    const int inc_x=1;
    const int inc_y= m1+1;
    //pointers
    const double *ptr;
#pragma omp parallel for  default(shared) private(ptr,i1, i2, r,id, idn,comp, vol, yA, yB, yC, yD, yM)
        for (i2=0; i2<m2; i2++){
            for (i1=0; i1<m1; i1++){
                // get edges and cell center (shared --> private)
                comp = 0;
                idn = i1+i2*inc_y;
                for (r=0; r<2; r++) {
                ptr = yc + idn+comp;
                yA[r] = *ptr; // yc[idn           +comp];
                ptr = ptr+inc_x;
                yB[r] = *ptr; //yc[idn+1         +comp];
                ptr = ptr+inc_y;
                yD[r] = *ptr; //yc[idn+1+(m1+1)  +comp];
                ptr = ptr-inc_x;
                yC[r] = *ptr; //* yc[idn  +(m1+1)  +comp];
                yM[r] = 0.25*(yA[r]+yB[r]+yC[r]+yD[r]);
                // init private result vectors (private)
                comp+=nn;
            }
            
            // compute volumes of triangles and sum up to pixel volume (private)
            vol  = volTriangle2D(yA, yB, yM);
            vol += volTriangle2D(yB, yD, yM);
            vol += volTriangle2D(yD, yC, yM);
            vol += volTriangle2D(yC, yA, yM);
            
            // write into result vector (private --> shared)
            id = i1+i2*m1;
            volu[id] = vol;
        }
    }
    return;
}

void volume3D(double *volu, const double *yc, const double *m, const int n) {
    const int m1 = m[0], m2 = m[1], m3 = m[2];
    double yA[3], yB[3], yC[3], yD[3], yE[3], yF[3], yG[3], yH[3];
    double yABDC[3], yBDHF[3], yABFE[3], yCAEG[3], yDCGH[3], yEFHG[3], yM[3];
    double vol[24];
    
    int id, i1, i2, i3, r, comp, idn;
    
    const int nn = (m1+1)*(m2+1)*(m3+1);
    const int incn_x = 1;
    const int incn_y = m1+1;
    const int incn_z = (m1+1) * (m2+1);
    const int inc_x = 1;
    const int inc_y = m1;
    const int inc_z = m1 * m2;
    
    const double *rptr;
    
for (i1=0; i1<m1; i1++) {
    for (i2=0; i2<m2; i2++) {
        for (i3=0; i3<m3; i3++) {
            id  = i1 + i2 * inc_y  + i3 * inc_z;
            idn = i1 + i2 * incn_y + i3 * incn_z;
            // get edges, face stg and cell center of yc (shared --> private)
            comp = 0;
            for (r=0; r<3; r++) {
                rptr = yc+idn+comp;
                // A
                yA[r] = *rptr;
                // B
                rptr += incn_x;
                yB[r] = *rptr;
                // D
                rptr += incn_y;
                yD[r] = *rptr;
                // C
                rptr -= incn_x;
                yC[r] = *rptr;
                // G
                rptr += incn_z;
                yG[r] = *rptr;
                // H
                rptr += incn_x;
                yH[r] = *rptr;
                // F
                rptr -= incn_y;
                yF[r] = *rptr;
                // E
                rptr -= incn_x;
                yE[r] = *rptr;
                
                yABDC[r] = 0.25*(yA[r]+yB[r]+yD[r]+yC[r]);
                yBDHF[r] = 0.25*(yB[r]+yD[r]+yH[r]+yF[r]);
                yABFE[r] = 0.25*(yA[r]+yB[r]+yF[r]+yE[r]);
                yCAEG[r] = 0.25*(yC[r]+yA[r]+yE[r]+yG[r]);
                yDCGH[r] = 0.25*(yD[r]+yC[r]+yG[r]+yH[r]);
                yEFHG[r] = 0.25*(yE[r]+yF[r]+yH[r]+yG[r]);
                yM[r]    = 0.5 *(yABDC[r]+yEFHG[r]);
                
                comp +=nn;
            }
            
            // compute volumes of tetrahedra (private)
            vol[0]  = -det3D(yA, yB, yABDC, yM);
            vol[1]  = -det3D(yB, yD, yABDC, yM);
            vol[2]  = -det3D(yD, yC, yABDC, yM);
            vol[3]  = -det3D(yC, yA, yABDC, yM);
            
            vol[4]  =  det3D(yB, yD, yBDHF, yM);
            vol[5]  =  det3D(yD, yH, yBDHF, yM);
            vol[6]  =  det3D(yH, yF, yBDHF, yM);
            vol[7]  =  det3D(yF, yB, yBDHF, yM);
            
            vol[8]  =  det3D(yA, yB, yABFE, yM);
            vol[9]  =  det3D(yB, yF, yABFE, yM);
            vol[10] =  det3D(yF, yE, yABFE, yM);
            vol[11] =  det3D(yE, yA, yABFE, yM);
            
            vol[12] =  det3D(yC, yA, yCAEG, yM);
            vol[13] =  det3D(yA, yE, yCAEG, yM);
            vol[14] =  det3D(yE, yG, yCAEG, yM);
            vol[15] =  det3D(yG, yC, yCAEG, yM);
            
            vol[16] =  det3D(yD, yC, yDCGH, yM);
            vol[17] =  det3D(yC, yG, yDCGH, yM);
            vol[18] =  det3D(yG, yH, yDCGH, yM);
            vol[19] =  det3D(yH, yD, yDCGH, yM);
            
            vol[20] =  det3D(yE, yF, yEFHG, yM);
            vol[21] =  det3D(yF, yH, yEFHG, yM);
            vol[22] =  det3D(yH, yG, yEFHG, yM);
            vol[23] =  det3D(yG, yE, yEFHG, yM);
            
            // write into result vector (private --> shared)
            for (r=0; r<24; r++) {
                volu[id] = vol[r];
                id += n;
            }
        }
    }
}
return;
}

void volumeRange3D( double *range, const double *yc, const double *m,const int n){
    const int m1 = m[0], m2 = m[1],m3 = m[2];
    double yA[3], yB[3], yC[3], yD[3], yE[3], yF[3], yG[3], yH[3], yABDC[3], yBDHF[3], yABFE[3], yCAEG[3], yDCGH[3], yEFHG[3], yM[3];
    double vol,mini=1000.0, maxi=0.0;
    
    int id, i1, i2, i3, r, comp,idn;
    
    const int nn = (m1+1)*(m2+1)*(m3+1);
    const int inc_x = 1;
    const int inc_y = m1+1;
    const int inc_z = (m1+1)*(m2+1);
    
    const double *rptr;
    for (i3=0; i3<m3; i3++){
        for (i2=0; i2<m2; i2++){
            for (i1=0; i1<m1; i1++){
                idn = i1+ i2* inc_y+ i3*inc_z;
                // get edges, face stg and cell center of yc (shared --> private)
                comp = 0;
                for(r=0;r<3;r++){
                    rptr = yc+idn+comp;
                    // A
                    yA[r] = *rptr;
                    // B
                    rptr += inc_x;
                    yB[r] = *rptr;
                    // D
                    rptr += inc_y;
                    yD[r] = *rptr;
                    // C
                    rptr -= inc_x;
                    yC[r] = *rptr;
                    // G
                    rptr += inc_z;
                    yG[r] = *rptr;
                    // H
                    rptr += inc_x;
                    yH[r] = *rptr;
                    // F
                    rptr -= inc_y;
                    yF[r] = *rptr;
                    // E
                    rptr -= inc_x;
                    yE[r] = *rptr;
                    
                    yABDC[r] = 0.25*(yA[r]+yB[r]+yD[r]+yC[r]);
                    yBDHF[r] = 0.25*(yB[r]+yD[r]+yH[r]+yF[r]);
                    yABFE[r] = 0.25*(yA[r]+yB[r]+yF[r]+yE[r]);
                    yCAEG[r] = 0.25*(yC[r]+yA[r]+yE[r]+yG[r]);
                    yDCGH[r] = 0.25*(yD[r]+yC[r]+yG[r]+yH[r]);
                    yEFHG[r] = 0.25*(yE[r]+yF[r]+yH[r]+yG[r]);
                    yM[r]    = 0.5 *(yABDC[r]+yEFHG[r]);
                    
                    comp +=nn;
                }
                
                
                // compute volumes of tetrahedra (private)
                vol = -det3D(yA, yB, yABDC, yM);
                if (mini>vol) {    mini=vol;} if (maxi<vol) {    maxi=vol;}
                vol = -det3D(yB, yD, yABDC, yM);
                if (mini>vol) {    mini=vol;} if (maxi<vol) {    maxi=vol;}
                vol = -det3D(yD, yC, yABDC, yM);
                if (mini>vol) {    mini=vol;} if (maxi<vol) {    maxi=vol;}
                vol = -det3D(yC, yA, yABDC, yM);
                if (mini>vol) {    mini=vol;} if (maxi<vol) {    maxi=vol;}
                
                vol =  det3D(yB, yD, yBDHF, yM);
                if (mini>vol) {    mini=vol;} if (maxi<vol) {    maxi=vol;}
                vol =  det3D(yD, yH, yBDHF, yM);
                if (mini>vol) {    mini=vol;} if (maxi<vol) {    maxi=vol;}
                vol =  det3D(yH, yF, yBDHF, yM);
                if (mini>vol) {    mini=vol;} if (maxi<vol) {    maxi=vol;}
                vol =  det3D(yF, yB, yBDHF, yM);
                if (mini>vol) {    mini=vol;} if (maxi<vol) {    maxi=vol;}
                
                vol =  det3D(yA, yB, yABFE, yM);
                if (mini>vol) {    mini=vol;} if (maxi<vol) {    maxi=vol;}
                vol =  det3D(yB, yF, yABFE, yM);
                if (mini>vol) {    mini=vol;} if (maxi<vol) {    maxi=vol;}
                vol =  det3D(yF, yE, yABFE, yM);
                if (mini>vol) {    mini=vol;} if (maxi<vol) {    maxi=vol;}
                vol =  det3D(yE, yA, yABFE, yM);
                if (mini>vol) {    mini=vol;} if (maxi<vol) {    maxi=vol;}
                
                vol =  det3D(yC, yA, yCAEG, yM);
                if (mini>vol) {    mini=vol;} if (maxi<vol) {    maxi=vol;}
                vol =  det3D(yA, yE, yCAEG, yM);
                if (mini>vol) {    mini=vol;} if (maxi<vol) {    maxi=vol;}
                vol =  det3D(yE, yG, yCAEG, yM);
                if (mini>vol) {    mini=vol;} if (maxi<vol) {    maxi=vol;}
                vol =  det3D(yG, yC, yCAEG, yM);
                if (mini>vol) {    mini=vol;} if (maxi<vol) {    maxi=vol;}
                
                vol =  det3D(yD, yC, yDCGH, yM);
                if (mini>vol) {    mini=vol;} if (maxi<vol) {    maxi=vol;}
                vol =  det3D(yC, yG, yDCGH, yM);
                if (mini>vol) {    mini=vol;} if (maxi<vol) {    maxi=vol;}
                vol =  det3D(yG, yH, yDCGH, yM);
                if (mini>vol) {    mini=vol;} if (maxi<vol) {    maxi=vol;}
                vol =  det3D(yH, yD, yDCGH, yM);
                if (mini>vol) {    mini=vol;} if (maxi<vol) {    maxi=vol;}
                
                vol =  det3D(yE, yF, yEFHG, yM);
                if (mini>vol) {    mini=vol;} if (maxi<vol) {    maxi=vol;}
                vol =  det3D(yF, yH, yEFHG, yM);
                if (mini>vol) {    mini=vol;} if (maxi<vol) {    maxi=vol;}
                vol =  det3D(yH, yG, yEFHG, yM);
                if (mini>vol) {    mini=vol;} if (maxi<vol) {    maxi=vol;}
                vol =  det3D(yG, yE, yEFHG, yM);
                if (mini>vol) {    mini=vol;} if (maxi<vol) {    maxi=vol;}
            }
        }
    }
    *range = mini; range++;
    *range = maxi;
    return;
}

void Jac3D( double *volu, const double *yc,const double *m, int n){
    const int m1 = m[0], m2 = m[1],m3 = m[2];
    double yA[3], yB[3], yC[3], yD[3], yE[3], yF[3], yG[3], yH[3], yABDC[3], yBDHF[3], yABFE[3], yCAEG[3], yDCGH[3], yEFHG[3], yM[3];
    double vol;
    
    int id, i1, i2, i3, r, comp,idn;
    
    const int nn = (m1+1)*(m2+1)*(m3+1);
    const int inc_x = 1;
    const int inc_y = m1+1;
    const int inc_z = (m1+1)*(m2+1);
    
    const double *rptr;
#pragma omp parallel for  default(shared) private(rptr,id, i1, i2, i3, r, comp,idn, vol, yA, yB, yC, yD, yE, yF, yG, yH, yABDC, yBDHF, yABFE, yCAEG, yDCGH, yEFHG, yM)
            for (i3=0; i3<m3; i3++){
                for (i2=0; i2<m2; i2++){
                    for (i1=0; i1<m1; i1++){
                        idn = i1+ i2* inc_y+ i3*inc_z;
                        // get edges, face stg and cell center of yc (shared --> private)
                        comp = 0;
                        for(r=0;r<3;r++){
                            rptr = yc+idn+comp;
                            // A
                            yA[r] = *rptr;
                            // B
                            rptr += inc_x;
                            yB[r] = *rptr;
                            // D
                            rptr += inc_y;
                            yD[r] = *rptr;
                            // C
                            rptr -= inc_x;
                            yC[r] = *rptr;
                            // G
                            rptr += inc_z;
                            yG[r] = *rptr;
                            // H
                            rptr += inc_x;
                            yH[r] = *rptr;
                            // F
                            rptr -= inc_y;
                            yF[r] = *rptr;
                            // E
                            rptr -= inc_x;
                            yE[r] = *rptr;
                            
                            yABDC[r] = 0.25*(yA[r]+yB[r]+yD[r]+yC[r]);
                            yBDHF[r] = 0.25*(yB[r]+yD[r]+yH[r]+yF[r]);
                            yABFE[r] = 0.25*(yA[r]+yB[r]+yF[r]+yE[r]);
                            yCAEG[r] = 0.25*(yC[r]+yA[r]+yE[r]+yG[r]);
                            yDCGH[r] = 0.25*(yD[r]+yC[r]+yG[r]+yH[r]);
                            yEFHG[r] = 0.25*(yE[r]+yF[r]+yH[r]+yG[r]);
                            yM[r]    = 0.5 *(yABDC[r]+yEFHG[r]);
                            
                            comp +=nn;
                        }
                        
                        // compute volumes of tetrahedra and sum up to voxel volume (private)
                        vol   = -det3D(yA, yB, yABDC, yM);
                        vol  += -det3D(yB, yD, yABDC, yM);
                        vol  += -det3D(yD, yC, yABDC, yM);
                        vol  += -det3D(yC, yA, yABDC, yM);
                        
                        vol  +=  det3D(yB, yD, yBDHF, yM);
                        vol  +=  det3D(yD, yH, yBDHF, yM);
                        vol  +=  det3D(yH, yF, yBDHF, yM);
                        vol  +=  det3D(yF, yB, yBDHF, yM);
                        
                        vol  +=  det3D(yA, yB, yABFE, yM);
                        vol  +=  det3D(yB, yF, yABFE, yM);
                        vol  +=  det3D(yF, yE, yABFE, yM);
                        vol  +=  det3D(yE, yA, yABFE, yM);
                        
                        vol  +=  det3D(yC, yA, yCAEG, yM);
                        vol  +=  det3D(yA, yE, yCAEG, yM);
                        vol  +=  det3D(yE, yG, yCAEG, yM);
                        vol  +=  det3D(yG, yC, yCAEG, yM);
                        
                        vol  +=  det3D(yD, yC, yDCGH, yM);
                        vol  +=  det3D(yC, yG, yDCGH, yM);
                        vol  +=  det3D(yG, yH, yDCGH, yM);
                        vol  +=  det3D(yH, yD, yDCGH, yM);
                        
                        vol  +=  det3D(yE, yF, yEFHG, yM);
                        vol  +=  det3D(yF, yH, yEFHG, yM);
                        vol  +=  det3D(yH, yG, yEFHG, yM);
                        vol  +=  det3D(yG, yE, yEFHG, yM);
                        
                        // write into result vector (private --> shared)
                        id =  i1+m1 *(i2+i3*m2);
                        volu[id] = vol;
                    }
                }
            }
            return;
}

double dVxTriangle2D(double *y1, double *y2, double *yM, double *x1, double *x2, double *xM){
    /*
     *  dV * x = tr(cofDy^T * Dx)
     */
    double cofDy[4];
    cofactor2D(cofDy, y1, y2, yM, 1.0);
    return (cofDy[0] * (x1[0]-xM[0]) + cofDy[2] * (x1[1]-xM[1]) + cofDy[1] * (x2[0]-xM[0]) + cofDy[3] * (x2[1]-xM[1]));
}

double dVxTetra3D(double *cofDy, const double *y1, const double *y2,const double *y3, const double *yM, const double *x1,const double *x2, const double *x3,const double *xM){
    cofactor3D(cofDy, y1, y2, y3, yM,1.0);
    return   cofDy[0] * (x1[0]-xM[0]) + cofDy[3] * (x1[1]-xM[1]) + cofDy[6] * (x1[2]-xM[2])
    + cofDy[1] * (x2[0]-xM[0]) + cofDy[4] * (x2[1]-xM[1]) + cofDy[7] * (x2[2]-xM[2])
    + cofDy[2] * (x3[0]-xM[0]) + cofDy[5] * (x3[1]-xM[1]) + cofDy[8] * (x3[2]-xM[2]);
}

void dJacx2D( double *result, const double *yc,const double *x,const double *m, int n){
    double yA[2], yB[2], yC[2], yD[2], yM[2], vol;
    double xA[2], xB[2], xC[2], xD[2], xM[2];
    int i1, i2, r,id, idn,comp;
    double res;
    
    const int m1 = m[0], m2 = m[1];
    const int nn = (m1+1)*(m2+1);
    // increments
    const int inc_x=1;
    const int inc_y= m1+1;
    //pointers
    const double *ptr;

#pragma omp parallel for  default(shared) private(ptr,i1, i2, r,id, idn,comp,res, yA, yB, yC, yD, yM, xA, xB, xC, xD, xM)
        for (i2=0; i2<m2; i2++){
            for (i1=0; i1<m1; i1++){
                // get edges and cell center of yc (shared --> private)
                idn = i1+i2*inc_y;
                comp = 0;
                for (r=0; r<2; r++) {
                    // get y at edges
                    ptr = yc + idn+comp;
                    yA[r] = *ptr; // yc[idn           +comp];
                    ptr = ptr+inc_x;
                    yB[r] = *ptr; //yc[idn+1         +comp];
                    ptr = ptr+inc_y;
                yD[r] = *ptr; //yc[idn+1+(m1+1)  +comp];
                ptr = ptr-inc_x;
                yC[r] = *ptr; //* yc[idn  +(m1+1)  +comp];
                yM[r] = 0.25*(yA[r]+yB[r]+yC[r]+yD[r]);
                
                // get x at edges
                ptr = x + idn+comp;
                xA[r] = *ptr; // yc[idn           +comp];
                ptr = ptr+inc_x;
                xB[r] = *ptr; //yc[idn+1         +comp];
                ptr = ptr+inc_y;
                xD[r] = *ptr; //yc[idn+1+(m1+1)  +comp];
                ptr = ptr-inc_x;
                xC[r] = *ptr; //* yc[idn  +(m1+1)  +comp];
                xM[r] = 0.25*(xA[r]+xB[r]+xC[r]+xD[r]);
                
                comp+=nn;
                
            }
            // compute dV*x on triangles (private)
            res  = dVxTriangle2D(yA, yB, yM, xA, xB, xM);
            res += dVxTriangle2D(yB, yD, yM, xB, xD, xM);
            res += dVxTriangle2D(yD, yC, yM, xD, xC, xM);
            res += dVxTriangle2D(yC, yA, yM, xC, xA, xM);
            
            // write into result vector (private --> shared)
            id = i1+i2*m1;
            result[id] = res;
        }
    }
}

void dJacx3D( double *result, const double *yc,const double *x,const double *m,const int n){
    const int m1 = m[0], m2 = m[1],m3 = m[2];
    double yA[3], yB[3], yC[3], yD[3], yE[3], yF[3], yG[3], yH[3], yABDC[3], yBDHF[3], yABFE[3], yCAEG[3], yDCGH[3], yEFHG[3], yM[3];
    double xA[3], xB[3], xC[3], xD[3], xE[3], xF[3], xG[3], xH[3], xABDC[3], xBDHF[3], xABFE[3], xCAEG[3], xDCGH[3], xEFHG[3], xM[3];
    double cofDy[9];
    double res;
    
    int id, i1, i2, i3, r, comp,idn;
    
    const int nn = (m1+1)*(m2+1)*(m3+1);
    const int inc_x = 1;
    const int inc_y = m1+1;
    const int inc_z = (m1+1)*(m2+1);
    
    const double *rptr;
    
#pragma omp parallel for  default(shared) private(rptr, id, i1, i2, i3, r, comp,idn, res, cofDy, yA, yB, yC, yD, yE, yF, yG, yH, yABDC, yBDHF, yABFE, yCAEG, yDCGH, yEFHG, yM, xA, xB, xC, xD, xE, xF, xG, xH, xABDC, xBDHF, xABFE, xCAEG, xDCGH, xEFHG, xM)
            for (i3=0; i3<m3; i3++){
                for (i2=0; i2<m2; i2++){
                    for (i1=0; i1<m1; i1++){
                        idn = i1+ i2* inc_y+ i3*inc_z;
                        // get edges, face stg and cell center of yc (shared --> private)
                        comp = 0;
                        for(r=0;r<3;r++){
                            rptr = yc+idn+comp;
                            // A
                            yA[r] = *rptr;
                            // B
                            rptr += inc_x;
                            yB[r] = *rptr;
                            // D
                            rptr += inc_y;
                            yD[r] = *rptr;
                            // C
                            rptr -= inc_x;
                            yC[r] = *rptr;
                            // G
                            rptr += inc_z;
                            yG[r] = *rptr;
                            // H
                            rptr += inc_x;
                            yH[r] = *rptr;
                            // F
                            rptr -= inc_y;
                            yF[r] = *rptr;
                            // E
                            rptr -= inc_x;
                            yE[r] = *rptr;
                            
                            yABDC[r] = 0.25*(yA[r]+yB[r]+yD[r]+yC[r]);
                            yBDHF[r] = 0.25*(yB[r]+yD[r]+yH[r]+yF[r]);
                            yABFE[r] = 0.25*(yA[r]+yB[r]+yF[r]+yE[r]);
                            yCAEG[r] = 0.25*(yC[r]+yA[r]+yE[r]+yG[r]);
                            yDCGH[r] = 0.25*(yD[r]+yC[r]+yG[r]+yH[r]);
                            yEFHG[r] = 0.25*(yE[r]+yF[r]+yH[r]+yG[r]);
                            yM[r]    = 0.5 *(yABDC[r]+yEFHG[r]);
                            
                            rptr = x+idn+comp;
                            // A
                            xA[r] = *rptr;
                            // B
                            rptr += inc_x;
                            xB[r] = *rptr;
                            // D
                            rptr += inc_y;
                            xD[r] = *rptr;
                            // C
                            rptr -= inc_x;
                            xC[r] = *rptr;
                            // G
                            rptr += inc_z;
                            xG[r] = *rptr;
                            // H
                            rptr += inc_x;
                            xH[r] = *rptr;
                            // F
                            rptr -= inc_y;
                            xF[r] = *rptr;
                            // E
                            rptr -= inc_x;
                            xE[r] = *rptr;
                            
                            xABDC[r] = 0.25*(xA[r]+xB[r]+xD[r]+xC[r]);
                            xBDHF[r] = 0.25*(xB[r]+xD[r]+xH[r]+xF[r]);
                            xABFE[r] = 0.25*(xA[r]+xB[r]+xF[r]+xE[r]);
                            xCAEG[r] = 0.25*(xC[r]+xA[r]+xE[r]+xG[r]);
                            xDCGH[r] = 0.25*(xD[r]+xC[r]+xG[r]+xH[r]);
                            xEFHG[r] = 0.25*(xE[r]+xF[r]+xH[r]+xG[r]);
                            xM[r]    = 0.5 *(xABDC[r]+xEFHG[r]);
                            comp +=nn;
                        }
                        
                        
                        // compute dV*x on tetrahedra (private)
                        res = -dVxTetra3D(cofDy, yA, yB, yABDC, yM, xA, xB, xABDC, xM);
                        res+= -dVxTetra3D(cofDy, yB, yD, yABDC, yM, xB, xD, xABDC, xM);
                        res+= -dVxTetra3D(cofDy, yD, yC, yABDC, yM, xD, xC, xABDC, xM);
                        res+= -dVxTetra3D(cofDy, yC, yA, yABDC, yM, xC, xA, xABDC, xM);
                        
                        res+=  dVxTetra3D(cofDy, yB, yD, yBDHF, yM, xB, xD, xBDHF, xM);
                        res+=  dVxTetra3D(cofDy, yD, yH, yBDHF, yM, xD, xH, xBDHF, xM);
                        res+=  dVxTetra3D(cofDy, yH, yF, yBDHF, yM, xH, xF, xBDHF, xM);
                        res+=  dVxTetra3D(cofDy, yF, yB, yBDHF, yM, xF, xB, xBDHF, xM);
                        
                        res+=  dVxTetra3D(cofDy, yA, yB, yABFE, yM, xA, xB, xABFE, xM);
                        res+=  dVxTetra3D(cofDy, yB, yF, yABFE, yM, xB, xF, xABFE, xM);
                        res+=  dVxTetra3D(cofDy, yF, yE, yABFE, yM, xF, xE, xABFE, xM);
                        res+=  dVxTetra3D(cofDy, yE, yA, yABFE, yM, xE, xA, xABFE, xM);
                        
                        res+=  dVxTetra3D(cofDy, yC, yA, yCAEG, yM, xC, xA, xCAEG, xM);
                        res+=  dVxTetra3D(cofDy, yA, yE, yCAEG, yM, xA, xE, xCAEG, xM);
                        res+=  dVxTetra3D(cofDy, yE, yG, yCAEG, yM, xE, xG, xCAEG, xM);
                        res+=  dVxTetra3D(cofDy, yG, yC, yCAEG, yM, xG, xC, xCAEG, xM);
                        
                        res+=  dVxTetra3D(cofDy, yD, yC, yDCGH, yM, xD, xC, xDCGH, xM);
                        res+=  dVxTetra3D(cofDy, yC, yG, yDCGH, yM, xC, xG, xDCGH, xM);
                        res+=  dVxTetra3D(cofDy, yG, yH, yDCGH, yM, xG, xH, xDCGH, xM);
                        res+=  dVxTetra3D(cofDy, yH, yD, yDCGH, yM, xH, xD, xDCGH, xM);
                        
                        res+=  dVxTetra3D(cofDy, yE, yF, yEFHG, yM, xE, xF, xEFHG, xM);
                        res+=  dVxTetra3D(cofDy, yF, yH, yEFHG, yM, xF, xH, xEFHG, xM);
                        res+=  dVxTetra3D(cofDy, yH, yG, yEFHG, yM, xH, xG, xEFHG, xM);
                        res+=  dVxTetra3D(cofDy, yG, yE, yEFHG, yM, xG, xE, xEFHG, xM);
                        
                        // write into result vector (private --> shared)
                        id =  i1+m1 *(i2+i3*m2);
                        result[id]= res;
                    }
                }
            }
}


inline void dVTTriangle2D(double *y1, double *y2, double *yM, double x, double *r1, double *r2, double *rM){
    double cof[4];
    cofactor2D(cof, y1, y2, yM, x);
    // update result vectors
    r1[0] += cof[0];
    r1[1] += cof[2];
    
    r2[0] += cof[1];
    r2[1] += cof[3];
    
    rM[0] -= (cof[0]+cof[1]);
    rM[1] -= (cof[2]+cof[3]);
    
}

void dVTTetra3D(const double *yA,const double *yB,const double *yC,const double *yM,const double x, double *rA, double *rB, double *rC, double *rM){
    int i, j,k;
    double cof[9];
    cofactor3D(cof, yA, yB, yC, yM, x);
    
    /* update result vectors */
    k=0;
//    rM[0] -= (cof[0]+cof[1]+cof[2]);
//    rM[1] -= (cof[3]+cof[4]+cof[5]);
//    rM[2] -= (cof[6]+cof[7]+cof[8]);
    
    for (i=0; i<3; i++){
        for (j=0; j<3; j++){
            rM[i] -= cof[j+k];
        }
        rA[i] += cof[0+k];
        rB[i] += cof[1+k];
        rC[i] += cof[2+k];
        k+=3;
    }
    
}

void dJacTx2D( double *prod, const double *yc,const double *x,const double *m, int n){
    double yA[2], yB[2], yC[2], yD[2], yM[2], vol;
    double rA[2], rB[2], rC[2], rD[2], rM[2];
    double xp;
    int i1, i2, r,id, idn,comp,off1,off2;
    
    const int m1 = m[0], m2 = m[1];
    const int nn = (m1+1)*(m2+1);
    // increments
    const int inc_x=1;
    const int inc_y= m1+1;
    //pointers
    const double *ptr;
    double *wptr;
    for (off1=0; off1<2; off1++) {
        for (off2=0; off2<2; off2++) {
#pragma omp parallel for default(shared) private(ptr,wptr,i1, i2, r,id, idn,comp, xp, yA, yB, yC, yD, yM, rA, rB, rC, rD, rM)
            for (i2=off2; i2<m2; i2+=2){
                for (i1=off1; i1<m1; i1+=2){
                    idn = i1+i2*inc_y;
                    id = i1+i2*m1;
                    comp = 0;
                    for (r=0; r<2; r++) {
                        // get y at edges
                        ptr = yc + idn + comp;
                        yA[r] = *ptr; // yc[idn           +comp];
                        ptr = ptr+inc_x;
                        yB[r] = *ptr; //yc[idn+1         +comp];
                        ptr = ptr+inc_y;
                        yD[r] = *ptr; //yc[idn+1+(m1+1)  +comp];
                        ptr = ptr-inc_x;
                        yC[r] = *ptr; //* yc[idn  +(m1+1)  +comp];
                        
                        yM[r] = 0.25*(yA[r]+yB[r]+yC[r]+yD[r]);
                        
                        rA[r] = 0.0;
                        rB[r] = 0.0;
                        rC[r] = 0.0;
                        rD[r] = 0.0;
                        rM[r] = 0.0;
                        
                        comp+=nn;
                    }
                    
                    xp = x[id];
                    
                    // compute on triangles (private)
                    dVTTriangle2D(yA, yB, yM, xp, rA, rB, rM);
                    dVTTriangle2D(yB, yD, yM, xp, rB, rD, rM);
                    dVTTriangle2D(yD, yC, yM, xp, rD, rC, rM);
                    dVTTriangle2D(yC, yA, yM, xp, rC, rA, rM);
                    
                    // write into result vector (private --> shared)
                    comp = 0;
                    for (r=0; r < 2; r++){
                        wptr = prod + idn + comp;
                        rM[r] *= .25;
                        //dS[idn     +comp] += rA[r]+0.25*rM[r];
                        *wptr += rA[r]+rM[r];
                        wptr += inc_x;
                        //dS[idn +1   +comp] += rB[r]+rM[r];
                        *wptr += rB[r]+rM[r];
                        wptr +=  inc_y;
                        //dS[idn    +(m1+1)   +comp] += rC[r]+rM[r];
                        *wptr += rD[r]+rM[r];
                        wptr -= inc_x;
                        *wptr += rC[r]+rM[r];
                        //dS[idn +1 +(m1+1)   +comp] += rD[r]+rM[r];
                        comp+=nn;
                    }
                }
            }
        }
    }

}

void  dJacTx3D( double *prod, const double *yc,const double *x,const double *m, int n){
    const int m1 = m[0], m2 = m[1],m3 = m[2];
    double yA[3], yB[3], yC[3], yD[3], yE[3], yF[3], yG[3], yH[3], yABDC[3], yBDHF[3], yABFE[3], yCAEG[3], yDCGH[3], yEFHG[3], yM[3];
    double rA[3], rB[3], rC[3], rD[3], rE[3], rF[3], rG[3], rH[3], rABDC[3], rBDHF[3], rABFE[3], rCAEG[3], rDCGH[3], rEFHG[3], rM[3];
    double xp;
    
    int id, i1, i2, i3, r, comp,idn,off1,off2,off3;
    
    const int nn = (m1+1)*(m2+1)*(m3+1);
    const int inc_x = 1;
    const int inc_y = m1+1;
    const int inc_z = (m1+1)*(m2+1);
    
    const double *rptr;
    double *wptr;
    
    for (off1=0; off1<2; off1++) {
        for (off2=0; off2<2; off2++) {
            for (off3=0; off3<2; off3++) {
#pragma omp parallel for  default(shared) private(rptr,wptr,id, i1, i2, i3, r, comp,idn, xp, yA, yB, yC, yD, yE, yF, yG, yH, yABDC, yBDHF, yABFE, yCAEG, yDCGH, yEFHG, yM, rA, rB, rC, rD, rE, rF, rG, rH, rABDC, rBDHF, rABFE, rCAEG, rDCGH, rEFHG, rM)
                for (i3=off3; i3<m3; i3+=2){
                    for (i2=off2; i2<m2; i2+=2){
                        for (i1=off1; i1<m1; i1+=2){
                            idn = i1+ i2* inc_y+ i3*inc_z;
                            id =  i1+m1 *(i2+i3*m2);
                            // get edges, face stg and cell center of yc (shared --> private)
                            comp = 0;
                            for(r=0;r<3;r++){
                                rptr = yc+idn+comp;
                                // A
                                yA[r] = *rptr;
                                // B
                                rptr += inc_x;
                                yB[r] = *rptr;
                                // D
                                rptr += inc_y;
                                yD[r] = *rptr;
                                // C
                                rptr -= inc_x;
                                yC[r] = *rptr;
                                // G
                                rptr += inc_z;
                                yG[r] = *rptr;
                                // H
                                rptr += inc_x;
                                yH[r] = *rptr;
                                // F
                                rptr -= inc_y;
                                yF[r] = *rptr;
                                // E
                                rptr -= inc_x;
                                yE[r] = *rptr;
                                
                                yABDC[r] = 0.25*(yA[r]+yB[r]+yD[r]+yC[r]);
                                yBDHF[r] = 0.25*(yB[r]+yD[r]+yH[r]+yF[r]);
                                yABFE[r] = 0.25*(yA[r]+yB[r]+yF[r]+yE[r]);
                                yCAEG[r] = 0.25*(yC[r]+yA[r]+yE[r]+yG[r]);
                                yDCGH[r] = 0.25*(yD[r]+yC[r]+yG[r]+yH[r]);
                                yEFHG[r] = 0.25*(yE[r]+yF[r]+yH[r]+yG[r]);
                                yM[r]    = 0.5 *(yABDC[r]+yEFHG[r]);
                                
                                comp +=nn;
                                
                                rA[r] = 0.0;
                                rB[r] = 0.0;
                                rC[r] = 0.0;
                                rD[r] = 0.0;
                                rE[r] = 0.0;
                                rF[r] = 0.0;
                                rG[r] = 0.0;
                                rH[r] = 0.0;
                                rABDC[r] = 0.0;
                                rBDHF[r] = 0.0;
                                rABFE[r] = 0.0;
                                rCAEG[r] = 0.0;
                                rDCGH[r] = 0.0;
                                rEFHG[r] = 0.0;
                                rM[r] = 0.0;
                            }
                            xp = x[id];
                            
                            
                            /* compute */
                            dVTTetra3D(yA, yB, yABDC, yM, -xp, rA, rB, rABDC, rM);
                            dVTTetra3D(yB, yD, yABDC, yM, -xp, rB, rD, rABDC, rM);
                            dVTTetra3D(yD, yC, yABDC, yM, -xp, rD, rC, rABDC, rM);
                            dVTTetra3D(yC, yA, yABDC, yM, -xp, rC, rA, rABDC, rM);
                            
                            dVTTetra3D(yB, yD, yBDHF, yM, xp, rB, rD, rBDHF, rM);
                            dVTTetra3D(yD, yH, yBDHF, yM, xp, rD, rH, rBDHF, rM);
                            dVTTetra3D(yH, yF, yBDHF, yM, xp, rH, rF, rBDHF, rM);
                            dVTTetra3D(yF, yB, yBDHF, yM, xp, rF, rB, rBDHF, rM);
                            
                            dVTTetra3D(yA, yB, yABFE, yM, xp, rA, rB, rABFE, rM);
                            dVTTetra3D(yB, yF, yABFE, yM, xp, rB, rF, rABFE, rM);
                            dVTTetra3D(yF, yE, yABFE, yM, xp, rF, rE, rABFE, rM);
                            dVTTetra3D(yE, yA, yABFE, yM, xp, rE, rA, rABFE, rM);
                            
                            dVTTetra3D(yC, yA, yCAEG, yM, xp, rC, rA, rCAEG, rM);
                            dVTTetra3D(yA, yE, yCAEG, yM, xp, rA, rE, rCAEG, rM);
                            dVTTetra3D(yE, yG, yCAEG, yM, xp, rE, rG, rCAEG, rM);
                            dVTTetra3D(yG, yC, yCAEG, yM, xp, rG, rC, rCAEG, rM);
                            
                            dVTTetra3D(yD, yC, yDCGH, yM, xp, rD, rC, rDCGH, rM);
                            dVTTetra3D(yC, yG, yDCGH, yM, xp, rC, rG, rDCGH, rM);
                            dVTTetra3D(yG, yH, yDCGH, yM, xp, rG, rH, rDCGH, rM);
                            dVTTetra3D(yH, yD, yDCGH, yM, xp, rH, rD, rDCGH, rM);
                            
                            dVTTetra3D(yE, yF, yEFHG, yM, xp, rE, rF, rEFHG, rM);
                            dVTTetra3D(yF, yH, yEFHG, yM, xp, rF, rH, rEFHG, rM);
                            dVTTetra3D(yH, yG, yEFHG, yM, xp, rH, rG, rEFHG, rM);
                            dVTTetra3D(yG, yE, yEFHG, yM, xp, rG, rE, rEFHG, rM);
                            
                            // write into result vector (private --> shared)
                            comp = 0;
                            for (r=0; r < 3; r++){
                                // scale vecs
                                rABDC[r] *= 0.25;
                                rBDHF[r] *= 0.25;
                                rABFE[r] *= 0.25;
                                rCAEG[r] *= 0.25;
                                rDCGH[r] *= 0.25;
                                rEFHG[r] *= 0.25;
                                rM[r]    *= 0.125;
                                // get pointer to A
                                wptr = prod + idn + comp;
                                *wptr +=  rA[r]+ rABDC[r]+rABFE[r]+rCAEG[r] + rM[r];
                                // B
                                wptr += inc_x;
                                *wptr +=  rB[r]+ rABDC[r]+rBDHF[r]+rABFE[r]+rM[r];
                                // D
                                wptr += inc_y;
                                *wptr += rD[r]+rABDC[r]+rBDHF[r]+rDCGH[r]+rM[r];
                                // C
                                wptr -= inc_x;
                                *wptr += rC[r]+rABDC[r]+rCAEG[r]+rDCGH[r]+rM[r];
                                // G
                                wptr += inc_z;
                                *wptr += rG[r]+rCAEG[r]+rDCGH[r]+rEFHG[r]+rM[r];
                                // H
                                wptr += inc_x;
                                *wptr += rH[r]+rBDHF[r]+rDCGH[r]+rEFHG[r]+rM[r];
                                // F
                                wptr -= inc_y;
                                *wptr += rF[r]+rBDHF[r]+rABFE[r]+rEFHG[r]+rM[r];
                                // E
                                wptr -= inc_x;
                                *wptr += rE[r]+rABFE[r]+rCAEG[r]+rEFHG[r]+rM[r];
                                comp +=nn;
                            }
                        }
                    }
                }
            }
        }
    }
    

}

#pragma mark -
#pragma mark hyperElastic two-D
inline double SvolTriangle2D(double *y1, double *y2, double *yM,const double VRef, bool doDerivative, double *r1, double *r2, double *rM,double alphaVolume=1.0){
    double vol =  volTriangle2D(y1, y2, yM)/VRef;
    double S = alphaVolume * phicaley(vol);
    if (doDerivative) {
        double dS = dphicaley(vol);
        double cof[4];
        cofactor2D(cof, y1, y2, yM, alphaVolume * dS/(VRef));
        
        // update result vectors: back projection
        r1[0] += cof[0];
        r1[1] += cof[2];
        
        r2[0] += cof[1];
        r2[1] += cof[3];
        
        rM[0] -= (cof[0]+cof[1]);
        rM[1] -= (cof[2]+cof[3]);
        
    }
    return S;
}
void S2D(double *S, double *dS, const double *yc,const double *VRef, const double *alpha ,const double *h,const double *m,const int n, const bool doDerivative){
    int id,idn, i1, i2, r, comp,off1,off2;
    double yA[2], yB[2], yC[2], yD[2], yM[2];
    double rA[2], rB[2], rC[2], rD[2], rM[2];
    double sum=0.0,sumvoxel;
    
    const int m1 = m[0], m2 = m[1];
    const int nn = (m1+1)*(m2+1);
    const double hd = h[0] * h[1];
    const double alphaVolume = alpha[2]*hd;
    const int inc_x=1;
    const int inc_y= m1+1;
    const double *ptr;
    double *wptr;
    if (alphaVolume >0 ) {
        
        for (off1=0; off1<2; off1++) {
            for (off2=0; off2<2; off2++) {
#pragma omp parallel for default(shared) private(wptr,ptr,id,idn, i1, i2, r, comp,sumvoxel, yA, yB, yC, yD, yM, rA, rB, rC, rD, rM) reduction(+: sum)
                    for (i2=off2; i2<m2; i2+=2){
                        for (i1=off1; i1<m1; i1+=2){
                            // get edges and cell center of yc (shared --> private)
                            idn = i1+i2*inc_y;
                            comp = 0;
                            for (r=0; r<2; r++) {
                                ptr = yc + idn+comp;
                                yA[r] = *ptr; // yc[idn           +comp];
                                ptr = ptr+inc_x;
                                yB[r] = *ptr; //yc[idn+1         +comp];
                                ptr = ptr+inc_y;
                                yD[r] = *ptr; //yc[idn+1+(m1+1)  +comp];
                                ptr = ptr-inc_x;
                                yC[r] = *ptr; //* yc[idn  +(m1+1)  +comp];
                                yM[r] = 0.25*(yA[r]+yB[r]+yC[r]+yD[r]);
                                comp+=nn;
                                // init private result vectors (private)
                                rA[r] = 0.0;
                                rB[r] = 0.0;
                                rC[r] = 0.0;
                                rD[r] = 0.0;
                                rM[r] = 0.0;
                                
                            }
                            // do the work (private)
                            sumvoxel  = SvolTriangle2D(yA, yB, yM, VRef[0], doDerivative, rA, rB, rM,alphaVolume);
                            sumvoxel += SvolTriangle2D(yB, yD, yM, VRef[1], doDerivative, rB, rD, rM,alphaVolume);
                            sumvoxel += SvolTriangle2D(yD, yC, yM, VRef[2], doDerivative, rD, rC, rM,alphaVolume);
                            sumvoxel += SvolTriangle2D(yC, yA, yM, VRef[3], doDerivative, rC, rA, rM,alphaVolume);
                            
                            // write into result vector (private --> shared)
                            sum += sumvoxel;
                            if (doDerivative) {
                                /* write result into prod */
                                comp = 0;
                                for (r=0; r < 2; r++){
                                    wptr = dS + idn + comp;
                                    rM[r] *= .25;
                                    //dS[idn     +comp] += rA[r]+0.25*rM[r];
                                    *wptr += rA[r]+rM[r];
                                    wptr += inc_x;
                                    //dS[idn +1   +comp] += rB[r]+rM[r];
                                    *wptr += rB[r]+rM[r];
                                    wptr +=  inc_y;
                                    //dS[idn    +(m1+1)   +comp] += rC[r]+rM[r];
                                    *wptr += rD[r]+rM[r];
                                    wptr -= inc_x;
                                    *wptr += rC[r]+rM[r];
                                    //dS[idn +1 +(m1+1)   +comp] += rD[r]+rM[r];
                                    comp+=nn;
                                }
                            }
                        }
                    }
                }             
            }
    }
    *S = sum;
}
void d2SvolTriangle2D(double *y1, double *y2, double *yM, double *x1, double *x2, double *xM,
        const double VRef, double *r1, double *r2, double *rM,double alphaVolume){
    double dVx, vol, d2S;
    int i, j,k;
    dVx = dVxTriangle2D(y1, y2, yM, x1, x2, xM);
    vol =  volTriangle2D(y1, y2, yM);
    d2S = d2phicaley(vol/VRef)/(VRef*VRef);
    
    double cof[4];
    cofactor2D(cof, y1, y2, yM,alphaVolume* (d2S*dVx));
    
    /* update result vectors */
    k=0;
    for (i=0; i<2; i++){
        for (j=0; j<2; j++){
            rM[i] -= cof[j+k];
        }
        r1[i] += cof[0+k];
        r2[i] += cof[1+k];
        k+=2;
    }
}
void d2S2D(double *prod, const double *yc,const double *x,const double *VRef,const double *alpha,const double *h,const double *m, const int n){
    int idn, i1, i2, r, comp,off1,off2;
    double yA[2], yB[2], yC[2], yD[2], yM[2];
    double xA[2], xB[2], xC[2], xD[2], xM[2];
    double rA[2], rB[2], rC[2], rD[2], rM[2];
    
    const int m1 = m[0], m2 = m[1];
    const int nn = (m1+1)*(m2+1);
    const double hd = h[0] * h[1];
    const double alphaVolume = alpha[2]*hd;
    // increments along x and y
    const int inc_x=1;
    const int inc_y= m1+1;
    // pointers for reading / writing
    const double *ptr;
    double *wptr;
    
    for (off1=0; off1<2; off1++) {
        for (off2=0; off2<2; off2++) {
#pragma omp parallel private(wptr,ptr,idn, i1, i2, r, comp,yA,yB,yC,yD,yM,xA,xB,xC,xD,xM,rA,rB,rC,rD,rM)
            {
#pragma omp for
                for (i2=off2; i2<m2; i2+=2){
                    for (i1=off1; i1<m1; i1+=2){
                        idn = i1+i2*inc_y;
                        comp = 0;
                        for (r=0; r<2; r++) {
                            // get y at edges
                            ptr = yc + idn+comp;
                            yA[r] = *ptr; // yc[idn           +comp];
                            ptr = ptr+inc_x;
                            yB[r] = *ptr; //yc[idn+1         +comp];
                            ptr = ptr+inc_y;
                            yD[r] = *ptr; //yc[idn+1+(m1+1)  +comp];
                            ptr = ptr-inc_x;
                            yC[r] = *ptr; //* yc[idn  +(m1+1)  +comp];
                            yM[r] = 0.25*(yA[r]+yB[r]+yC[r]+yD[r]);
                            
                            // get x at edges
                            ptr = x + idn+comp;
                            xA[r] = *ptr; // yc[idn           +comp];
                            ptr = ptr+inc_x;
                            xB[r] = *ptr; //yc[idn+1         +comp];
                            ptr = ptr+inc_y;
                            xD[r] = *ptr; //yc[idn+1+(m1+1)  +comp];
                            ptr = ptr-inc_x;
                            xC[r] = *ptr; //* yc[idn  +(m1+1)  +comp];
                            xM[r] = 0.25*(xA[r]+xB[r]+xC[r]+xD[r]);
                            
                            // init private result vectors (private)
                            rA[r] = 0.0;
                            rB[r] = 0.0;
                            rC[r] = 0.0;
                            rD[r] = 0.0;
                            rM[r] = 0.0;
                            comp+=nn;
                            
                        }
                        
                        // do the work (private)
                        d2SvolTriangle2D(yA, yB, yM, xA, xB, xM, VRef[0], rA, rB, rM,alphaVolume);
                        d2SvolTriangle2D(yB, yD, yM, xB, xD, xM, VRef[1], rB, rD, rM,alphaVolume);
                        d2SvolTriangle2D(yD, yC, yM, xD, xC, xM, VRef[2], rD, rC, rM,alphaVolume);
                        d2SvolTriangle2D(yC, yA, yM, xC, xA, xM, VRef[3], rC, rA, rM,alphaVolume);
                        
                        // write into result vector (private --> shared)
                        comp = 0;
                        for (r=0; r < 2; r++){
                            wptr = prod + idn + comp;
                            rM[r] *= .25;
                            //dS[idn     +comp] += rA[r]+0.25*rM[r];
                            *wptr += rA[r]+rM[r];
                            wptr += inc_x;
                            //dS[idn +1   +comp] += rB[r]+rM[r];
                            *wptr += rB[r]+rM[r];
                            wptr +=  inc_y;
                            //dS[idn    +(m1+1)   +comp] += rC[r]+rM[r];
                            *wptr += rD[r]+rM[r];
                            wptr -= inc_x;
                            *wptr += rC[r]+rM[r];
                            //dS[idn +1 +(m1+1)   +comp] += rD[r]+rM[r];
                            comp+=nn;
                        }
                    }
                }
        }
    }

}

}
void d2SdiagTriangle(double *y1, double *y2, double *yM, const double VRef, double *r1, double *r2, double *r3, double *r4, double alphaVolume){
    double dVx, vol, d2S;
    // things we need anyway for this triangle
    vol =  volTriangle2D(y1, y2, yM);
    d2S = alphaVolume* d2phicaley(vol/VRef)/(VRef*VRef);
    double cof[4];
    cofactor2D(cof, y1, y2, yM,1.0);
    
    //dVx = (cof[0] * (x1[0]-xM[0]) + cof[2] * (x1[1]-xM[1]) + cof[1] * (x2[0]-xM[0]) + cof[3] * (x2[1]-xM[1]));
    
    // now we have 4 cases inside this triangle
    //  1) x1[0] = 1 --> xM[0] = 1/4 --> write to r1[0]
    dVx = cof[0] * (0.75) + cof[1] * (-0.25);
    r1[0] +=  dVx * d2S * dVx;
    
    //  2) x1[1] = 1 --> xM[1] = 1/4 --> write to r1[1]
    dVx = cof[2] * (0.75) + cof[3] * (-0.25);
    r1[1] += dVx * d2S * dVx;
    
    //  3) x2[0] = 1 --> xM[0] = 1/4 --> write to r2[0]
    dVx = cof[0] * (-0.25) +  cof[1] * (0.75);
    r2[0] += dVx * d2S * dVx;
    
    //  4) x2[1] = 1 --> xM[1] = 1/4 --> write to r2[1]
    dVx =  cof[2] * (-0.25) + cof[3] * (0.75) ;
    r2[1] += dVx * d2S * dVx;;
    
    
    // and we have 2 cases outside this triangle
    //  1) x3[0] = 1 or x4[0] = 1--> xM[0] = 1/4 --> write to r3[0] and r4[0]
    dVx = (cof[0] + cof[1]) * (-0.25);
    dVx = dVx * d2S * dVx;
    r3[0] +=  dVx;
    r4[0] += dVx;
    
    //  2) x3[1] = 1 or x4[1] --> xM[1] = 1/4 --> write to r3[1] and r4[1]
    dVx = (cof[2]  + cof[3]) * (-0.25);
    dVx = dVx * d2S * dVx;
    r3[1] += dVx ;
    r4[1] += dVx;
    
}

void d2Sdiag2D(double *prod, const double *yc,const double *VRef,const double *alpha,const double *h,const double *m, const int n){
    int idn, i1, i2, r, comp,off1,off2;
    double yA[2], yB[2], yC[2], yD[2], yM[2];
    double rA[2], rB[2], rC[2], rD[2];
    
    const int m1 = m[0], m2 = m[1];
    const int nn = (m1+1)*(m2+1);
    const double hd = h[0] * h[1];
    const double alphaVolume = alpha[2]*hd;
    // increments along x and y
    const int inc_x=1;
    const int inc_y= m1+1;
    // pointers for reading / writing
    const double *ptr;
    double *wptr;
    
    for (off1=0; off1<2; off1++) {
        for (off2=0; off2<2; off2++) {
#pragma omp parallel private(wptr,ptr,idn, i1, i2, r, comp,yA,yB,yC,yD,yM,rA,rB,rC,rD)
            {
#pragma omp for
                for (i2=off2; i2<m2; i2+=2){
                    for (i1=off1; i1<m1; i1+=2){
                        idn = i1+i2*(m1+1);
                        comp = 0;
                        for (r=0; r<2; r++) {
                            // get y at edges
                            ptr = yc + idn+comp;
                            yA[r] = *ptr; // yc[idn           +comp];
                            ptr = ptr+inc_x;
                            yB[r] = *ptr; //yc[idn+1         +comp];
                            ptr = ptr+inc_y;
                            yD[r] = *ptr; //yc[idn+1+(m1+1)  +comp];
                            ptr = ptr-inc_x;
                            yC[r] = *ptr; //* yc[idn  +(m1+1)  +comp];
                            yM[r] = 0.25*(yA[r]+yB[r]+yC[r]+yD[r]);
                            
                            // init private result vectors (private)
                            rA[r] = 0.0;
                            rB[r] = 0.0;
                            rC[r] = 0.0;
                            rD[r] = 0.0;
                            comp+=nn;
                            
                        }
                        
                        // do the work (private)
                        d2SdiagTriangle(yA, yB, yM, VRef[0], rA, rB,rD,rC,alphaVolume);
                        d2SdiagTriangle(yB, yD, yM, VRef[1], rB, rD,rC,rA,alphaVolume);
                        d2SdiagTriangle(yD, yC, yM, VRef[2], rD, rC,rA,rB,alphaVolume);
                        d2SdiagTriangle(yC, yA, yM, VRef[3], rC, rA,rB,rD,alphaVolume);
                        
                        // write into result vector (private --> shared)
                        comp = 0;
                        for (r=0; r < 2; r++){
                            wptr = prod + idn + comp;
                            //dS[idn     +comp] += rA[r]+0.25*rM[r];
                            *wptr += rA[r];
                            wptr += inc_x;
                            //dS[idn +1   +comp] += rB[r]+rM[r];
                            *wptr += rB[r];
                            wptr +=  inc_y;
                            //dS[idn    +(m1+1)   +comp] += rC[r]+rM[r];
                            *wptr += rD[r];
                            wptr -= inc_x;
                            *wptr += rC[r];
                            //dS[idn +1 +(m1+1)   +comp] += rD[r]+rM[r];
                            comp+=nn;
                        }
                    }
                }
            }
        }
    }
 
    
}

#pragma mark -
#pragma mark hyperElastic three-D
double arTriangle3D(double *y1, double *y2, double *yM,const double ARef, double *r1, double *r2,
        double *rM, bool doDerivative,const double alphaArea){
    /*
     * The triangle is given by the nodal edges y1 and y2 and the cell-center yM
     *
     * The area is given as the euclidean norm of the crossproduct
     *
     * cprod = (y1-yM) x (y2-yM) := v x w
     *
     */
    double v[3], w[3], cprod[3], A;
    vecminus(v, y1, yM, 3);
    vecminus(w, y2, yM, 3);
    
    // cross-product
    cprod[0] =  v[1] * w[2] - v[2] * w[1];
    cprod[1] =  v[0] * w[2] - v[2] * w[0];
    cprod[2] =  v[0] * w[1] - v[1] * w[0];
    
    
    A= cprod[0]*cprod[0] + cprod[1]*cprod[1] + cprod[2]*cprod[2]; // area of triangle
    
    A = A/ARef; // change of surface area
    
    if (doDerivative) {
        double dA = alphaArea *  dphiarea(A)/ARef;
        double dwA, dvA;
        
        // first component
        dvA =  2.0 * (cprod[1] * w[2] + cprod[2] * w[1]); // d_{v[0]} A
        dwA = -2.0 * (cprod[1] * v[2] + cprod[2] * v[1]); // d_{w[0]} A
        
        rM[0] += dA * (-dvA - dwA);
        r1[0] += dA * ( dvA      );
        r2[0] += dA * (       dwA);
        
        // second component
        dvA = 2.0 * (cprod[0] * w[2] - cprod[2] * w[0]); // d_{v[1]} A
        dwA = 2.0 * (-cprod[0] * v[2] + cprod[2] * v[0]); //d_{w[1]} A
        
        rM[1] += dA * (-dvA - dwA);
        r1[1] += dA * ( dvA      );
        r2[1] += dA * (     + dwA);
        
        // third component
        dvA = 2.0 * (-cprod[0] * w[1] - cprod[1] * w[0]); // d_{v[2]} A
        dwA = 2.0 * (cprod[0] * v[1]  + cprod[1] * v[0]); // d_{w[2]} A
        rM[2] += dA * (-dvA - dwA);
        r1[2] += dA * ( dvA      );
        r2[2] += dA * (     + dwA);
        
    }
    return alphaArea * phiarea(A);
}

double areaTriangle3D(double *y1, double *y2, double *yM){
    double v[3], w[3], cprod[3];
    vecminus(v, y1, yM, 3);
    vecminus(w, y2, yM, 3);
    
    cprod[0] =  v[1] * w[2] - v[2] * w[1];
    cprod[1] =  v[0] * w[2] - v[2] * w[0];
    cprod[2] =  v[0] * w[1] - v[1] * w[0];
    return  cprod[0]*cprod[0] + cprod[1]*cprod[1] + cprod[2]*cprod[2];
}

double SvolTetra3D(double *yA, double *yB, double *yC, double *yM, const double VRef, bool doDerivative, double *rA, double *rB, double *rC, double *rM, double inv,const double alphaVolume){
    double vol = inv* det3D(yA, yB, yC, yM)/VRef;
    double S = alphaVolume *  phicaley(vol);
    if (doDerivative) {
        int i, j,k;
        double cof[9];
        cofactor3D(cof, yA, yB, yC, yM,inv* alphaVolume*(dphicaley(vol)/(VRef)));
        // update result vectors
        k=0;
        for (i=0; i<3; i++){
            for (j=0; j<3; j++){
                rM[i] -= cof[j+k];
            }
            rA[i] += cof[0+k];
            rB[i] += cof[1+k];
            rC[i] += cof[2+k];
            k+=3;
        }
    }
    return S;
}
inline void d2SvolTetra3D(double *yA, double *yB, double *yC, double *yM, double *xA, double *xB, double *xC, double *xM,const double VRef, double *rA, double *rB, double *rC, double *rM, const double inv,const double alphaVolume){
    double dVx, vol,vol2, d2S, cof[9];
    int i, j,k;
    dVx = dVxTetra3D(cof, yA, yB, yC, yM, xA, xB, xC, xM);
    vol = inv*( (yA[0]-yM[0]) * cof[0] + (yB[0]-yM[0]) * cof[1] + (yC[0]-yM[0]) * cof[2]);
    d2S = alphaVolume * d2phicaley(vol/VRef)/(VRef*VRef) * dVx;
    // update result vectors
    k=0;
    for (i=0; i<3; i++){
        for (j=0; j<3; j++){
            rM[i] -= cof[j+k]*d2S;
        }
        rA[i] += cof[k]*d2S;
        rB[i] += cof[1+k]*d2S;
        rC[i] += cof[2+k]*d2S;
        k+=3;
    }
}
inline void d2SarTriangle3D( double *y1,double *y2, double *yM, double *x1, double *x2, double *xM,const double ARef, double *r1, double *r2, double *rM,const double alphaArea){
    double vy[3], wy[3], wx[3], vx[3],cprod[3],A;
    
    vecminus(vy, y1, yM, 3);  // vy = y1 - yM;
    vecminus(wy, y2, yM, 3);  // wx = y2 - yM;
    vecminus(vx, x1, xM, 3);  // vx = x1 - xM;
    vecminus(wx, x2, xM, 3);  // wx = x2 - xM;
    
    // compute cross-product for y
    cprod[0] =  vy[1] * wy[2] - vy[2] * wy[1];
    cprod[1] =  vy[0] * wy[2] - vy[2] * wy[0];
    cprod[2] =  vy[0] * wy[1] - vy[1] * wy[0];
    
    A= cprod[0]*cprod[0] + cprod[1]*cprod[1] + cprod[2]*cprod[2];
    
    double dwA[3], dvA[3];
    // first component
    dvA[0] =  2.0 * (cprod[1] * wy[2] + cprod[2] * wy[1]); // d_{v[0]} A
    dwA[0] = -2.0 * (cprod[1] * vy[2] + cprod[2] * vy[1]); // d_{w[0]} A
    // second component
    dvA[1] = 2.0 * (cprod[0] * wy[2] - cprod[2] * wy[0]); // d_{v[1]} A
    dwA[1] = -2.0 * (cprod[0] *vy[2] - cprod[2] * vy[0]); //d_{w[1]} A
    // third component
    dvA[2] = -2.0 * (cprod[0] * wy[1] + cprod[1] * wy[0]); // d_{v[2]} A
    dwA[2] = 2.0 * (cprod[0] *  vy[1]  + cprod[1] * vy[0]); // d_{w[2]} A
    
    /*  double x = 2.0 * ((cprod[1] * wy[2] + cprod[2] * wy[1]) * vx[0]
                      -(cprod[1] * vy[2] + cprod[2] * vy[1]) * wx[0]
                      +(cprod[0] * wy[2] - cprod[2] * wy[0]) * vx[1]
                      -(cprod[0] * vy[2] - cprod[2] * vy[0]) * wx[1]
                      +(-cprod[0] * wy[1] - cprod[1] * wy[0]) * vx[2]
                      +( cprod[0] * vy[1] + cprod[1] * vy[0]) * wx[2]);*/
    
    double x = dvA[0]*vx[0] + dwA[0]*wx[0] + dvA[1]*vx[1]+ dwA[1]*wx[1]+dvA[2]*vx[2]+dwA[2]*wx[2];
    
    double d2A = d2phiarea(A/ARef) * x/(ARef*ARef) * alphaArea;
    
    
    
    rM[0] += d2A * (-dvA[0] - dwA[0]);
    r1[0] += d2A * ( dvA[0]         );
    r2[0] += d2A * (       dwA[0]   );
    
    
    rM[1] += d2A * (-dvA[1] - dwA[1]);
    r1[1] += d2A * ( dvA[1]      );
    r2[1] += d2A * (     + dwA[1]);
    
    rM[2] += d2A * (-dvA[2] - dwA[2]);
    r1[2] += d2A * ( dvA[2]      );
    r2[2] += d2A * (     + dwA[2]);
    
}

void S3D( double *S, double *dS, const double *yc,const double *ARef,const double *VRef, const double *alpha,const double *h,const double *m, int n, bool doDerivative){
    const int m1 = m[0], m2 = m[1], m3 = m[2];
    double yA[3], yB[3], yC[3], yD[3], yE[3], yF[3], yG[3], yH[3], yABDC[3], yBDHF[3], yABFE[3], yCAEG[3], yDCGH[3], yEFHG[3], yM[3];
    double rA[3], rB[3], rC[3], rD[3], rE[3], rF[3], rG[3], rH[3], rABDC[3], rBDHF[3], rABFE[3], rCAEG[3], rDCGH[3], rEFHG[3], rM[3];
    double sum = 0.0;
    int i1, i2, i3, r, comp,idn,off1,off2,off3;
    const double hd = h[0] * h[1] * h[2];
    const double alphaArea = alpha[1]*hd;
    const double alphaVolume = alpha[2]*hd;
    
    const int nn = (m1+1)*(m2+1)*(m3+1);
    const int inc_x = 1;
    const int inc_y = m1+1;
    const int inc_z = (m1+1)*(m2+1);
    
    const double *rptr;
    double *wptr;
    
    for (off1=0; off1<2; off1++) {
        for (off2=0; off2<2; off2++) {
            for (off3=0; off3<2; off3++) {
#pragma omp parallel for  default(shared) private(rptr,wptr,i1, i2, i3, r, comp,idn,  yA, yB, yC, yD, yE, yF, yG, yH, yABDC, yBDHF, yABFE, yCAEG, yDCGH, yEFHG, yM, rA, rB, rC, rD, rE, rF, rG, rH, rABDC, rBDHF, rABFE, rCAEG, rDCGH, rEFHG, rM) reduction(+ : sum)
                for (i3=off3; i3<m3; i3+=2){
                    for (i2=off2; i2<m2; i2+=2){
                        for (i1=off1; i1<m1; i1+=2){
                            idn = i1+ i2* inc_y+ i3*inc_z;
                            comp = 0;
                            for(r=0;r<3;r++){
                                rptr = yc+idn+comp;
                                // A
                                yA[r] = *rptr;
                                // B
                                rptr += inc_x;
                                yB[r] = *rptr;
                                // D
                                rptr += inc_y;
                                yD[r] = *rptr;
                                // C
                                rptr -= inc_x;
                                yC[r] = *rptr;
                                // G
                                rptr += inc_z;
                                yG[r] = *rptr;
                                // H
                                rptr += inc_x;
                                yH[r] = *rptr;
                                // F
                                rptr -= inc_y;
                                yF[r] = *rptr;
                                // E
                                rptr -= inc_x;
                                yE[r] = *rptr;
                                comp +=nn;
                                
                                yABDC[r] = 0.25*(yA[r]+yB[r]+yD[r]+yC[r]);
                                yBDHF[r] = 0.25*(yB[r]+yD[r]+yH[r]+yF[r]);
                                yABFE[r] = 0.25*(yA[r]+yB[r]+yF[r]+yE[r]);
                                yCAEG[r] = 0.25*(yC[r]+yA[r]+yE[r]+yG[r]);
                                yDCGH[r] = 0.25*(yD[r]+yC[r]+yG[r]+yH[r]);
                                yEFHG[r] = 0.25*(yE[r]+yF[r]+yH[r]+yG[r]);
                                yM[r]    = 0.5 *(yABDC[r]+yEFHG[r]);
                                if(doDerivative) {
                                    rA[r] = 0.0;
                                    rB[r] = 0.0;
                                    rC[r] = 0.0;
                                    rD[r] = 0.0;
                                    rE[r] = 0.0;
                                    rF[r] = 0.0;
                                    rG[r] = 0.0;
                                    rH[r] = 0.0;
                                    rABDC[r] = 0.0;
                                    rBDHF[r] = 0.0;
                                    rABFE[r] = 0.0;
                                    rCAEG[r] = 0.0;
                                    rDCGH[r] = 0.0;
                                    rEFHG[r] = 0.0;
                                    rM[r] = 0.0;
                                }
                                
                                
                            }
                            
                            if (alphaArea>0) {
                                sum += arTriangle3D(yABDC, yA, yB, ARef[0], rABDC, rA, rB, doDerivative,alphaArea);
                                sum += arTriangle3D(yABDC, yB, yD, ARef[1], rABDC, rB, rD, doDerivative,alphaArea);
                                sum += arTriangle3D(yABDC, yD, yC, ARef[2], rABDC, rD, rC, doDerivative,alphaArea);
                                sum += arTriangle3D(yABDC, yC, yA, ARef[3], rABDC, rC, rA, doDerivative,alphaArea);
                                
                                sum += arTriangle3D(yBDHF, yB, yD, ARef[4], rBDHF, rB, rD, doDerivative,alphaArea);
                                sum += arTriangle3D(yBDHF, yD, yH, ARef[5], rBDHF, rD, rH, doDerivative,alphaArea);
                                sum += arTriangle3D(yBDHF, yH, yF, ARef[6], rBDHF, rH, rF, doDerivative,alphaArea);
                                sum += arTriangle3D(yBDHF, yF, yB, ARef[7], rBDHF, rF, rB, doDerivative,alphaArea);
                                
                                sum += arTriangle3D(yABFE, yA, yB, ARef[8], rABFE, rA, rB, doDerivative,alphaArea);
                                sum += arTriangle3D(yABFE, yB, yF, ARef[9], rABFE, rB, rF, doDerivative,alphaArea);
                                sum += arTriangle3D(yABFE, yF, yE, ARef[10], rABFE, rF, rE, doDerivative,alphaArea);
                                sum += arTriangle3D(yABFE, yE, yA, ARef[11], rABFE, rE, rA, doDerivative,alphaArea);
                                
                                sum += arTriangle3D(yCAEG, yC, yA, ARef[12], rCAEG, rC, rA, doDerivative,alphaArea);
                                sum += arTriangle3D(yCAEG, yA, yE, ARef[13], rCAEG, rA, rE, doDerivative,alphaArea);
                                sum += arTriangle3D(yCAEG, yE, yG, ARef[14], rCAEG, rE, rG, doDerivative,alphaArea);
                                sum += arTriangle3D(yCAEG, yG, yC, ARef[15], rCAEG, rG, rC, doDerivative,alphaArea);
                                
                                sum += arTriangle3D(yDCGH, yD, yC, ARef[16], rDCGH, rD, rC, doDerivative,alphaArea);
                                sum += arTriangle3D(yDCGH, yC, yG, ARef[17], rDCGH, rC, rG, doDerivative,alphaArea);
                                sum += arTriangle3D(yDCGH, yG, yH, ARef[18], rDCGH, rG, rH, doDerivative,alphaArea);
                                sum += arTriangle3D(yDCGH, yH, yD, ARef[19], rDCGH, rH, rD, doDerivative,alphaArea);
                                
                                sum += arTriangle3D(yEFHG, yE, yF, ARef[20], rEFHG, rE, rF, doDerivative,alphaArea);
                                sum += arTriangle3D(yEFHG, yF, yH, ARef[21], rEFHG, rF, rH, doDerivative,alphaArea);
                                sum += arTriangle3D(yEFHG, yH, yG, ARef[22], rEFHG, rH, rG, doDerivative,alphaArea);
                                sum += arTriangle3D(yEFHG, yG, yE, ARef[23], rEFHG, rG, rE, doDerivative,alphaArea);
                            }
                            
                            if (alphaVolume>0) {
                                // do the work (private)
                                sum   +=  SvolTetra3D(yA, yB, yABDC, yM, VRef[0], doDerivative, rA, rB, rABDC, rM, -1.0,alphaVolume);
                                sum   +=  SvolTetra3D(yB, yD, yABDC, yM, VRef[1], doDerivative, rB, rD, rABDC, rM, -1.0,alphaVolume);
                                sum   +=  SvolTetra3D(yD, yC, yABDC, yM, VRef[2], doDerivative, rD, rC, rABDC, rM, -1.0,alphaVolume);
                                sum   +=  SvolTetra3D(yC, yA, yABDC, yM, VRef[3], doDerivative, rC, rA, rABDC, rM, -1.0,alphaVolume);
                                
                                sum   +=  SvolTetra3D(yB, yD, yBDHF, yM, VRef[4], doDerivative, rB, rD, rBDHF, rM, 1.0,alphaVolume);
                                sum   +=  SvolTetra3D(yD, yH, yBDHF, yM, VRef[5], doDerivative, rD, rH, rBDHF, rM, 1.0,alphaVolume);
                                sum   +=  SvolTetra3D(yH, yF, yBDHF, yM, VRef[6], doDerivative, rH, rF, rBDHF, rM, 1.0,alphaVolume);
                                sum   +=  SvolTetra3D(yF, yB, yBDHF, yM, VRef[7], doDerivative, rF, rB, rBDHF, rM, 1.0,alphaVolume);
                                
                                sum   +=  SvolTetra3D(yA, yB, yABFE, yM, VRef[8], doDerivative, rA, rB, rABFE, rM, 1.0,alphaVolume);
                                sum   +=  SvolTetra3D(yB, yF, yABFE, yM, VRef[9], doDerivative, rB, rF, rABFE, rM, 1.0,alphaVolume);
                                sum   +=  SvolTetra3D(yF, yE, yABFE, yM, VRef[10], doDerivative, rF, rE, rABFE, rM, 1.0,alphaVolume);
                                sum   +=  SvolTetra3D(yE, yA, yABFE, yM, VRef[11], doDerivative, rE, rA, rABFE, rM, 1.0,alphaVolume);
                                
                                sum   +=  SvolTetra3D(yC, yA, yCAEG, yM, VRef[12], doDerivative, rC, rA, rCAEG, rM, 1.0,alphaVolume);
                                sum   +=  SvolTetra3D(yA, yE, yCAEG, yM, VRef[13], doDerivative, rA, rE, rCAEG, rM, 1.0,alphaVolume);
                                sum   +=  SvolTetra3D(yE, yG, yCAEG, yM, VRef[14], doDerivative, rE, rG, rCAEG, rM, 1.0,alphaVolume);
                                sum   +=  SvolTetra3D(yG, yC, yCAEG, yM, VRef[15], doDerivative, rG, rC, rCAEG, rM, 1.0,alphaVolume);
                                
                                sum   +=  SvolTetra3D(yD, yC, yDCGH, yM, VRef[16], doDerivative, rD, rC, rDCGH, rM, 1.0,alphaVolume);
                                sum   +=  SvolTetra3D(yC, yG, yDCGH, yM, VRef[17], doDerivative, rC, rG, rDCGH, rM, 1.0,alphaVolume);
                                sum   +=  SvolTetra3D(yG, yH, yDCGH, yM, VRef[18], doDerivative, rG, rH, rDCGH, rM, 1.0,alphaVolume);
                                sum   +=  SvolTetra3D(yH, yD, yDCGH, yM, VRef[19], doDerivative, rH, rD, rDCGH, rM, 1.0,alphaVolume);
                                
                                sum   +=  SvolTetra3D(yE, yF, yEFHG, yM, VRef[20], doDerivative, rE, rF, rEFHG, rM, 1.0,alphaVolume);
                                sum   +=  SvolTetra3D(yF, yH, yEFHG, yM, VRef[21], doDerivative, rF, rH, rEFHG, rM, 1.0,alphaVolume);
                                sum   +=  SvolTetra3D(yH, yG, yEFHG, yM, VRef[22], doDerivative, rH, rG, rEFHG, rM, 1.0,alphaVolume);
                                sum   +=  SvolTetra3D(yG, yE, yEFHG, yM, VRef[23], doDerivative, rG, rE, rEFHG, rM, 1.0,alphaVolume);
                                
                            }
                            
                            // write into result vector (private --> shared)
                            if (doDerivative) {
                                comp = 0;
                                for (r=0; r < 3; r++){
                                    // scale vecs
                                    rABDC[r] *= 0.25;
                                    rBDHF[r] *= 0.25;
                                    rABFE[r] *= 0.25;
                                    rCAEG[r] *= 0.25;
                                    rDCGH[r] *= 0.25;
                                    rEFHG[r] *= 0.25;
                                    rM[r]    *= 0.125;
                                    // get pointer to A
                                    wptr = dS + idn + comp;
                                    *wptr +=  rA[r]+ rABDC[r]+rABFE[r]+rCAEG[r] + rM[r];
                                    // B
                                    wptr += inc_x;
                                    *wptr +=  rB[r]+ rABDC[r]+rBDHF[r]+rABFE[r]+rM[r];
                                    // D
                                    wptr += inc_y;
                                    *wptr += rD[r]+rABDC[r]+rBDHF[r]+rDCGH[r]+rM[r];
                                    // C
                                    wptr -= inc_x;
                                    *wptr += rC[r]+rABDC[r]+rCAEG[r]+rDCGH[r]+rM[r];
                                    // G
                                    wptr += inc_z;
                                    *wptr += rG[r]+rCAEG[r]+rDCGH[r]+rEFHG[r]+rM[r];
                                    // H
                                    wptr += inc_x;
                                    *wptr += rH[r]+rBDHF[r]+rDCGH[r]+rEFHG[r]+rM[r];
                                    // F
                                    wptr -= inc_y;
                                    *wptr += rF[r]+rBDHF[r]+rABFE[r]+rEFHG[r]+rM[r];
                                    // E
                                    wptr -= inc_x;
                                    *wptr += rE[r]+rABFE[r]+rCAEG[r]+rEFHG[r]+rM[r];
                                    comp +=nn;
                                }
                            }
                        }
                    }
                }
            }
        }
    }


*S = sum;
return;
}



void d2S3D(double *prod, const double *yc,const double *x,const double *ARef,const double *VRef,const  double *alpha,const double *h,const double *m, int n){
    const int m1 = m[0], m2 = m[1],m3 = m[2];
    double yA[3], yB[3], yC[3], yD[3], yE[3], yF[3], yG[3], yH[3], yABDC[3], yBDHF[3], yABFE[3], yCAEG[3], yDCGH[3], yEFHG[3], yM[3];
    double xA[3], xB[3], xC[3], xD[3], xE[3], xF[3], xG[3], xH[3], xABDC[3], xBDHF[3], xABFE[3], xCAEG[3], xDCGH[3], xEFHG[3], xM[3];
    double rA[3], rB[3], rC[3], rD[3], rE[3], rF[3], rG[3], rH[3], rABDC[3], rBDHF[3], rABFE[3], rCAEG[3], rDCGH[3], rEFHG[3], rM[3];
    
    int id, i1, i2, i3, r, comp,idn,off1,off2,off3;
    const double hd = h[0] * h[1] * h[2];
    const double alphaLength = alpha[0]*hd;
    const double alphaArea = alpha[1]*hd;
    const double alphaVolume = alpha[2]*hd;
    
    const int nn = (m1+1)*(m2+1)*(m3+1);
    const int inc_x = 1;
    const int inc_y = m1+1;
    const int inc_z = (m1+1)*(m2+1);
    
    const double *rptr;
    double *wptr;
    
    for (off1=0; off1<2; off1++) {
        for (off2=0; off2<2; off2++) {
            for (off3=0; off3<2; off3++) {
#pragma omp parallel for  default(shared) private(rptr,wptr,id, i1, i2, i3, r, comp,idn, yA, yB, yC, yD, yE, yF, yG, yH, yABDC, yBDHF, yABFE, yCAEG, yDCGH, yEFHG, yM, rA, rB, rC, rD, rE, rF, rG, rH, rABDC, rBDHF, rABFE, rCAEG, rDCGH, rEFHG, rM, xA, xB, xC, xD, xE, xF, xG, xH, xABDC, xBDHF, xABFE, xCAEG, xDCGH, xEFHG, xM)
                for (i3=off3; i3<m3; i3+=2){
                    for (i2=off2; i2<m2; i2+=2){
                        for (i1=off1; i1<m1; i1+=2){
                            idn = i1+ i2* inc_y+ i3*inc_z;
                            // get edges, face stg and cell center of yc (shared --> private)
                            comp = 0;
                            for(r=0;r<3;r++){
                                rptr = yc+idn+comp;
                                // A
                                yA[r] = *rptr;
                                // B
                                rptr += inc_x;
                                yB[r] = *rptr;
                                // D
                                rptr += inc_y;
                                yD[r] = *rptr;
                                // C
                                rptr -= inc_x;
                                yC[r] = *rptr;
                                // G
                                rptr += inc_z;
                                yG[r] = *rptr;
                                // H
                                rptr += inc_x;
                                yH[r] = *rptr;
                                // F
                                rptr -= inc_y;
                                yF[r] = *rptr;
                                // E
                                rptr -= inc_x;
                                yE[r] = *rptr;
                                
                                yABDC[r] = 0.25*(yA[r]+yB[r]+yD[r]+yC[r]);
                                yBDHF[r] = 0.25*(yB[r]+yD[r]+yH[r]+yF[r]);
                                yABFE[r] = 0.25*(yA[r]+yB[r]+yF[r]+yE[r]);
                                yCAEG[r] = 0.25*(yC[r]+yA[r]+yE[r]+yG[r]);
                                yDCGH[r] = 0.25*(yD[r]+yC[r]+yG[r]+yH[r]);
                                yEFHG[r] = 0.25*(yE[r]+yF[r]+yH[r]+yG[r]);
                                yM[r]    = 0.5 *(yABDC[r]+yEFHG[r]);
                                
                                rptr = x+idn+comp;
                                // A
                                xA[r] = *rptr;
                                // B
                                rptr += inc_x;
                                xB[r] = *rptr;
                                // D
                                rptr += inc_y;
                                xD[r] = *rptr;
                                // C
                                rptr -= inc_x;
                                xC[r] = *rptr;
                                // G
                                rptr += inc_z;
                                xG[r] = *rptr;
                                // H
                                rptr += inc_x;
                                xH[r] = *rptr;
                                // F
                                rptr -= inc_y;
                                xF[r] = *rptr;
                                // E
                                rptr -= inc_x;
                                xE[r] = *rptr;
                                
                                xABDC[r] = 0.25*(xA[r]+xB[r]+xD[r]+xC[r]);
                                xBDHF[r] = 0.25*(xB[r]+xD[r]+xH[r]+xF[r]);
                                xABFE[r] = 0.25*(xA[r]+xB[r]+xF[r]+xE[r]);
                                xCAEG[r] = 0.25*(xC[r]+xA[r]+xE[r]+xG[r]);
                                xDCGH[r] = 0.25*(xD[r]+xC[r]+xG[r]+xH[r]);
                                xEFHG[r] = 0.25*(xE[r]+xF[r]+xH[r]+xG[r]);
                                xM[r]    = 0.5 *(xABDC[r]+xEFHG[r]);
                                comp +=nn;
                                
                                rA[r] = 0.0;
                                rB[r] = 0.0;
                                rC[r] = 0.0;
                                rD[r] = 0.0;
                                rE[r] = 0.0;
                                rF[r] = 0.0;
                                rG[r] = 0.0;
                                rH[r] = 0.0;
                                rABDC[r] = 0.0;
                                rBDHF[r] = 0.0;
                                rABFE[r] = 0.0;
                                rCAEG[r] = 0.0;
                                rDCGH[r] = 0.0;
                                rEFHG[r] = 0.0;
                                rM[r] = 0.0;
                            }
                            
                            
                            if (alphaArea>0) {
                                d2SarTriangle3D(yA, yB, yABDC, xA, xB, xABDC, ARef[0], rA, rB, rABDC,alphaArea);
                                d2SarTriangle3D(yB, yD, yABDC, xB, xD, xABDC, ARef[1], rB, rD, rABDC,alphaArea);
                                d2SarTriangle3D(yD, yC, yABDC, xD, xC, xABDC, ARef[2], rD, rC, rABDC,alphaArea);
                                d2SarTriangle3D(yC, yA, yABDC, xC, xA, xABDC, ARef[3], rC, rA, rABDC,alphaArea);
                                
                                d2SarTriangle3D(yB, yD, yBDHF, xB, xD, xBDHF, ARef[4], rB, rD, rBDHF,alphaArea);
                                d2SarTriangle3D(yD, yH, yBDHF, xD, xH, xBDHF, ARef[5], rD, rH, rBDHF,alphaArea);
                                d2SarTriangle3D(yH, yF, yBDHF, xH, xF, xBDHF, ARef[6], rH, rF, rBDHF,alphaArea);
                                d2SarTriangle3D(yF, yB, yBDHF, xF, xB, xBDHF, ARef[7], rF, rB, rBDHF,alphaArea);
                                
                                d2SarTriangle3D(yA, yB, yABFE, xA, xB, xABFE, ARef[8], rA, rB, rABFE,alphaArea);
                                d2SarTriangle3D(yB, yF, yABFE, xB, xF, xABFE, ARef[9], rB, rF, rABFE,alphaArea);
                                d2SarTriangle3D(yF, yE, yABFE, xF, xE, xABFE, ARef[10], rF, rE, rABFE,alphaArea);
                                d2SarTriangle3D(yE, yA, yABFE, xE, xA, xABFE, ARef[11], rE, rA, rABFE,alphaArea);
                                
                                d2SarTriangle3D(yC, yA, yCAEG, xC, xA, xCAEG, ARef[12], rC, rA, rCAEG,alphaArea);
                                d2SarTriangle3D(yA, yE, yCAEG, xA, xE, xCAEG, ARef[13], rA, rE, rCAEG,alphaArea);
                                d2SarTriangle3D(yE, yG, yCAEG, xE, xG, xCAEG, ARef[14], rE, rG, rCAEG,alphaArea);
                                d2SarTriangle3D(yG, yC, yCAEG, xG, xC, xCAEG, ARef[15], rG, rC, rCAEG,alphaArea);
                                
                                d2SarTriangle3D(yD, yC, yDCGH, xD, xC, xDCGH, ARef[16], rD, rC, rDCGH,alphaArea);
                                d2SarTriangle3D(yC, yG, yDCGH, xC, xG, xDCGH, ARef[17], rC, rG, rDCGH,alphaArea);
                                d2SarTriangle3D(yG, yH, yDCGH, xG, xH, xDCGH, ARef[18], rG, rH, rDCGH,alphaArea);
                                d2SarTriangle3D(yH, yD, yDCGH, xH, xD, xDCGH, ARef[19], rH, rD, rDCGH,alphaArea);
                                
                                d2SarTriangle3D(yE, yF, yEFHG, xE, xF, xEFHG, ARef[20], rE, rF, rEFHG,alphaArea);
                                d2SarTriangle3D(yF, yH, yEFHG, xF, xH, xEFHG, ARef[21], rF, rH, rEFHG,alphaArea);
                                d2SarTriangle3D(yH, yG, yEFHG, xH, xG, xEFHG, ARef[22], rH, rG, rEFHG,alphaArea);
                                d2SarTriangle3D(yG, yE, yEFHG, xG, xE, xEFHG, ARef[23], rG, rE, rEFHG,alphaArea);
                                
                            }
                            
                            if (alphaVolume>0) {
                                // do the work (private)
                                d2SvolTetra3D(yA, yB, yABDC, yM, xA, xB, xABDC, xM, VRef[0], rA, rB, rABDC, rM, -1.0,alphaVolume);
                                d2SvolTetra3D(yB, yD, yABDC, yM, xB, xD, xABDC, xM, VRef[1], rB, rD, rABDC, rM, -1.0,alphaVolume);
                                d2SvolTetra3D(yD, yC, yABDC, yM, xD, xC, xABDC, xM, VRef[2], rD, rC, rABDC, rM, -1.0,alphaVolume);
                                d2SvolTetra3D(yC, yA, yABDC, yM, xC, xA, xABDC, xM, VRef[3], rC, rA, rABDC, rM, -1.0,alphaVolume);
                                
                                d2SvolTetra3D(yB, yD, yBDHF, yM, xB, xD, xBDHF, xM, VRef[4], rB, rD, rBDHF, rM, 1.0,alphaVolume);
                                d2SvolTetra3D(yD, yH, yBDHF, yM, xD, xH, xBDHF, xM, VRef[5], rD, rH, rBDHF, rM, 1.0,alphaVolume);
                                d2SvolTetra3D(yH, yF, yBDHF, yM, xH, xF, xBDHF, xM, VRef[6], rH, rF, rBDHF, rM, 1.0,alphaVolume);
                                d2SvolTetra3D(yF, yB, yBDHF, yM, xF, xB, xBDHF, xM, VRef[7], rF, rB, rBDHF, rM, 1.0,alphaVolume);
                                
                                d2SvolTetra3D(yA, yB, yABFE, yM, xA, xB, xABFE, xM, VRef[8], rA, rB, rABFE, rM, 1.0,alphaVolume);
                                d2SvolTetra3D(yB, yF, yABFE, yM, xB, xF, xABFE, xM, VRef[9], rB, rF, rABFE, rM, 1.0,alphaVolume);
                                d2SvolTetra3D(yF, yE, yABFE, yM, xF, xE, xABFE, xM, VRef[10], rF, rE, rABFE, rM, 1.0,alphaVolume);
                                d2SvolTetra3D(yE, yA, yABFE, yM, xE, xA, xABFE, xM, VRef[11], rE, rA, rABFE, rM, 1.0,alphaVolume);
                                
                                d2SvolTetra3D(yC, yA, yCAEG, yM, xC, xA, xCAEG, xM, VRef[12], rC, rA, rCAEG, rM, 1.0,alphaVolume);
                                d2SvolTetra3D(yA, yE, yCAEG, yM, xA, xE, xCAEG, xM, VRef[13], rA, rE, rCAEG, rM, 1.0,alphaVolume);
                                d2SvolTetra3D(yE, yG, yCAEG, yM, xE, xG, xCAEG, xM, VRef[14], rE, rG, rCAEG, rM, 1.0,alphaVolume);
                                d2SvolTetra3D(yG, yC, yCAEG, yM, xG, xC, xCAEG, xM, VRef[15], rG, rC, rCAEG, rM, 1.0,alphaVolume);
                                
                                d2SvolTetra3D(yD, yC, yDCGH, yM, xD, xC, xDCGH, xM, VRef[16], rD, rC, rDCGH, rM, 1.0,alphaVolume);
                                d2SvolTetra3D(yC, yG, yDCGH, yM, xC, xG, xDCGH, xM, VRef[17], rC, rG, rDCGH, rM, 1.0,alphaVolume);
                                d2SvolTetra3D(yG, yH, yDCGH, yM, xG, xH, xDCGH, xM, VRef[18], rG, rH, rDCGH, rM, 1.0,alphaVolume);
                                d2SvolTetra3D(yH, yD, yDCGH, yM, xH, xD, xDCGH, xM, VRef[19], rH, rD, rDCGH, rM, 1.0,alphaVolume);
                                
                                d2SvolTetra3D(yE, yF, yEFHG, yM, xE, xF, xEFHG, xM, VRef[20], rE, rF, rEFHG, rM, 1.0,alphaVolume);
                                d2SvolTetra3D(yF, yH, yEFHG, yM, xF, xH, xEFHG, xM, VRef[21], rF, rH, rEFHG, rM, 1.0,alphaVolume);
                                d2SvolTetra3D(yH, yG, yEFHG, yM, xH, xG, xEFHG, xM, VRef[22], rH, rG, rEFHG, rM, 1.0,alphaVolume);
                                d2SvolTetra3D(yG, yE, yEFHG, yM, xG, xE, xEFHG, xM, VRef[23], rG, rE, rEFHG, rM, 1.0,alphaVolume);
                                
                            }
                            
                            // write into result vector (private --> shared)
                            comp = 0;
                            for (r=0; r < 3; r++){
                                // scale vecs
                                rABDC[r] *= 0.25;
                                rBDHF[r] *= 0.25;
                                rABFE[r] *= 0.25;
                                rCAEG[r] *= 0.25;
                                rDCGH[r] *= 0.25;
                                rEFHG[r] *= 0.25;
                                rM[r]    *= 0.125;
                                // get pointer to A
                                wptr = prod + idn + comp;
                                *wptr +=  rA[r]+ rABDC[r]+rABFE[r]+rCAEG[r] + rM[r];
                                // B
                                wptr += inc_x;
                                *wptr +=  rB[r]+ rABDC[r]+rBDHF[r]+rABFE[r]+rM[r];
                                // D
                                wptr += inc_y;
                                *wptr += rD[r]+rABDC[r]+rBDHF[r]+rDCGH[r]+rM[r];
                                // C
                                wptr -= inc_x;
                                *wptr += rC[r]+rABDC[r]+rCAEG[r]+rDCGH[r]+rM[r];
                                // G
                                wptr += inc_z;
                                *wptr += rG[r]+rCAEG[r]+rDCGH[r]+rEFHG[r]+rM[r];
                                // H
                                wptr += inc_x;
                                *wptr += rH[r]+rBDHF[r]+rDCGH[r]+rEFHG[r]+rM[r];
                                // F
                                wptr -= inc_y;
                                *wptr += rF[r]+rBDHF[r]+rABFE[r]+rEFHG[r]+rM[r];
                                // E
                                wptr -= inc_x;
                                *wptr += rE[r]+rABFE[r]+rCAEG[r]+rEFHG[r]+rM[r];
                                comp +=nn;
                            }
                        }
                    }
                }
            }
        }
    }
  

}
inline void d2SdiagTetra(double *y1, double *y2, double *yF, double *yM,const double ARef,
        const double VRef, double *r1, double *r2, double *r3, double *r4,
        double *r5, double *r6, double *r7, double *r8,const double inv , double alphaArea,
        double alphaVolume){
    double dVx, vol,vol2, d2S, cof[9],x,A;
    double vy[3], wy[3];
    int i, j,k;
    
    // things we need anyway
    vol = inv* det3D(y1, y2, yF, yM);
    cofactor3D(cof, y1, y2, yF, yM,1.0);
    
    if (alphaArea>0) {
        double vy[3], wy[3], cprod[3], dvA[3], dwA[3];
        
        vecminus(vy, y1, yF, 3);  // vy = y1 - yM;
        vecminus(wy, y2, yF, 3);  // wx = y2 - yM;
        
        // compute cross-product for y
        cprod[0] =  vy[1] * wy[2] - vy[2] * wy[1];
        cprod[1] =  vy[0] * wy[2] - vy[2] * wy[0];
        cprod[2] =  vy[0] * wy[1] - vy[1] * wy[0];
        
        A= cprod[0]*cprod[0] + cprod[1]*cprod[1] + cprod[2]*cprod[2]; // area of triangle
        
        dvA[0] = 2.0 * (cprod[1] * wy[2] + cprod[2] * wy[1]); // d_{v[0]} A
        dvA[1] = 2.0 * (cprod[0] * wy[2] - cprod[2] * wy[0]); // d_{v[1]} A
        dvA[2] = 2.0 * (-cprod[0] * wy[1] - cprod[1] * wy[0]);// d_{v[2]} A
        
        dwA[0] = -2.0 * (cprod[1] * vy[2] + cprod[2] * vy[1]); // d_{w[0]} A
        dwA[1] = 2.0 * (-cprod[0] *vy[2] + cprod[2] * vy[0]); //d_{w[1]} A
        dwA[2] = 2.0 * (cprod[0] *  vy[1]  + cprod[1] * vy[0]); // d_{w[2]} A
        
        // now the cases
        double dAx,
                
                d2A= d2phiarea(A/ARef)/(ARef*ARef);
        
        // Two edges are part of the triangle, namely x1 and x2
        /* template
         * dAx =   (dvA[0] * vx[0] + dwA[0] * wx[0])
         * + (dvA[1] * vx[1] + dwA[1] * wx[1])
         * + (dvA[2] * vx[2] + dwA[2] * wx[2]);
         *
         */
        // 1) x1[0] = 1.0 --> xF[0] = 0.25 ==> vx[0] = 0.75 and wx[0]=-0.25 ==> write in r1[0]
        dAx =     (dvA[0] * 0.75 + dwA[0] * (-0.25));
        r1[0] += alphaArea * dAx * d2A * dAx;
        
        // 2) x1[1] = 1.0 --> xF[1] = 0.25 ==> vx[1] = 0.75 and wx[1]=-0.25 ==> write in r1[1]
        dAx =    (dvA[1] * 0.75 + dwA[1] * (-0.25));
        r1[1] += alphaArea * dAx * d2A * dAx;
        
        // 3) x1[2] = 1.0 --> xF[2] = 0.25 ==> vx[2] = 0.75 and wx[2]=-0.25 ==> write in r1[2]
        dAx =    (dvA[2] * 0.75 + dwA[2] * (-0.25));
        r1[2] += alphaArea * dAx * d2A * dAx;
        
        // 4) x2[0] = 1.0 --> xF[0] = 0.25 ==> vx[0] = -0.25 and wx[0]= 0.75 ==> write in r2[0]
        dAx =   (dvA[0] *(-0.25) + dwA[0] * 0.75);
        r2[0] += alphaArea * dAx * d2A * dAx;
        
        // 5) x2[1] = 1.0 --> xF[1] = 0.25 ==> vx[1] = -0.25 and wx[1]= 0.75 ==> write in r2[1]
        dAx = (dvA[1] * (-0.25) + dwA[1] * 0.75);
        r2[1] += alphaArea * dAx * d2A * dAx;
        
        // 6) x2[2] = 1.0 --> xF[2] = 0.25 ==> vx[2] = -0.25 and wx[2]= 0.75 ==> write in r2[2]
        dAx = (dvA[2] * (-0.25) + dwA[2] * 0.75);
        r2[2] += alphaArea * dAx * d2A * dAx;
        
        // ----- now the edges not part of the triangle (x3 and x4) interacting via xF
        // 7) x3[0] = 1.0 --> xF[0]= 0.25 ==> vx[0] = -0.25 and wx[0] = -0.25 ==> write in r3[0]
        dAx =   (dvA[0] * (-0.25) + dwA[0] * (-0.25));
        dAx = alphaArea * dAx * d2A * dAx;
        r3[0] += dAx;
        r4[0] += dAx;
        // 8) x3[1] = 1.0 --> xF[1]= 0.25 ==> vx[1] = -0.25 and wx[1] = -0.25 ==> write in r3[1]
        dAx =   (dvA[1] * (-0.25) + dwA[1] * (-0.25));
        dAx = alphaArea * dAx * d2A * dAx;
        r3[1] += dAx;
        r4[1] += dAx;
        // 9) x3[2] = 1.0 --> xF[2]= 0.25 ==> vx[2] = -0.25 and wx[2] = -0.25 ==> write in r3[2]
        dAx =   (dvA[2] * (-0.25) + dwA[2] * (-0.25));
        dAx = alphaArea * dAx * d2A * dAx;
        r3[2] += dAx;
        r4[2] += dAx;
        
    }
    if (alphaVolume>0) {
        // template for
        // dVx =   cof[0] * (x1[0]-xM[0]) + cof[3] * (x1[1]-xM[1]) + cof[6] * (x1[2]-xM[2])
        //       + cof[1] * (x2[0]-xM[0]) + cof[4] * (x2[1]-xM[1]) + cof[7] * (x2[2]-xM[2])
        //       + cof[2] * (xF[0]-xM[0]) + cof[5] * (xF[1]-xM[1]) + cof[8] * (xF[2]-xM[2]);
        d2S = alphaVolume* d2phicaley(vol/VRef)/(VRef*VRef);
        
        // now there are 6 cases directly on this tetrahedra
        // 1) x1[0] = 1 --> xF[0] = 0.25 and xM[0] = 0.125 ==> write to r1[0]
        dVx = cof[0] * (1-0.125) + cof[1] * (-0.125) + cof[2] * (0.25-0.125);
        r1[0] += dVx * d2S * dVx;
        // 2) x1[1] = 1 --> xF[1] = 0.25 and xM[1] = 0.125 ==> write to r1[1]
        dVx = cof[3] * (1.0-0.125) + cof[4] * (-0.125) + cof[5] * (0.25-0.125);
        r1[1] += dVx * d2S * dVx;
        // 3) x1[2] = 1 --> xF[2] = 0.25 and xM[2] = 0.125 ==> write to r1[2]
        dVx = cof[6] * (1.0-0.125) + cof[7] * (-0.125) + cof[8] * (0.25-0.125);
        r1[2] += dVx * d2S * dVx;
        
        // 4) x2[0] = 1 --> xF[0] = 0.25 and xM[0] = 0.125 ==> write to r2[0]
        dVx = cof[0] * (-0.125) + cof[1] * (1.0-0.125) + cof[2] * (0.25-0.125);
        r2[0] += dVx * d2S * dVx;
        // 5) x2[1] = 1 --> xF[1] = 0.25 and xM[1] = 0.125 ==> write to r2[1]
        dVx = cof[3] * (-0.125) + cof[4] * (1.0-0.125) + cof[5] * (0.25-0.125);
        r2[1] += dVx * d2S * dVx;
        // 6) x2[2] = 1 --> xF[2] = 0.25 and xM[2] = 0.125 ==> write to r2[2]
        dVx = cof[6] * (-0.125) + cof[7] * (1.0-0.125) + cof[8] * (0.25-0.125);
        r2[2] += dVx * d2S * dVx;
        
        // there are also 6 cases not on this terahedra, but on the same face
        // 7) x3[0] = 1 --> xF[0] = 0.25 and xM[0] = 0.125 ==> write to r3[0]
        dVx = cof[0] * (-0.125) + cof[1] * (-0.125) + cof[2] * (0.25-0.125);
        dVx = dVx * d2S * dVx;
        r3[0] += dVx;
        r4[0] += dVx;
        // 8) x3[1] = 1 --> xF[1] = 0.25 and xM[1] = 0.125 ==> write to r3[1]
        dVx = cof[3] * (-0.125) + cof[4] * (-0.125) + cof[5] * (0.25-0.125);
        dVx = dVx * d2S * dVx;
        r3[1] += dVx;
        r4[1] += dVx;
        // 9) x3[2] = 1 --> xF[2] = 0.25 and xM[2] = 0.125 ==> write to r3[2]
        dVx = cof[6] * (-0.125) + cof[7] * (-0.125) + cof[8] * (0.25-0.125);
        dVx = dVx * d2S * dVx;
        r3[2] += dVx;
        r4[2] += dVx;
        
        
        // finally there are  12 cases left, all on the opposite surface, so interaction just via the cell-center
        // 13) x5[0] = 1 -->  xM[0] = 0.125 ==> write to r5[0]
        dVx = cof[0] * (-0.125) + cof[1] * (-0.125) + cof[2] * (-0.125);
        dVx = dVx * d2S * dVx;
        r5[0] += dVx;
        r6[0] += dVx;
        r7[0] += dVx;
        r8[0] += dVx;
        // 14) x5[1] = 1 -->  xM[1] = 0.125 ==> write to r5[1]
        dVx = cof[3] * (-0.125) + cof[4] * (-0.125) + cof[5] * (-0.125);
        dVx = dVx * d2S * dVx;
        r5[1] += dVx;
        r6[1] += dVx;
        r7[1] += dVx;
        r8[1] += dVx;
        
        // 15) x5[2] = 1 -->  xM[2] = 0.125 ==> write to r5[2]
        dVx = cof[6] * (-0.125) + cof[7] * (-0.125) + cof[8] * (-0.125);
        dVx = dVx * d2S * dVx;
        r5[2] += dVx;
        r6[2] += dVx;
        r7[2] += dVx;
        r8[2] += dVx;
        
    }
    
}


void d2S3Ddiag(double *prod, const double *yc,const double *ARef,const double *VRef,const  double *alpha,const double *h,const double *m, int n){
    const int m1 = m[0], m2 = m[1],m3 = m[2];
    double yA[3], yB[3], yC[3], yD[3], yE[3], yF[3], yG[3], yH[3], yABDC[3], yBDHF[3], yABFE[3], yCAEG[3], yDCGH[3], yEFHG[3], yM[3];
    double rA[3], rB[3], rC[3], rD[3], rE[3], rF[3], rG[3], rH[3];
    
    int id, i1, i2, i3, r, comp,idn,off1,off2,off3;
    const double hd = h[0] * h[1] * h[2];
    const double alphaLength = alpha[0]*hd;
    const double alphaArea = alpha[1]*hd;
    const double alphaVolume = alpha[2]*hd;
    
    const int nn = (m1+1)*(m2+1)*(m3+1);
    const int inc_x = 1;
    const int inc_y = m1+1;
    const int inc_z = (m1+1)*(m2+1);
    
    const double *rptr;
    double *wptr;
    
    for (off1=0; off1<2; off1++) {
        for (off2=0; off2<2; off2++) {
            for (off3=0; off3<2; off3++) {
#pragma omp parallel for  default(shared) private(rptr,wptr,id, i1, i2, i3, r, comp,idn, yA, yB, yC, yD, yE, yF, yG, yH, yABDC, yBDHF, yABFE, yCAEG, yDCGH, yEFHG, yM, rA, rB, rC, rD, rE, rF, rG, rH)
                for (i3=off3; i3<m3; i3+=2){
                    for (i2=off2; i2<m2; i2+=2){
                        for (i1=off1; i1<m1; i1+=2){
                            idn = i1+ i2* inc_y+ i3*inc_z;
                            // get edges, face stg and cell center of yc (shared --> private)
                            comp = 0;
                            for(r=0;r<3;r++){
                                rptr = yc+idn+comp;
                                // A
                                yA[r] = *rptr;
                                // B
                                rptr += inc_x;
                                yB[r] = *rptr;
                                // D
                                rptr += inc_y;
                                yD[r] = *rptr;
                                // C
                                rptr -= inc_x;
                                yC[r] = *rptr;
                                // G
                                rptr += inc_z;
                                yG[r] = *rptr;
                                // H
                                rptr += inc_x;
                                yH[r] = *rptr;
                                // F
                                rptr -= inc_y;
                                yF[r] = *rptr;
                                // E
                                rptr -= inc_x;
                                yE[r] = *rptr;
                                
                                yABDC[r] = 0.25*(yA[r]+yB[r]+yD[r]+yC[r]);
                                yBDHF[r] = 0.25*(yB[r]+yD[r]+yH[r]+yF[r]);
                                yABFE[r] = 0.25*(yA[r]+yB[r]+yF[r]+yE[r]);
                                yCAEG[r] = 0.25*(yC[r]+yA[r]+yE[r]+yG[r]);
                                yDCGH[r] = 0.25*(yD[r]+yC[r]+yG[r]+yH[r]);
                                yEFHG[r] = 0.25*(yE[r]+yF[r]+yH[r]+yG[r]);
                                yM[r]    = 0.5 *(yABDC[r]+yEFHG[r]);
                                
                                comp +=nn;
                                
                                rA[r] = 0.0;
                                rB[r] = 0.0;
                                rC[r] = 0.0;
                                rD[r] = 0.0;
                                rE[r] = 0.0;
                                rF[r] = 0.0;
                                rG[r] = 0.0;
                                rH[r] = 0.0;
                            }
                            
                            
                            if ((alphaArea>0) || (alphaVolume>0)) {
                                d2SdiagTetra(yA, yB, yABDC, yM, ARef[0], VRef[0], rA, rB, rD, rC, rE, rF, rH, rG, -1.0, alphaArea, alphaVolume);
                                d2SdiagTetra(yB, yD, yABDC, yM, ARef[1], VRef[1], rB, rD, rC, rA, rF, rH, rG, rE, -1.0, alphaArea, alphaVolume);
                                d2SdiagTetra(yD, yC, yABDC, yM, ARef[2], VRef[2], rD, rC, rA, rB, rH, rG, rE, rF, -1.0, alphaArea, alphaVolume);
                                d2SdiagTetra(yC, yA, yABDC, yM, ARef[3], VRef[3], rC, rA, rB, rD, rG, rE, rF, rH, -1.0, alphaArea, alphaVolume);
                                
                                d2SdiagTetra(yB, yD, yBDHF, yM, ARef[4], VRef[4], rB, rD, rH, rF, rA, rC, rG, rE, 1.0,alphaArea, alphaVolume);
                                d2SdiagTetra(yD, yH, yBDHF, yM, ARef[5], VRef[5], rD, rH, rF, rB, rC, rG, rE, rA, 1.0,alphaArea, alphaVolume);
                                d2SdiagTetra(yH, yF, yBDHF, yM, ARef[6], VRef[6], rH, rF, rB, rD, rG, rE, rA, rC, 1.0,alphaArea, alphaVolume);
                                d2SdiagTetra(yF, yB, yBDHF, yM, ARef[7], VRef[7], rF, rB, rD, rH, rE, rA, rC, rG, 1.0,alphaArea, alphaVolume);
                                
                                d2SdiagTetra(yA, yB, yABFE, yM, ARef[8], VRef[8],   rA, rB, rF, rE, rC, rD, rH, rG, 1.0,alphaArea, alphaVolume);
                                d2SdiagTetra(yB, yF, yABFE, yM, ARef[9], VRef[9],   rB, rF, rE, rA, rD, rH, rG, rC, 1.0,alphaArea, alphaVolume);
                                d2SdiagTetra(yF, yE, yABFE, yM, ARef[10], VRef[10], rF, rE, rA, rB, rH, rG, rC, rD, 1.0,alphaArea, alphaVolume);
                                d2SdiagTetra(yE, yA, yABFE, yM, ARef[11], VRef[11], rE, rA, rB, rF, rG, rC, rD, rH, 1.0,alphaArea, alphaVolume);
                                
                                d2SdiagTetra(yC, yA, yCAEG, yM, ARef[12], VRef[12], rC, rA, rE, rG, rD, rB, rF, rH, 1.0,alphaArea, alphaVolume);
                                d2SdiagTetra(yA, yE, yCAEG, yM, ARef[13], VRef[13], rA, rE, rG, rC, rB, rF, rH, rD, 1.0,alphaArea, alphaVolume);
                                d2SdiagTetra(yE, yG, yCAEG, yM, ARef[14], VRef[14], rE, rG, rC, rA, rF, rH, rD, rB, 1.0,alphaArea, alphaVolume);
                                d2SdiagTetra(yG, yC, yCAEG, yM, ARef[15], VRef[15], rG, rC, rA, rE, rH, rD, rB, rF, 1.0,alphaArea, alphaVolume);
                                
                                d2SdiagTetra(yD, yC, yDCGH, yM, ARef[16], VRef[16], rD, rC, rG, rH, rB, rA, rE, rF, 1.0,alphaArea, alphaVolume);
                                d2SdiagTetra(yC, yG, yDCGH, yM, ARef[17], VRef[17], rC, rG, rH, rD, rA, rE, rF, rB, 1.0,alphaArea, alphaVolume);
                                d2SdiagTetra(yG, yH, yDCGH, yM, ARef[18], VRef[18], rG, rH, rD, rC, rE, rF, rB, rA, 1.0,alphaArea, alphaVolume);
                                d2SdiagTetra(yH, yD, yDCGH, yM, ARef[19], VRef[19], rH, rD, rC, rG, rF, rB, rA, rE, 1.0,alphaArea, alphaVolume);
                                
                                d2SdiagTetra(yE, yF, yEFHG, yM, ARef[20], VRef[20], rE, rF, rH, rG, rA, rB, rD, rC, 1.0,alphaArea, alphaVolume);
                                d2SdiagTetra(yF, yH, yEFHG, yM, ARef[21], VRef[21], rF, rH, rG, rE, rB, rD, rC, rA, 1.0,alphaArea, alphaVolume);
                                d2SdiagTetra(yH, yG, yEFHG, yM, ARef[22], VRef[22], rH, rG, rE, rF, rD, rC, rA, rB, 1.0,alphaArea, alphaVolume);
                                d2SdiagTetra(yG, yE, yEFHG, yM, ARef[23], VRef[23], rG, rE, rF, rH, rC, rA, rB, rD, 1.0,alphaArea, alphaVolume);
                                
                            }
                            // write into result vector (private --> shared)
                            comp = 0;
                            for (r=0; r < 3; r++){
                                // scale vecs
                                // get pointer to A
                                wptr = prod + idn + comp;
                                *wptr +=  rA[r];
                                // B
                                wptr += inc_x;
                                *wptr +=  rB[r];
                                // D
                                wptr += inc_y;
                                *wptr += rD[r];
                                // C
                                wptr -= inc_x;
                                *wptr += rC[r];
                                // G
                                wptr += inc_z;
                                *wptr += rG[r];
                                // H
                                wptr += inc_x;
                                *wptr += rH[r];
                                // F
                                wptr -= inc_y;
                                *wptr += rF[r];
                                // E
                                wptr -= inc_x;
                                *wptr += rE[r];
                                comp +=nn;
                            }
                        }
                    }
                }
            }
        }
    }
   

}

void area3D( double *A,const double *yc, const double *m, int n){
    const int m1 = m[0], m2 = m[1],m3 = m[2];
    double yA[3], yB[3], yC[3], yD[3], yE[3], yF[3], yG[3], yH[3], yABDC[3], yBDHF[3], yABFE[3], yCAEG[3], yDCGH[3], yEFHG[3], yM[3];
    double area[24];
    
    int id, i1, i2, i3, r, comp,idn;
    
    const int nn = (m1+1)*(m2+1)*(m3+1);
    const int inc_x = 1;
    const int inc_y = m1+1;
    const int inc_z = (m1+1)*(m2+1);
    
    const double *rptr;
#pragma omp parallel for  default(shared) private(rptr,id, i1, i2, i3, r, comp,idn,area, yA, yB, yC, yD, yE, yF, yG, yH, yABDC, yBDHF, yABFE, yCAEG, yDCGH, yEFHG)
            for (i3=0; i3<m3; i3++){
                for (i2=0; i2<m2; i2++){
                    for (i1=0; i1<m1; i1++){
                        idn = i1+ i2* inc_y+ i3*inc_z;
                        // get edges, face stg and cell center of yc (shared --> private)
                        comp = 0;
                for(r=0;r<3;r++){
                    rptr = yc+idn+comp;
                    // A
                    yA[r] = *rptr;
                    // B
                    rptr += inc_x;
                    yB[r] = *rptr;
                    // D
                    rptr += inc_y;
                    yD[r] = *rptr;
                    // C
                    rptr -= inc_x;
                    yC[r] = *rptr;
                    // G
                    rptr += inc_z;
                    yG[r] = *rptr;
                    // H
                    rptr += inc_x;
                    yH[r] = *rptr;
                    // F
                    rptr -= inc_y;
                    yF[r] = *rptr;
                    // E
                    rptr -= inc_x;
                    yE[r] = *rptr;
                    
                    yABDC[r] = 0.25*(yA[r]+yB[r]+yD[r]+yC[r]);
                    yBDHF[r] = 0.25*(yB[r]+yD[r]+yH[r]+yF[r]);
                    yABFE[r] = 0.25*(yA[r]+yB[r]+yF[r]+yE[r]);
                    yCAEG[r] = 0.25*(yC[r]+yA[r]+yE[r]+yG[r]);
                    yDCGH[r] = 0.25*(yD[r]+yC[r]+yG[r]+yH[r]);
                    yEFHG[r] = 0.25*(yE[r]+yF[r]+yH[r]+yG[r]);
                    yM[r]    = 0.5 *(yABDC[r]+yEFHG[r]);
                    
                    comp +=nn;
                    
                }
                
                
                
                // do the work (private)
                area[0] = areaTriangle3D(yABDC, yA, yB);
                area[1] = areaTriangle3D(yABDC, yB, yD);
                area[2] = areaTriangle3D(yABDC, yD, yC);
                area[3] = areaTriangle3D(yABDC, yC, yA);
                
                area[4] = areaTriangle3D(yBDHF, yB, yD);
                area[5] = areaTriangle3D(yBDHF, yD, yH);
                area[6] = areaTriangle3D(yBDHF, yH, yF);
                area[7] = areaTriangle3D(yBDHF, yF, yB);
                
                area[8] = areaTriangle3D(yABFE, yA, yB);
                area[9] = areaTriangle3D(yABFE, yB, yF);
                area[10] = areaTriangle3D(yABFE, yF, yE);
                area[11] = areaTriangle3D(yABFE, yE, yA);
                
                area[12] = areaTriangle3D(yCAEG, yC, yA);
                area[13] = areaTriangle3D(yCAEG, yA, yE);
                area[14] = areaTriangle3D(yCAEG, yE, yG);
                area[15] = areaTriangle3D(yCAEG, yG, yC);
                
                area[16] = areaTriangle3D(yDCGH, yD, yC);
                area[17] = areaTriangle3D(yDCGH, yC, yG);
                area[18] = areaTriangle3D(yDCGH, yG, yH);
                area[19] = areaTriangle3D(yDCGH, yH, yD);
                
                area[20] = areaTriangle3D(yEFHG, yE, yF);
                area[21] = areaTriangle3D(yEFHG, yF, yH);
                area[22] = areaTriangle3D(yEFHG, yH, yG);
                area[23] = areaTriangle3D(yEFHG, yG, yE);
                
                // write into result vector (private --> shared)
                id =  i1+m1 *(i2+i3*m2);
                for(r=0;r<24;r++){
                    A[id] = area[r];
                    id+=n;
                }
            }
        }
    }
    return;
    
}
# pragma mark -
# pragma mark VAMPIRE
void VAMPIREdiagTriangle(double *y1, double *y2, double *yM, double *r1, double *r2, double *r3, double *r4,
        double Ti){
    double dVx, cof[4];
    cofactor2D(cof, y1, y2, yM,1.0);
    
    //dVx = (cof[0] * (x1[0]-xM[0]) + cof[2] * (x1[1]-xM[1]) + cof[1] * (x2[0]-xM[0]) + cof[3] * (x2[1]-xM[1]));
    
    // now we have 4 cases inside this triangle
    //  1) x1[0] = 1 --> xM[0] = 1/4 --> write to r1[0]
    dVx = cof[0] * (0.75) + cof[1] * (-0.25);
    r1[0] += Ti * dVx;
    
    //  2) x1[1] = 1 --> xM[1] = 1/4 --> write to r1[1]
    dVx = cof[2] * (0.75) + cof[3] * (-0.25);
    r1[1] += Ti * dVx;
    
    //  3) x2[0] = 1 --> xM[0] = 1/4 --> write to r2[0]
    dVx = cof[0] * (-0.25) +  cof[1] * (0.75);
    r2[0] += Ti * dVx;
    
    //  4) x2[1] = 1 --> xM[1] = 1/4 --> write to r2[1]
    dVx =  cof[2] * (-0.25) + cof[3] * (0.75) ;
    r2[1] += Ti * dVx;;
    
    
    // and we have 2 cases outside this triangle
    //  1) x3[0] = 1 or x4[0] = 1--> xM[0] = 1/4 --> write to r3[0] and r4[0]
    dVx = (cof[0] + cof[1]) * (-0.25);
    r3[0] += Ti * dVx;
    r4[0] += Ti * dVx;
    
    //  2) x3[1] = 1 or x4[1] --> xM[1] = 1/4 --> write to r3[1] and r4[1]
    dVx = (cof[2]  + cof[3]) * (-0.25);
    r3[1] += Ti * dVx ;
    r4[1] += Ti * dVx;
    
}
void VAMPIREdiag2D(double *diag, const double *yc, const double *m, const double *Jac,
        const double *dT, const double *Tc,const double *h){
    int idn, i1, i2, r, comp,id,off1,off2;
    double yA[2], yB[2], yC[2], yD[2], yM[2];
    double dTi[2],Jaci,Ti;
    double rA[2], rB[2], rC[2], rD[2];
    
    const int m1 = m[0], m2 = m[1];
    const int nn = (m1+1)*(m2+1);
    const int n = (m1)*(m2);
    const double hd = h[0] * h[1];
    
    // increments along x and y
    const int inc_x=1;
    const int inc_y= m1+1;
    // pointers for reading / writing
    const double *ptr;
    double *wptr;
    for (off1=0; off1<2; off1++) {
        for (off2=0; off2<2; off2++) {
#pragma omp parallel private(wptr,ptr,idn, i1, i2, r, comp,dTi,Jaci,Ti,yA,yB,yC,yD,yM,rA,rB,rC,rD)
            {
#pragma omp for
                for (i2=off2; i2<m2; i2+=2){
                    for (i1=off1; i1<m1; i1+=2){
                        idn = i1+i2*(m1+1);
                        id = i1+i2*(m1);
                        comp = 0;
                        // get image and derivatives
                        Ti = Tc[id];
                        dTi[0] = dT[id];
                        dTi[1] = dT[id+n];
                        Jaci = Jac[id];
                        for (r=0; r<2; r++) {
                            // get y at edges
                            ptr = yc + idn+comp;
                            yA[r] = *ptr; // yc[idn           +comp];
                            ptr = ptr+inc_x;
                            yB[r] = *ptr; //yc[idn+1         +comp];
                            ptr = ptr+inc_y;
                            yD[r] = *ptr; //yc[idn+1+(m1+1)  +comp];
                            ptr = ptr-inc_x;
                            yC[r] = *ptr; //* yc[idn  +(m1+1)  +comp];
                            yM[r] = 0.25*(yA[r]+yB[r]+yC[r]+yD[r]);
                            
                            // init private result vectors (private)
                            rA[r] = Jaci * 0.25 * dTi[r];
                            rB[r] = Jaci * 0.25 * dTi[r];
                            rC[r] = Jaci * 0.25 * dTi[r];
                            rD[r] = Jaci * 0.25 * dTi[r];
                            comp+=nn;
                            
                        }
                        // do the work (private)
                        VAMPIREdiagTriangle(yA, yB, yM, rA, rB,rD,rC,Ti/hd);
                        VAMPIREdiagTriangle(yB, yD, yM, rB, rD,rC,rA,Ti/hd);
                        VAMPIREdiagTriangle(yD, yC, yM, rD, rC,rA,rB,Ti/hd);
                        VAMPIREdiagTriangle(yC, yA, yM, rC, rA,rB,rD,Ti/hd);
                        
                        // write into result vector (private --> shared)
                        comp = 0;
                        for (r=0; r < 2; r++){
                            wptr = diag + idn + comp;
                            //dS[idn     +comp] += rA[r]+0.25*rM[r];
                            *wptr += hd * rA[r] * rA[r];
                            wptr += inc_x;
                            //dS[idn +1   +comp] += rB[r]+rM[r];
                            *wptr += hd *rB[r]*rB[r];
                            wptr +=  inc_y;
                            //dS[idn    +(m1+1)   +comp] += rC[r]+rM[r];
                            *wptr += hd *rD[r]*rD[r];
                            wptr -= inc_x;
                            *wptr += hd *rC[r]*rC[r];
                            //dS[idn +1 +(m1+1)   +comp] += rD[r]+rM[r];
                            comp+=nn;
                        }
                    }
                }
            }
        }
    }


    
    
}
inline void VAMPIREdiagTetra(double *y1, double *y2, double *yF, double *yM,
        double *r1, double *r2, double *r3, double *r4,
        double *r5, double *r6, double *r7, double *r8,double Ti){
    double dVx,cof[9];
    
    cofactor3D(cof, y1, y2, yF, yM,1.0);
    
    // now there are 6 cases directly on this tetrahedra
    // 1) x1[0] = 1 --> xF[0] = 0.25 and xM[0] = 0.125 ==> write to r1[0]
    dVx = cof[0] * (1-0.125) + cof[1] * (-0.125) + cof[2] * (0.25-0.125);
    r1[0] += Ti * dVx;
    // 2) x1[1] = 1 --> xF[1] = 0.25 and xM[1] = 0.125 ==> write to r1[1]
    dVx = cof[3] * (1.0-0.125) + cof[4] * (-0.125) + cof[5] * (0.25-0.125);
    r1[1] += Ti * dVx;
    // 3) x1[2] = 1 --> xF[2] = 0.25 and xM[2] = 0.125 ==> write to r1[2]
    dVx = cof[6] * (1.0-0.125) + cof[7] * (-0.125) + cof[8] * (0.25-0.125);
    r1[2] += Ti * dVx;
    
    // 4) x2[0] = 1 --> xF[0] = 0.25 and xM[0] = 0.125 ==> write to r2[0]
    dVx = cof[0] * (-0.125) + cof[1] * (1.0-0.125) + cof[2] * (0.25-0.125);
    r2[0] +=  Ti * dVx;
    // 5) x2[1] = 1 --> xF[1] = 0.25 and xM[1] = 0.125 ==> write to r2[1]
    dVx = cof[3] * (-0.125) + cof[4] * (1.0-0.125) + cof[5] * (0.25-0.125);
    r2[1] +=  Ti * dVx;
    // 6) x2[2] = 1 --> xF[2] = 0.25 and xM[2] = 0.125 ==> write to r2[2]
    dVx = cof[6] * (-0.125) + cof[7] * (1.0-0.125) + cof[8] * (0.25-0.125);
    r2[2] +=  Ti * dVx;
    
    // there are also 6 cases not on this terahedra, but on the same face
    // 7) x3[0] = 1 --> xF[0] = 0.25 and xM[0] = 0.125 ==> write to r3[0]
    dVx = cof[0] * (-0.125) + cof[1] * (-0.125) + cof[2] * (0.25-0.125);
    r3[0] += Ti * dVx;
    r4[0] += Ti * dVx;
    // 8) x3[1] = 1 --> xF[1] = 0.25 and xM[1] = 0.125 ==> write to r3[1]
    dVx = cof[3] * (-0.125) + cof[4] * (-0.125) + cof[5] * (0.25-0.125);
    r3[1] += Ti * dVx;
    r4[1] += Ti * dVx;
    // 9) x3[2] = 1 --> xF[2] = 0.25 and xM[2] = 0.125 ==> write to r3[2]
    dVx = cof[6] * (-0.125) + cof[7] * (-0.125) + cof[8] * (0.25-0.125);
    r3[2] += Ti * dVx;
    r4[2] += Ti * dVx;
    
    
    // finally there are  12 cases left, all on the opposite surface, so interaction just via the cell-center
    // 13) x5[0] = 1 -->  xM[0] = 0.125 ==> write to r5[0]
    dVx = cof[0] * (-0.125) + cof[1] * (-0.125) + cof[2] * (-0.125);
    r5[0] += Ti * dVx;
    r6[0] += Ti * dVx;
    r7[0] += Ti * dVx;
    r8[0] += Ti * dVx;
    // 14) x5[1] = 1 -->  xM[1] = 0.125 ==> write to r5[1]
    dVx = cof[3] * (-0.125) + cof[4] * (-0.125) + cof[5] * (-0.125);
    r5[1] += Ti * dVx;
    r6[1] += Ti * dVx;
    r7[1] += Ti * dVx;
    r8[1] += Ti * dVx;
    
    // 15) x5[2] = 1 -->  xM[2] = 0.125 ==> write to r5[2]
    dVx = cof[6] * (-0.125) + cof[7] * (-0.125) + cof[8] * (-0.125);
    r5[2] += Ti * dVx;
    r6[2] += Ti * dVx;
    r7[2] += Ti * dVx;
    r8[2] += Ti * dVx;
    
    
}
void VAMPIREdiag3D(double *diag, const double *yc, const double *m, const double *Jac,
        const double *dT, const double *Tc,const double *h){
    const int m1 = m[0], m2 = m[1],m3 = m[2];
    double yA[3], yB[3], yC[3], yD[3], yE[3], yF[3], yG[3], yH[3], yABDC[3], yBDHF[3], yABFE[3], yCAEG[3], yDCGH[3], yEFHG[3], yM[3];
    double rA[3], rB[3], rC[3], rD[3], rE[3], rF[3], rG[3], rH[3];
    
    double dTi[3],Jaci,Ti;
    int id, i1, i2, i3, r, comp,idn,off1,off2,off3;
    const double hd = h[0] * h[1] * h[2];
    
    const int nn = (m1+1)*(m2+1)*(m3+1);
    const int n = (m1)*(m2)*(m3);
    
    const int inc_x = 1;
    const int inc_y = m1+1;
    const int inc_z = (m1+1)*(m2+1);
    
    const double *rptr;
    double *wptr;
    for (off1=0; off1<2; off1++) {
        for (off2=0; off2<2; off2++) {
            for (off3=0; off3<2; off3++) {
#pragma omp parallel for  default(shared) private(rptr,wptr,id, i1, i2, i3, r, comp,idn, dTi,Jaci,Ti, yA, yB, yC, yD, yE, yF, yG, yH, yABDC, yBDHF, yABFE, yCAEG, yDCGH, yEFHG, yM, rA, rB, rC, rD, rE, rF, rG, rH)
                for (i3=off3; i3<m3; i3+=2){
                    for (i2=off2; i2<m2; i2+=2){
                        for (i1=off1; i1<m1; i1+=2){
                            idn = i1+ i2* inc_y+ i3*inc_z;
                            id =  i1+m1 *(i2+i3*m2);
                            // get edges, face stg and cell center of yc (shared --> private)
                            comp = 0;
                            // get image and derivatives
                            Ti     = Tc[id]/hd;
                            dTi[0] = dT[id];
                            dTi[1] = dT[id+n];
                            dTi[2] = dT[id+2*n];
                            Jaci   = Jac[id];
                            for(r=0;r<3;r++){
                                rptr = yc+idn+comp;
                                // A
                                yA[r] = *rptr;
                                // B
                                rptr += inc_x;
                                yB[r] = *rptr;
                                // D
                                rptr += inc_y;
                                yD[r] = *rptr;
                                // C
                                rptr -= inc_x;
                                yC[r] = *rptr;
                                // G
                                rptr += inc_z;
                                yG[r] = *rptr;
                                // H
                                rptr += inc_x;
                                yH[r] = *rptr;
                                // F
                                rptr -= inc_y;
                                yF[r] = *rptr;
                                // E
                                rptr -= inc_x;
                                yE[r] = *rptr;
                                
                                yABDC[r] = 0.25*(yA[r]+yB[r]+yD[r]+yC[r]);
                                yBDHF[r] = 0.25*(yB[r]+yD[r]+yH[r]+yF[r]);
                                yABFE[r] = 0.25*(yA[r]+yB[r]+yF[r]+yE[r]);
                                yCAEG[r] = 0.25*(yC[r]+yA[r]+yE[r]+yG[r]);
                                yDCGH[r] = 0.25*(yD[r]+yC[r]+yG[r]+yH[r]);
                                yEFHG[r] = 0.25*(yE[r]+yF[r]+yH[r]+yG[r]);
                                yM[r]    = 0.5 *(yABDC[r]+yEFHG[r]);
                                
                                comp +=nn;
                                
                                rA[r] = Jaci * 0.125 * dTi[r];
                                rB[r] = Jaci * 0.125 * dTi[r];
                                rC[r] = Jaci * 0.125 * dTi[r];
                                rD[r] = Jaci * 0.125 * dTi[r];
                                rE[r] = Jaci * 0.125 * dTi[r];
                                rF[r] = Jaci * 0.125 * dTi[r];
                                rG[r] = Jaci * 0.125 * dTi[r];
                                rH[r] = Jaci * 0.125 * dTi[r];
                            }
                            
                            VAMPIREdiagTetra(yA, yB, yABDC, yM, rA, rB, rD, rC, rE, rF, rH, rG,-Ti);
                            VAMPIREdiagTetra(yB, yD, yABDC, yM, rB, rD, rC, rA, rF, rH, rG, rE,-Ti);
                            VAMPIREdiagTetra(yD, yC, yABDC, yM, rD, rC, rA, rB, rH, rG, rE, rF,-Ti);
                            VAMPIREdiagTetra(yC, yA, yABDC, yM, rC, rA, rB, rD, rG, rE, rF, rH,-Ti);
                            
                            VAMPIREdiagTetra(yB, yD, yBDHF, yM, rB, rD, rH, rF, rA, rC, rG, rE,Ti);
                            VAMPIREdiagTetra(yD, yH, yBDHF, yM, rD, rH, rF, rB, rC, rG, rE, rA,Ti);
                            VAMPIREdiagTetra(yH, yF, yBDHF, yM, rH, rF, rB, rD, rG, rE, rA, rC,Ti);
                            VAMPIREdiagTetra(yF, yB, yBDHF, yM, rF, rB, rD, rH, rE, rA, rC, rG,Ti);
                            
                            VAMPIREdiagTetra(yA, yB, yABFE, yM, rA, rB, rF, rE, rC, rD, rH, rG,Ti);
                            VAMPIREdiagTetra(yB, yF, yABFE, yM, rB, rF, rE, rA, rD, rH, rG, rC,Ti);
                            VAMPIREdiagTetra(yF, yE, yABFE, yM, rF, rE, rA, rB, rH, rG, rC, rD,Ti);
                            VAMPIREdiagTetra(yE, yA, yABFE, yM, rE, rA, rB, rF, rG, rC, rD, rH,Ti);
                            
                            VAMPIREdiagTetra(yC, yA, yCAEG, yM, rC, rA, rE, rG, rD, rB, rF, rH,Ti);
                            VAMPIREdiagTetra(yA, yE, yCAEG, yM, rA, rE, rG, rC, rB, rF, rH, rD,Ti);
                            VAMPIREdiagTetra(yE, yG, yCAEG, yM, rE, rG, rC, rA, rF, rH, rD, rB,Ti);
                            VAMPIREdiagTetra(yG, yC, yCAEG, yM, rG, rC, rA, rE, rH, rD, rB, rF,Ti);
                            
                            VAMPIREdiagTetra(yD, yC, yDCGH, yM, rD, rC, rG, rH, rB, rA, rE, rF,Ti);
                            VAMPIREdiagTetra(yC, yG, yDCGH, yM, rC, rG, rH, rD, rA, rE, rF, rB,Ti);
                            VAMPIREdiagTetra(yG, yH, yDCGH, yM, rG, rH, rD, rC, rE, rF, rB, rA,Ti);
                            VAMPIREdiagTetra(yH, yD, yDCGH, yM, rH, rD, rC, rG, rF, rB, rA, rE,Ti);
                            
                            VAMPIREdiagTetra(yE, yF, yEFHG, yM, rE, rF, rH, rG, rA, rB, rD, rC,Ti);
                            VAMPIREdiagTetra(yF, yH, yEFHG, yM, rF, rH, rG, rE, rB, rD, rC, rA,Ti);
                            VAMPIREdiagTetra(yH, yG, yEFHG, yM, rH, rG, rE, rF, rD, rC, rA, rB,Ti);
                            VAMPIREdiagTetra(yG, yE, yEFHG, yM, rG, rE, rF, rH, rC, rA, rB, rD,Ti);
                            
                            // write into result vector (private --> shared)
                            comp = 0;
                            for (r=0; r < 3; r++){
                                // scale vecs
                                // get pointer to A
                                wptr = diag + idn + comp;
                                *wptr +=  hd * rA[r] * rA[r];
                                // B
                                wptr += inc_x;
                                *wptr +=  hd * rB[r]* rB[r];
                                // D
                                wptr += inc_y;
                                *wptr += hd * rD[r]* rD[r];
                                // C
                                wptr -= inc_x;
                                *wptr += hd * rC[r]* rC[r];
                                // G
                                wptr += inc_z;
                                *wptr += hd * rG[r]* rG[r];
                                // H
                                wptr += inc_x;
                                *wptr += hd * rH[r]* rH[r];
                                // F
                                wptr -= inc_y;
                                *wptr += hd * rF[r]* rF[r];
                                // E
                                wptr -= inc_x;
                                *wptr += hd * rE[r]* rE[r];
                                comp +=nn;
                            }
                        }
                    }
                }
            }
        }
    }
    

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