# FAIR: Flexible Algorithms for Image Registration 
see https://github.com/C4IR/FAIR.m for details and license issues.

## What is it? 
FAIR stands for Flexible Algorithms for Image Registration and is a package written in MATLAB which can be used for solving the so-called image registration problem (alas co-registration, correspondence, fusion, matching, optical flow, , warping, ...). 
A documentation is the book 

    Jan Modersitzki: FAIR - Flexible Algorithms for Image Registration, SIAM 2009

    @BOOK{2009-FAIR,
      title = {{FAIR}: Flexible Algorithms for Image Registration},
      publisher = {SIAM},
      address = {Philadelphia},
      year = {2009},
      author = {J. Modersitzki},
    }

The objective is to automatically establish point-to-point correspondences between objects in two different scenes. Given are two (or more) images T and R. Wanted is a geometrical transformation y (correction, deformation, displacement, distortion, ...) such that the transformed image T(y) is similar to R. Given a similarity/distance measure D (see distances for options) and a regularization for y (see regularizer for options) the problem is phrased in a variational setting

    Minimize the joint energy J with respect to y
    J(y) = D(T(y),R) + S(y)
    
The toolbox provides a variety of 
    
    - example data             kernel/data
    - image models             kernel/imgModels
    - distance measures        kernel/diatnces
    - regularizer              kernel/regularizers
    - transfomation models     kernel/transformations
    - image viewers            kernel/viewers
    
Additionally, the toolbox provides

    - numerical tools          kernel/numeric
    - matrix-free operators    kernel/matrixFree
    - landmark registrations   kernel/landmarks
    - genral tools (i/o etc)   kernel/tools
    

##  Purpose? 
FAIR is primarily designed as an academic and teaching tool. It can be used to explore existing techniques or to invent new features (please let us know and contribute to the community). Though the focus is on exploring methods, FAIR can be - and already has been - used as a prototyping tool for practically relevant registration problems.

## What is new? 
The main differences to older versions of the codes are
- now available  via github and direct communication
- replaced the interpolation toolbox by an image model
- includes a nonlinear hyperelastic regularization energy
- the package is supplemented by parallelized C versions of otherwise computationally costly mfiles

## What is next?
- Provide even better tools (enhance the pieces)
- Provide add-ons for specific application




### Developers
- Jan Modersitzki Institute of Mathematics and Image Computing, University of Lübeck, Germany
- Lars Ruthotto, Dept. of Mathematics and Computer Science, Emory University, Atlanta, Trump Country
- Fabian Gigengack, European Institute for Molecular Imaging (EIMI), University of Münster, Germany

### Demos and Screenshots 
The following list presents all intermediate results and the m-file for a non-parametric multi-level image registration with an affine linear pre-registration for 2D PET-CT data using the normalized gradient field as a distance measure and elastic regularization; a Gauss-Newton scheme is used as optimizer in a matrix free setting (using multigrid for the linear solve). The time was about 10sec for 128x128 data on a standard PC including all graphics and data pre-processing.

[mfile] [data] [multi-level representation] [PIR results level 4] [NPIR results level 4] [NPIR results level 5] [NPIR results level 6] [NPIR results level 7] [final results] [output]


## Available Apps
FAIR not only provides opportunities to develop specific building blocks but also to redesign fundamental structures. As modifications can be severe, we decided to split the toolbox into a kernel and add-ons. The stand alone kernel provides all the basic functionality. The add-ons can be tuned to even overlay standard behaviour. In the following, we present some examples. We hope to extend this list with the help of your contributions.

Following apps will be available soon:
- EPI please contact Fabian Gigengack or Lars Ruthotto
- DTI Correction please contact Fabian Gigengack or Lars Ruthotto
- VAMPIRE please contact Fabian Gigengack or Lars Ruthotto


