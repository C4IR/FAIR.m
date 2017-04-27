# FAIRFEM - Image Registration using Finite Elements

this is an extension for the Matlab based toolbox [FAIR.m](https://github.com/C4IR/FAIR.m) which must be installed before using this code. 

## What's new?

This package provides fininte element discretizations of elastic and hyperelastic image registration methods on triangular and tetrahedral meshes. Using these meshes is advantageous for hyperelastic registration to ensure that the discrete transformation is free of foldings (i.e., is a diffeomorphism). It also simplifies investigating the impact of triangulations on the computed transformation. 

## References

The numerical method is described in detail in:

- Ruthotto, L., & Modersitzki, J. (2015). 
   [Non-linear Image Registration](http://doi.org/10.1007/978-1-4939-0790-8_39). 
   In Handbook of Mathematical Methods in Imaging (pp. 2005–2051). New York, NY: Springer New York. 

The mass-preserving registration model is developed in:

- Gigengack, F., Ruthotto, L., Burger, M., Wolters, C. H., Jiang, X., & Schafers, K. P. (2012). 
  [Motion Correction in Dual Gated Cardiac PET Using Mass-Preserving Image Registration](http://doi.org/10.1109/TMI.2011.2175402). Medical Imaging, IEEE Transactions on, 31(3),   
  698–712. 


