# FAIRFEM - Image Registration using Finite Elements

this is an add-on for the Matlab based toolbox FAIR; see https://github.com/C4IR/FAIR.m for details and license issues and FAIRFEM/contents.m for an overview.  

## What's new?

This package provides fininte element discretizations of elastic and hyperelastic image registration methods on triangular and tetrahedral meshes. Using these meshes is advantageous for hyperelastic registration to ensure that the discrete transformation is free of foldings (i.e., is a diffeomorphism). It also simplifies investigating the impact of triangulations on the computed transformation. 

## References

The numerical method is described in detail in:

       @inbook{RuthottoModersitzki2015,
                author = {Ruthotto, L and Modersitzki J},
                title = {{Non-linear Image Registration}},
                book = {Handbook of Mathematical Methods in Imaging},
                editor = {Scherzer O},
                year = {2015},
                publisher = {Springer New York},
                pages = {2005--2051},
                keywords = {Image Registration},
                doi = {10.1007/978-1-4939-0790-8_39}. 
                }


The mass-preserving registration model is developed in:

       @article{GigengackEtAl2012,
                author = {Gigengack, F and Ruthotto, L and Burger, M and Wolters, C H and Jiang, Xiaoyi and Schafers, K P},
                title = {{Motion Correction in Dual Gated Cardiac PET Using Mass-Preserving Image Registration}},
                journal = {Medical Imaging, IEEE Transactions on},
                year = {2012},
                volume = {31},
                number = {3},
                pages = {698--712},
                keywords = {Image Registration},
                doi = {10.1109/TMI.2011.2175402},
                }


