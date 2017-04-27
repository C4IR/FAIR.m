# VAMPIRE.m
see https://github.com/C4IR/VAMPIRE.m for details and license issues.

## What is it? 
VAMPIRE stands for Variational Algorithm for Mass-Preserving Image REgistration and is an extension to the FAIR toolbox written in MATLAB (https://de.mathworks.com/). The VAMPIRE model is mass-preserving by design and thus adequate in situations where the overall integral of image intensities needs to be preserved in the registration. For example, VAMPIRE has been used for motion correction in cardiac Positron Emission Tomography (PET).

The approach is described in 

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

Please cite this work when using VAMPIRE in your research.

## Getting started

VAMPIRE is an extension to the FAIR toolbox. Please obtain a copy at https://github.com/C4IR/FAIR.m and make sure both packages are in your MATLAB path. Finally, take a look at some examples:

```
EV_3Dmouse_VAMPIRE
```

## Main Developers until 2017
- Lars Ruthotto, Dept. of Mathematics and Computer Science, Emory University, Atlanta, USA
- Fabian Gigengack, European Institute for Molecular Imaging (EIMI), University of Münster, Germany
