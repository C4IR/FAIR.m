# VAMPIRE.m
VAMPIRE stands for Variational Algorithm for Mass-Preserving Image REgistration and is an add-on to the FAIR toolbox; see https://github.com/C4IR/FAIR.m for details and license issues and VAMPIRE/contents.m for overview and advice on VAMPIRE.

## What is it? 
VAMPIRE is a smart implementation of a mass-preserving registration scheme and thus adequate in situations where the overall integral of image intensities needs to be preserved in the registration. VAMPIRE has been used for motion correction in cardiac Positron Emission Tomography (PET) and other applications.

The mathematics and the algorithm are described in 

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

## Main Developers until 2017
- Lars Ruthotto, Dept. of Mathematics and Computer Science, Emory University, Atlanta, USA
- Fabian Gigengack, European Institute for Molecular Imaging (EIMI), University of Münster, Germany
