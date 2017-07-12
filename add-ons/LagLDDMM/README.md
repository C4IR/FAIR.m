# A Lagrangian Gauss–Newton–Krylov solver for mass- and intensity-preserving diffeomorphic image registration
see https://github.com/C4IR/FAIR.m/master/tree/add-ons/LagLDDMM for details and license issues.

| <img src="https://github.com/C4IR/FAIR.m/blob/master/pictures/LDDMM.jpg"  width="300"  /> |
    
## What is it? 

This code adds new options for diffeomorphic registration to the [FAIR](https://github.com/C4IR/FAIR.m) toolbox. The mathematical framework used here is based on the [Large Deformation Diffeomorphic Metric Mapping (LDDMM)](https://en.wikipedia.org/wiki/Large_deformation_diffeomorphic_metric_mapping#Hamiltonian_LDDMM_for_Dense_Image_Matching) approach. What sets  our approach apart is the numerical realization (using a Lagrangian PDE solver to eliminate the constraint) and the flexibility (we support both intensity- and mass-preserving registration).

A solid documentation and reference is

    Mang, Ruthotto: A Lagrangian Gauss–Newton–Krylov Solver For Mass- And Intensity-Preserving Diffeomorphic Image Registration

    @article{MangRuthotto2017,
      Title = {A {L}agrangian {G}auss--{N}ewton--{K}rylov solver for mass- and intensity-preserving diffeomorphic image registration},
      Year = {2017},
      Journal = {SIAM Journal on Scientific Computing},
      Author = {A. Mang, L. Ruthotto},
    }

see also https://arxiv.org/abs/1703.04446 for a free preprint.

## Acknowledgements

This material is based upon work supported by the National Science Foundation under Grant Number 1522599 and U.S. Department of Energy, Office of Science, Office of Advanced Scientific Computing Research, Applied Mathematics program under Award Numbers DE-SC0010518 and DE-SC0009286; and by NIH grant 10042242. Any opinions, findings, and conclusions or recommendations expressed in this material are those of the author(s) and do not necessarily reflect the views of funding agencies.
