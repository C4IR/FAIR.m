# A Lagrangian Gauss–Newton–Krylov solver for mass- and intensity-preserving diffeomorphic image registration
see https://github.com/C4IR/LagLDDMM.m for details and license issues.
    
## What is it? 

This code adds new options for diffeomorphic registration to the [FAIR](https://github.com/C4IR/LagLDDMM.m) toolbox. The mathematical framework used here is based on the Large Deformation Diffeomorphic Metric Mapping (LDDMM) approach. What's making our approach special is the numerical realization (using a Lagrangian PDE solver to eliminate the constraint) and the flexibility (we support both intensity- and mass-preservation).

A solid documentation and reference is

    Mang, Ruthotto: A Lagrangian Gauss–Newton–Krylov Solver For Mass- And Intensity-Preserving Diffeomorphic Image Registration

    @article{MangRuthotto2017,
      Title = {A {L}agrangian {G}auss--{N}ewton--{K}rylov solver for mass- and intensity-preserving diffeomorphic image registration},
      Year = {2017},
      Journal = {SIAM Journal on Scientific Computing},
      Author = {A. Mang, L. Ruthotto},
    }
see also https://arxiv.org/abs/1703.04446 for a free preprint.

## Getting started

Before running this code, make sure you have downloaded FAIR and added it to your MATLAB path. Having done this, add this folder and subfolders to the path. The examples from the paper can be found in the `examples` folder.

## Acknowledgements

This material is in part based upon work supported by the National Science Foundation under Grant Number 1522599. Any opinions, findings, and conclusions or recommendations expressed in this material are those of the author(s) and do not necessarily reflect the views of the National Science Foundation.
