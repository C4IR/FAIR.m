%==============================================================================
% This code is part of the VAMPIRE app for the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR/FAIR.m 
%==============================================================================
%
% VAMPIRE - Variational Algorithm for Mass-Preserving Image REgistration
% 
% 
%         xxxxxxxx                                           xxxx          
%      xxxxxxxxxxxxxx                                     xxxxxxxxxxx      
%    xxxxxxxxxxxxxxxxx                                 xxxxxxxxxxxxxxxx    
%   xxxxxxxxxxxxxxxxxxxx                             xxxxxxxxxxxxxxxxxxxx  
%  xxxxxxxxxxxxxxxxxxxxxxx                          xxxxxxxxxxxxxxxxxxxxxx 
%       xxxxxxxxxxxxxxxxxxx                       xxxxxxxxxxxxxxxxxxxxx    
%         xxxxxxxxxxxxxxxxxxx     xxx   xxx      xxxxxxxxxxxxxxxxxxx       
%         xxxxxxxxxxxxxxxxxxxx    xxxxxxxxx    xxxxxxxxxxxxxxxxxxxx        
%          xxxxxxxxxxxxxxxxxxxxx  xxxxxxxxx  xxxxxxxxxxxxxxxxxxxxx         
%          xxxx  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx         
%                   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx       xx         
%                     xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx                    
%                      xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx                     
%                      x       xxxxxxxxxxxxxxxx    xxx                     
%                               xxxxxxxxxxxxxx                             
%                               xxx  xxxxxxxx                              
%                               x          xx                              
% 
%
% For a quick example, run: EV_3Dmouse_VAMPIRE.m
%
% Contents of VAMPIRE toolbox:
%
% VAMPIRENPIRobjFctn.m   - objective function for non-parametric registration
% VAMPIREPIRobjFctn.m    - objective function for non-parametric registration
% VAMPIREsolveGN_PCG.m   - matrix-free PCG solver for Gauss-Newton system
%
% Examples:
%
% EV_2DGaussian_VAMPIRE.m        - 2D example, non-parametric mass-preserving 
% EV_2DGaussian_VAMPIRE_MLPIR.m  - 2D example, parametric mass-preserving 
% EV_2DGaussianNoise_VAMPIRE.m   - 2D example, non-parametric mass-preserving 
% EV_3Dmouse_VAMPIRE.m           - 3D example, non-parametric mass-preserving
% EV_3Dmouse_VAMPIRE_MLPIR.m     - 3D example, parametric mass-preserving
%
% ==================================================================================

function debit = contents
if nargout == 0, help(mfilename); return; end;

debit = {
          'VAMPIRENPIRobjFctn.m'
          'VAMPIREPIRobjFctn.m'
          'VAMPIREsolveGN_PCG.m'
          'contents.m'
          'README.md'
          'LICENSE'
          'examples'
        };
