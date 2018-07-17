%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
% ##2
%==============================================================================
%
% Tutorial for FAIR: regulariztion
%
% given a certain force field f on a 2D domain, this tutorial 
% computes and visualizes the curvature regularized displacement u, such that 
%
%  B'*B u = f
%
% where B is the discrete curvature operator 
% see also getCurvatureMatrix
%==============================================================================

clear, close all, help(mfilename);

                    % 2D example
omega  = [0,1,0,1]; % physical domain
m      = [16,12];   % number of discretization points
n      = prod(m);   % number of cells, and staggered dimensions
hd     = prod((omega(2:2:end)-omega(1:2:end))./m);
alpha  = 10;

% generate cell centered and nodal grids
% note: computation is cell centered, others used for visualization
xc = getCellCenteredGrid(omega,m);
xn = getNodalGrid(omega,m);

% build elasticity operator on a staggered grid and visualize
B = getCurvatureMatrix(omega,m);
FAIRfigure(1); clf; subplot(2,2,1); spy(B); title('B curcature on cell-centered grid')

% create a force field o cell centered grid and visualize
fc  = [0*xc(1:n);-10*sin(pi*xc(1:n))];

% compute the displacement from (B'*B)*us = fs
% note B is derivative operator and has null space (constants)

% remove constant parts: E'*(f-E*w) = 0
E = [ones(n,1)*[1,0];ones(n,1)*[0,1]];
fc = fc - E*((E'*E)\(E'*fc));

% solve for displacements
warning off % matrix is singular, but MATLAB doesn't care
uc = hd*alpha*(B'*B)\fc;
warning on

% visualization on cell centered and nodal grids
un = grid2grid(uc,m,'centered','nodal');

% shortcuts for visualization, 
vecNorm   = @(v) sqrt(sum(reshape(v,[],2).^2,2));
viewField = @(v) viewImage2Dsc(vecNorm(v),omega,m); 
PG        = @(x) plotGrid(x,omega,m);

% visualize force field
FAIRfigure(1); subplot(2,2,3);
viewField(fc); hold on; colormap(gray); colorbar; PG(xn); 
qh = quiver(xc(1:n),xc(n+1:end),fc(1:n),fc(n+1:end),1);
set(qh,'color','r','linewidth',1)
title('a force field')

% visualize displacement field
FAIRfigure(1); subplot(2,2,4);
viewField(uc); hold on; colormap(gray); colorbar; PG(xn);
qh = quiver(xc(1:n),xc(n+1:end),uc(1:n),uc(n+1:end),1);
set(qh,'color','g','linewidth',1)
title('the displacement field')

% visuaize the displaced grid
FAIRfigure(1); subplot(2,2,2);
xh = PG(xn); hold on; axis equal image; yh = PG(xn+un); 
set(xh,'linewidth',2); set(xh,'linewidth',2,'color','g');
title('initial and displaced grids')

%==============================================================================
