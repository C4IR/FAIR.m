%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
% ##2
%==============================================================================
%
% Tutorial for FAIR: regulariztion
%
% given a certain force field f on a 2D domain, this tutorial 
% computes and visualizes the elastic displacement u, such that 
%
%  B'*B u = f
%
% where B is the discrete elastic operator 
% see also getElasticMatrixStg
%==============================================================================

clear, close all, help(mfilename);

                    % 2D example
omega  = [0,1,0,1]; % physical domain
m      = [16,12];   % number of discretization points
mu     = 1;         % Lame constants, control elasticity properties
lambda = 0;         % Youngs modulus and Poisson ratio
n      = prod(m);   % number of cells, and staggered dimensions
ns     = [(m(1)+1)*m(2),m(1)*(m(2)+1)];

% generate cell centered, staggered and nodal grids
% note: computation is staggered, others used for visualization
xc = getCellCenteredGrid(omega,m);
xn = getNodalGrid(omega,m);

% build elasticity operator on a staggered grid and visualize
B = getElasticMatrixStg(omega,m,mu,lambda);
FAIRfigure(1); clf; subplot(2,2,1); spy(B); title('B elastic on staggered grid')

% create a force field o cell centered grid and visualize
fc  = [0*xc(1:n);-10*sin(pi*xc(1:n))];
% computations on staggered grid
fs = grid2grid(fc,m,'centered','staggered');

% compute the displacement from (B'*B)*us = fs
% note B is derivative operator and has null space (constants)

% remove constant parts: E'*(f-E*w) = 0
E = [ones(ns(1),1)*[1,0];ones(ns(2),1)*[0,1]];
fs = fs - E*((E'*E)\(E'*fs));

% solve for displacements
warning off % matrix is singular, but MATLAB doesn't care
us = (B'*B)\fs;
warning on

% visualization on cell centered and nodal grids
uc = grid2grid(us,m,'staggered','centered');
un = grid2grid(us,m,'staggered','nodal');

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
