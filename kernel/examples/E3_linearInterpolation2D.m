%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
% 
% Tutorial for FAIR: linear interpolation in 2D
%
%==============================================================================

dataT = flipud([1,2,3,4;1,2,3,4;4,4,4,4])'; 
m     = size(dataT); 
omega = [0,m(1),0,m(2)]; 
Tcoef = dataT;  % for linear interpolation, coefficients=data


% define a grid and an image model
grid = @(m) getCellCenteredGrid(omega,m);
img  = @(x) linearInter(Tcoef,omega,x(:));

M     = {m,10*m}; % two resolutions, coarse and fine
xc = grid(M{1});     % coarse resolution
xf = grid(M{2});     % fine   resolution
Tc = img(xc);
Tf = img(xf);

xc = reshape(xc,[],2);
xf = reshape(xf,[M{2},2]);
FAIRfigure(2); clf; 
ph = plot3(xc(:,1),xc(:,2),Tc(:),'ro'); hold on;
qh = surf(xf(:,:,1),xf(:,:,2),reshape(Tf,M{2}));
%==============================================================================
