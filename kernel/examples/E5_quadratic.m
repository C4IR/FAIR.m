%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% Tutorial for FAIR:  Quadratic LM registration for hand example
%
% - load data (see setup2DhandData)
% - compute least squares solution and plot
% - this file for demonstration, see also LMreg for a convenient variant
%==============================================================================

clear, close all, help(mfilename)

setup2DhandData; 

xc = reshape(getCellCenteredGrid(omega,m),[],2); %LM(7,:) = [];
Q  = [ones(size(LM,1),1),LM(:,[3:4]),LM(:,3).^2,LM(:,4).^2,LM(:,3).*LM(:,4)];
wc = (Q'*Q)\(Q'*LM(:,1:2));

quad = @(w,x1,x2) (w(1)+w(2)*x1+w(3)*x2+w(4)*x1.^2+w(5)*x2.^2+w(6)*x1.*x2); 
yc = [quad(wc(:,1),xc(:,1),xc(:,2));quad(wc(:,2),xc(:,1),xc(:,2))];
LM(:,[5,6]) = [quad(wc(:,1),LM(:,3),LM(:,4)),quad(wc(:,2),LM(:,3),LM(:,4))];

P5_LM; % for nice plots
%==============================================================================
