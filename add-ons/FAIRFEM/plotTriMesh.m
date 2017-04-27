%==============================================================================
% This code is part of the Finite Element Method app for the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR/FAIRFEM 
%==============================================================================
%
% function ph = plotTriMesh(Mesh,xn,varargin)
%
% Input:
%  Mesh     - description
%  xn       - description
%  varargin - description
%
% Output:
%  ph       - description
%
% see also plotGrid
% =========================================================================
function ph = plotTriMesh(Mesh,yc,varargin)

if nargin==0
    runMinimalExample
    return
end

zdata  = zeros(Mesh.nnodes,1);
cdata   =zeros(Mesh.ntri ,1);
angle  = [0 90];
for k=1:2:length(varargin), % overwrites default parameter
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

yc = reshape(yc,[],2);
ph = trisurf(Mesh.tri,yc(:,1),yc(:,2),zdata,cdata); 
view(angle)    ; axis(Mesh.omega);axis xy;
axis image

function runMinimalExample

omega = [0 1 0 2];
m     = [4 7];
Mesh = TriMesh3(omega,m);
figure(1); clf;
subplot(1,2,1);
feval(mfilename,Mesh,Mesh.xn);
title('identity');
subplot(1,2,2);
feval(mfilename,Mesh,Mesh.xn+randn(size(Mesh.xn))*2e-2);
title('deformed');