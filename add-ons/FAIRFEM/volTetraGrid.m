%==============================================================================
% This code is part of the Finite Element Method app for the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR/FAIRFEM 
%==============================================================================
%
% function V = volTetraGrid(Mesh,yc,varargin)
%
% computes volumes of tetrahedral grid
%
% Input:
%   Mesh - description of mesh (struct)
%   yc   - nodes of finite element grid
%
% Output:
%   V    - volumes of tetrahedra
%
%==============================================================================
function [V,dV] = volTetraGrid(Mesh,yc,varargin)
dV = [];
if nargin==0, help(mfilename); runMinimalExample; return; end;

matrixFree   = 0;
doDerivative = (nargout>1);
for k=1:2:length(varargin),     % overwrites default parameter
    eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

dim = Mesh.dim;

yc = reshape(yc,[],dim);

switch dim
    case 2
        % project to points
        y1 = Mesh.mfPi(yc,1);
        y2 = Mesh.mfPi(yc,2);
        y3 = Mesh.mfPi(yc,3);
        
        % compute volume
        V = ((y1(:,1)-y3(:,1)).*(y2(:,2)-y3(:,2)) ...
            - (y2(:,1)-y3(:,1)).*(y1(:,2)-y3(:,2)))/2;
        if doDerivative,
            w = [y2(:,2)-y3(:,2), y1(:,2)-y3(:,2), ...
                y1(:,1)-y3(:,1), y2(:,1)-y3(:,1)];
            if not(matrixFree),
                P1 = Mesh.P1;   P2 = Mesh.P2;    P3 = Mesh.P3;
                dV =[
                    sdiag(w(:,1))*(P1-P3) ...
                    - sdiag(w(:,2))*(P2-P3),...
                    sdiag(w(:,3))*(P2-P3) ...
                    - sdiag(w(:,4))*(P1-P3)
                    ]/2;
            else
                P1 = @(v,i) Mesh.mfPi(v(:,i),1);
                P2 = @(v,i) Mesh.mfPi(v(:,i),2);
                P3 = @(v,i) Mesh.mfPi(v(:,i),3);
                rshp = @(v) reshape(v,[],2);
                
                dV.dV    = @(v)(   w(:,1).*(P1(rshp(v),1)-P3(rshp(v),1)) ...
                    - w(:,2).*(P2(rshp(v),1)-P3(rshp(v),1)) ...
                    + w(:,3).*(P2(rshp(v),2)-P3(rshp(v),2)) ...
                    - w(:,4).*(P1(rshp(v),2)-P3(rshp(v),2)))/2;
                
                dV.dVadj = @(v)[  (P1(w(:,1).*v(:),1)-P3(w(:,1).*v(:),1)) ...
                    - (P2(w(:,2).*v(:),1)-P3(w(:,2).*v(:),1)); ...
                    + (P2(w(:,3).*v(:),1)-P3(w(:,3).*v(:),1)) ...
                    - (P1(w(:,4).*v(:),1)-P3(w(:,4).*v(:),1))]/2;
            end
        end
    case 3
        P1 = @(x) Mesh.mfPi(x,1);   P2 = @(x) Mesh.mfPi(x,2);    
        P3 = @(x) Mesh.mfPi(x,3);   P4 = @(x) Mesh.mfPi(x,4);
        % project to edges
        e1 = P1(yc)-P4(yc);
        e2 = P2(yc)-P4(yc);
        e3 = P3(yc)-P4(yc);
        
        V  = (1/6)*(e1(:,1) .* e2(:,2) .* e3(:,3) ...
            + e2(:,1) .* e3(:,2) .* e1(:,3) ...
            + e3(:,1) .* e1(:,2) .* e2(:,3) ...
            - e1(:,3) .* e2(:,2) .* e3(:,1) ...
            - e2(:,3) .* e3(:,2) .* e1(:,1) ...
            - e3(:,3) .* e1(:,2) .* e2(:,1));
        if doDerivative,
            % cofactor matrix
            cof = [
                    e2(:,2) .* e3(:,3) - e2(:,3) .* e3(:,2),...
                    e3(:,2) .* e1(:,3) - e3(:,3) .* e1(:,2),...
                    e1(:,2) .* e2(:,3) - e1(:,3) .* e2(:,2),...
                    e3(:,1) .* e2(:,3) - e3(:,3) .* e2(:,1),...
                    e1(:,1) .* e3(:,3) - e1(:,3) .* e3(:,1),...
                    e2(:,1) .* e1(:,3) - e2(:,3) .* e1(:,1),...
                    e2(:,1) .* e3(:,2) - e2(:,2) .* e3(:,1),...
                    e3(:,1) .* e1(:,2) - e3(:,2) .* e1(:,1),...
                    e1(:,1) .* e2(:,2) - e1(:,2) .* e2(:,1),...
                  ];
          

              if not(matrixFree),
                P1 = Mesh.P1;   P2 = Mesh.P2;    P3 = Mesh.P3; P4 = Mesh.P4;
                dV = (1/6)*[ ...
                              sdiag(cof(:,1))* (P1-P4) ...
                            + sdiag(cof(:,2))* (P2-P4) ...
                            + sdiag(cof(:,3))* (P3-P4), ...
                            + sdiag(cof(:,4))* (P1-P4) ...
                            + sdiag(cof(:,5))* (P2-P4) ...
                            + sdiag(cof(:,6))* (P3-P4), ...
                            + sdiag(cof(:,7))* (P1-P4) ...
                            + sdiag(cof(:,8))* (P2-P4) ...
                            + sdiag(cof(:,9))* (P3-P4)];
              else
                P1 = @(v,i) P1(v(:,i));
                P2 = @(v,i) P2(v(:,i));
                P3 = @(v,i) P3(v(:,i));
                P4 = @(v,i) P4(v(:,i));
                rshp = @(v) reshape(v,[],3);
                
                dV.dV = @(v)(1/6)*( ...
                                     cof(:,1).* (P1(rshp(v),1)-P4(rshp(v),1)) ...
                                   + cof(:,2).* (P2(rshp(v),1)-P4(rshp(v),1)) ...
                                   + cof(:,3).* (P3(rshp(v),1)-P4(rshp(v),1)) ...
                                   + cof(:,4).* (P1(rshp(v),2)-P4(rshp(v),2)) ...
                                   + cof(:,5).* (P2(rshp(v),2)-P4(rshp(v),2)) ...
                                   + cof(:,6).* (P3(rshp(v),2)-P4(rshp(v),2)) ...
                                   + cof(:,7).* (P1(rshp(v),3)-P4(rshp(v),3)) ...
                                   + cof(:,8).* (P2(rshp(v),3)-P4(rshp(v),3)) ...
                                   + cof(:,9).* (P3(rshp(v),3)-P4(rshp(v),3)));
                dV.dVadj = @(v)(1/6)*[ ...
                                         (P1(cof(:,1).*v(:),1)-P4(cof(:,1).*v(:),1)) ...
                                      +  (P2(cof(:,2).*v(:),1)-P4(cof(:,2).*v(:),1)) ...
                                      +  (P3(cof(:,3).*v(:),1)-P4(cof(:,3).*v(:),1)); ...
                                      +  (P1(cof(:,4).*v(:),1)-P4(cof(:,4).*v(:),1)) ...
                                      +  (P2(cof(:,5).*v(:),1)-P4(cof(:,5).*v(:),1)) ...
                                      +  (P3(cof(:,6).*v(:),1)-P4(cof(:,6).*v(:),1)); ...
                                      +  (P1(cof(:,7).*v(:),1)-P4(cof(:,7).*v(:),1)) ...
                                      +  (P2(cof(:,8).*v(:),1)-P4(cof(:,8).*v(:),1)) ...
                                      +  (P3(cof(:,9).*v(:),1)-P4(cof(:,9).*v(:),1))];
              end
            
            
            
            
        end
    otherwise
        error('only dim=2,3.');
end
function D = sdiag(v)
D =  spdiags(v(:),0,numel(v),numel(v));

function runMinimalExample
omega = [0 3 0 5]; m = [2 4]; type = 1;
Mesh  = TriMesh2(omega,m);
xc    = Mesh.xn + 1e-1*randn(size(Mesh.xn));

[V,dV]     = volTetraGrid(Mesh,xc,'matrixFree',0);
[Vmf,dVmf] =  volTetraGrid(Mesh,xc,'matrixFree',1);

fctn = @(x)volTetraGrid(Mesh,x);
checkDerivative(fctn,xc(:));

fig = figure(1); clf;
subplot(1,2,1);
trisurf(Mesh.tri,xc(:,1),xc(:,2),0*xc(:,1),V); view([0 90]);
axis(omega);
title(sprintf('sum(V)=%f',sum(V)));
err1 = norm(V(:)-Vmf(:))/norm(V(:));
fprintf('%s, mesh type %d,    MFerror V: %1.3e OK? %d\n',mfilename,type,err1,err1<1e-14);

v1 = randn(size(Mesh.tri,1),1);
v2 = randn(numel(Mesh.xn),1);
err2 = norm(dV*v2-dVmf.dV(v2))/norm(dV*v2);
err3 = norm(dVmf.dV(v2)'*v1 - dVmf.dVadj(v1)'*v2);
fprintf('%s, mesh type %d,    MFerror dV: %1.3e OK? %d\n',mfilename,type,err2,err2<1e-14);
fprintf('%s, mesh type %d, MFadjerror dV: %1.3e OK? %d\n',mfilename,type,err3,err3<1e-14);

omega = [0 3 0 5 0 2]; m = [2 4 3]; type = 1;
Mesh  = TetraMesh1(omega,m)
xc    = Mesh.xn + 1e-1*randn(size(Mesh.xn));

[V,dV]   = volTetraGrid(Mesh,xc,'matrixFree',0);
[Vmf,dVmf] =  volTetraGrid(Mesh,xc,'matrixFree',1);

fctn = @(x)volTetraGrid(Mesh,x);
checkDerivative(fctn,xc(:));

figure(fig)
subplot(1,2,2);
tetramesh(Mesh.tri,xc,V,'FaceAlpha',.3);

axis(omega);
title(sprintf('sum(V)=%f',sum(V)));
err1 = norm(V(:)-Vmf(:))/norm(V(:));
fprintf('%s, mesh type %d, MFerror: %1.3e OK? %d\n',mfilename,type,err1,err1<1e-14);

v1 = randn(size(Mesh.tri,1),1);
v2 = randn(numel(Mesh.xn),1);
err2 = norm(dV*v2-dVmf.dV(v2))/norm(dV*v2);
err3 = norm(dVmf.dV(v2)'*v1 - dVmf.dVadj(v1)'*v2);
fprintf('%s, mesh type %d,    MFerror dV: %1.3e OK? %d\n',mfilename,type,err2,err2<1e-14);
fprintf('%s, mesh type %d, MFadjerror dV: %1.3e OK? %d\n',mfilename,type,err3,err3<1e-14);





