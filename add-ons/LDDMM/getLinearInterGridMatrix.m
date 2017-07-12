%==========================================================================
% ##1
%
% function Ty = getLinearInterGridMatrix(omega,m,yc)
%
% computes sparse matrix Ty such that
%
% Ty*v = linearInterGrid(vc,omega,m,y)
%
% Input:
%  omega   - domain of grid to be interpolated
%  m       - number of cells in grid to be interpolated
%  yc      - interpolation points
%
% Output:
%  Ty     - interpolation matrix (sparse_
%
% see also linearInterGrid, getLinearInterMatrix
% =========================================================================
function Ty = getLinearInterGridMatrix(omega,m,yc)

dim = numel(omega)/2;
switch regularizer('get','grid')
    case 'nodal'
        h    = (omega(2:2:end)-omega(1:2:end))./m;
        omegaNodal = omega + reshape([-h;h]/2,1,[]);
        Ty = getLinearInterMatrix(omegaNodal,m+1,yc);
        Ty = kron(speye(dim),Ty);
    case 'staggered'
        h    = (omega(2:2:end)-omega(1:2:end))./m;
        omegaStg1 = omega; omegaStg1(1:2) = omegaStg1(1:2) + .5*[-h(1) h(1)];
        omegaStg2 = omega; omegaStg2(3:4) = omegaStg2(3:4) + .5*[-h(2) h(2)];
        
        if dim==2
            Ty1 = getLinearInterMatrix(omegaStg1,m+[1,0],yc);
            Ty2 = getLinearInterMatrix(omegaStg2,m+[0,1],yc);
            Ty = blkdiag(Ty1,Ty2);
        else
            omegaStg3 = omega; omegaStg3(5:6) = omegaStg3(5:6) + .5*[-h(3) h(3)];
            Ty1 = getLinearInterMatrix(omegaStg1,m+[1,0,0],yc);
            Ty2 = getLinearInterMatrix(omegaStg2,m+[0,1,0],yc);
            Ty3 = getLinearInterMatrix(omegaStg3,m+[0,0,1],yc);
            Ty = blkdiag(Ty1,Ty2,Ty3);
        end
    otherwise % assume velocities are cell-centered
        Ty  = getLinearInterMatrix(omega,m,yc);
        Ty = kron(speye(dim),Ty);      
end