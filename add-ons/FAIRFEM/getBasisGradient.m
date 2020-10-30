% For details see 
% - https://github.com/C4IR/FAIRFEM 
%==============================================================================
%
% function [dphi] = getBasisGradient(Mesh)
%
% builds gradient operator for piecewise linear basis functions
%
% Input:
%
% Mesh          - description of mesh, struct
%
% Output:
%
% [dphi] - discrete partial derivative operators
%
%==============================================================================

function [dphi] = getBasisGradient(Mesh)

xn = Mesh.xn;
vol = Mesh.vol;
if Mesh.dim == 3
% compute edges
    v1   = Mesh.mfPi(xn,1);
    v2   = Mesh.mfPi(xn,2);
    v3   = Mesh.mfPi(xn,3);
    v4   = Mesh.mfPi(xn,4);
    e1   =  v1 - v4;
    e2   =  v2 - v4;
    e3   =  v3 - v4;
    % compute inverse transformation to reference element
    cofA =  [
          e2(:,2).*e3(:,3)-e2(:,3).*e3(:,2), ...
        -(e1(:,2).*e3(:,3)-e1(:,3).*e3(:,2)), ...
          e1(:,2).*e2(:,3)-e1(:,3).*e2(:,2), ...
        -(e2(:,1).*e3(:,3)-e2(:,3).*e3(:,1)), ...
          e1(:,1).*e3(:,3)-e1(:,3).*e3(:,1), ...
        -(e1(:,1).*e2(:,3)-e1(:,3).*e2(:,1)), ...
          e2(:,1).*e3(:,2)-e2(:,2).*e3(:,1), ...
        -(e1(:,1).*e3(:,2)-e1(:,2).*e3(:,1)), ...
          e1(:,1).*e2(:,2)-e1(:,2).*e2(:,1), ...
        ];
    detA = e1(:,1).*cofA(:,1) + e2(:,1).*cofA(:,2)+ e3(:,1).*cofA(:,3);
    
    
    % compute gradients of basis functions
    dphi{1} =   cofA(:,[1 4 7])./repmat(detA,[1 3]);
    dphi{2} =   cofA(:,[2 5 8])./repmat(detA,[1 3]);
    dphi{3} =   cofA(:,[3 6 9])./repmat(detA,[1 3]);
    dphi{4} =   -dphi{1} - dphi{2} - dphi{3};

else

    % compute edges
    x1 =  Mesh.mfPi(xn,1);
    x2 =  Mesh.mfPi(xn,2);
    x3 =  Mesh.mfPi(xn,3);
    e1 =  x1 - x3;
    e2 =  x2 - x3;
    
    % compute gradients of basis functions
    dphi{1} =   [ e2(:,2) -e2(:,1)]./[2*vol,2*vol];
    dphi{2} =   [-e1(:,2)  e1(:,1)]./[2*vol,2*vol];
    dphi{3} = -dphi{1} - dphi{2};

end