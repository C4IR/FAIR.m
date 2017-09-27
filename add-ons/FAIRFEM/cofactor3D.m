% For details see 
% - https://github.com/C4IR/FAIRFEM 
%==============================================================================
%
% function [cof,dCof] = cofactor3D(By,Mesh,doDerivative,matrixFree)
%
% cofactor matrix, and derivative
%
% Input:
%
% By            - partial derivative of y wrt x
% Mesh          - description of mesh, struct
% 
% Output:
%
% [cof, dcof] - cofactor and partial derivative operators
%
%==============================================================================

function [cof,dCof] = cofactor3D(By,Mesh,doDerivative,matrixFree)
     
if ~matrixFree
    
    dCof = [];
    D = @(i,j) By(:,(j-1)*3+i); % gets partial_i y^j
    cof = [
        D(2,2).*D(3,3)-D(3,2).*D(2,3),...%ok
       -D(1,2).*D(3,3)+D(3,2).*D(1,3),...%trouble
        D(1,2).*D(2,3)-D(2,2).*D(1,3),...%ok
       -D(2,1).*D(3,3)+D(3,1).*D(2,3),...%ok
        D(1,1).*D(3,3)-D(3,1).*D(1,3),...%ok
       -D(1,1).*D(2,3)+D(2,1).*D(1,3),...%ok
        D(2,1).*D(3,2)-D(3,1).*D(2,2),...%ok
       -D(1,1).*D(3,2)+D(3,1).*D(1,2),...%ok
        D(1,1).*D(2,2)-D(2,1).*D(1,2)... %ok
        ];

    if doDerivative
        D1 = Mesh.dx1; D2 = Mesh.dx2; D3 = Mesh.dx3;
        Z = sparse(size(D1,1), size(D2,2));
        D = @(i,j) sdiag(By(:,(j-1)*3+i));

        dCof{1} = [  Z                 , D(3,3)*D2-D(2,3)*D3, D(2,2)*D3-D(3,2)*D2];
        dCof{2} = [  Z                 , D(1,3)*D3-D(3,3)*D1, D(3,2)*D1-D(1,2)*D3];
        dCof{3} = [  Z                 , D(2,3)*D1-D(1,3)*D2, D(1,2)*D2-D(2,2)*D1];
        dCof{4} = [ D(2,3)*D3-D(3,3)*D2, Z                  , D(3,1)*D2-D(2,1)*D3];
        dCof{5} = [ D(3,3)*D1-D(1,3)*D3, Z                  , D(1,1)*D3-D(3,1)*D1];
        dCof{6} = [ D(1,3)*D2-D(2,3)*D1, Z                  , D(2,1)*D1-D(1,1)*D2];
        dCof{7} = [ D(3,2)*D2-D(2,2)*D3, D(2,1)*D3-D(3,1)*D2, Z];
        dCof{8} = [ D(1,2)*D3-D(3,2)*D1, D(3,1)*D1-D(1,1)*D3, Z];
        dCof{9} = [ D(2,2)*D1-D(1,2)*D2, D(1,1)*D2-D(2,1)*D1, Z];
    end

else

    dCof = [];
    dx1 = Mesh.mfdx1; dx2 = Mesh.mfdx2;  dx3 = Mesh.mfdx3;
            
    cof{1} = @(yc)  dx2.D(yc(:,2)).*dx3.D(yc(:,3)) - dx3.D(yc(:,2)).*dx2.D(yc(:,3));
    cof{2} = @(yc) -dx1.D(yc(:,2)).*dx3.D(yc(:,3)) + dx3.D(yc(:,2)).*dx1.D(yc(:,3));
    cof{3} = @(yc)  dx1.D(yc(:,2)).*dx2.D(yc(:,3)) - dx2.D(yc(:,2)).*dx1.D(yc(:,3));
    cof{4} = @(yc) -dx2.D(yc(:,1)).*dx3.D(yc(:,3)) + dx3.D(yc(:,1)).*dx2.D(yc(:,3));
    cof{5} = @(yc)  dx1.D(yc(:,1)).*dx3.D(yc(:,3)) - dx3.D(yc(:,1)).*dx1.D(yc(:,3));
    cof{6} = @(yc) -dx1.D(yc(:,1)).*dx2.D(yc(:,3)) + dx2.D(yc(:,1)).*dx1.D(yc(:,3));
    cof{7} = @(yc)  dx2.D(yc(:,1)).*dx3.D(yc(:,2)) - dx3.D(yc(:,1)).*dx2.D(yc(:,2));
    cof{8} = @(yc) -dx1.D(yc(:,1)).*dx3.D(yc(:,2)) + dx3.D(yc(:,1)).*dx1.D(yc(:,2));
    cof{9} = @(yc)  dx1.D(yc(:,1)).*dx2.D(yc(:,2)) - dx2.D(yc(:,1)).*dx1.D(yc(:,2));
    
    if doDerivative

        dx{1}.D = Mesh.mfdx1.D; dx{1}.Dadj = Mesh.mfdx1.Dadj;
        dx{2}.D = Mesh.mfdx2.D; dx{2}.Dadj = Mesh.mfdx2.Dadj;
        dx{3}.D = Mesh.mfdx3.D; dx{3}.Dadj = Mesh.mfdx3.Dadj;

        % adjoint/transpose of an element of dCof matrix
        dCofadjele = @(i,j,k,x,yc) dx{i}.Dadj(dx{j}.D(yc(:,k)).*x);

		dCof.dCofadj{1} = @(x,yc) [zeros(Mesh.nnodes,1); dCofadjele(2,3,3,x,yc)-dCofadjele(3,2,3,x,yc); dCofadjele(3,2,2,x,yc)-dCofadjele(2,3,2,x,yc)];
		dCof.dCofadj{2} = @(x,yc) [zeros(Mesh.nnodes,1); dCofadjele(3,1,3,x,yc)-dCofadjele(1,3,3,x,yc); dCofadjele(1,3,2,x,yc)-dCofadjele(3,1,2,x,yc)];
		dCof.dCofadj{3} = @(x,yc) [zeros(Mesh.nnodes,1); dCofadjele(1,2,3,x,yc)-dCofadjele(2,1,3,x,yc); dCofadjele(2,1,2,x,yc)-dCofadjele(1,2,2,x,yc)];

		dCof.dCofadj{4} = @(x,yc) [dCofadjele(3,2,3,x,yc)-dCofadjele(2,3,3,x,yc); zeros(Mesh.nnodes,1); dCofadjele(2,3,1,x,yc)-dCofadjele(3,2,1,x,yc)];
		dCof.dCofadj{5} = @(x,yc) [dCofadjele(1,3,3,x,yc)-dCofadjele(3,1,3,x,yc); zeros(Mesh.nnodes,1); dCofadjele(3,1,1,x,yc)-dCofadjele(1,3,1,x,yc)];
		dCof.dCofadj{6} = @(x,yc) [dCofadjele(2,1,3,x,yc)-dCofadjele(1,2,3,x,yc); zeros(Mesh.nnodes,1); dCofadjele(1,2,1,x,yc)-dCofadjele(2,1,1,x,yc)];

		dCof.dCofadj{7} = @(x,yc) [dCofadjele(2,3,2,x,yc)-dCofadjele(3,2,2,x,yc); dCofadjele(3,2,1,x,yc)-dCofadjele(2,3,1,x,yc); zeros(Mesh.nnodes,1)];
		dCof.dCofadj{8} = @(x,yc) [dCofadjele(3,1,2,x,yc)-dCofadjele(1,3,2,x,yc); dCofadjele(1,3,1,x,yc)-dCofadjele(3,1,1,x,yc); zeros(Mesh.nnodes,1)];
		dCof.dCofadj{9} = @(x,yc) [dCofadjele(1,2,2,x,yc)-dCofadjele(2,1,2,x,yc); dCofadjele(2,1,1,x,yc)-dCofadjele(1,2,1,x,yc); zeros(Mesh.nnodes,1)];

		% element of dCof matrix
		dCofele = @(i,j,k,l,x,yc) dx{i}.D(yc(:,j)).*dx{k}.D(x(:,l));

		dCof.dCof{1} = @(x,yc) dCofele(3,3,2,2,x,yc)-dCofele(2,3,3,2,x,yc) + dCofele(2,2,3,3,x,yc) - dCofele(3,2,2,3,x,yc);
		dCof.dCof{2} = @(x,yc) dCofele(1,3,3,2,x,yc)-dCofele(3,3,1,2,x,yc) + dCofele(3,2,1,3,x,yc) - dCofele(1,2,3,3,x,yc);
		dCof.dCof{3} = @(x,yc) dCofele(2,3,1,2,x,yc)-dCofele(1,3,2,2,x,yc) + dCofele(1,2,2,3,x,yc) - dCofele(2,2,1,3,x,yc);

		dCof.dCof{4} = @(x,yc) dCofele(2,3,3,1,x,yc)-dCofele(3,3,2,1,x,yc) + dCofele(3,1,2,3,x,yc) - dCofele(2,1,3,3,x,yc);
		dCof.dCof{5} = @(x,yc) dCofele(3,3,1,1,x,yc)-dCofele(1,3,3,1,x,yc) + dCofele(1,1,3,3,x,yc) - dCofele(3,1,1,3,x,yc);
		dCof.dCof{6} = @(x,yc) dCofele(1,3,2,1,x,yc)-dCofele(2,3,1,1,x,yc) + dCofele(2,1,1,3,x,yc) - dCofele(1,1,2,3,x,yc);

		dCof.dCof{7} = @(x,yc) dCofele(3,2,2,1,x,yc)-dCofele(2,2,3,1,x,yc) + dCofele(2,1,3,2,x,yc) - dCofele(3,1,2,2,x,yc);
		dCof.dCof{8} = @(x,yc) dCofele(1,2,3,1,x,yc)-dCofele(3,2,1,1,x,yc) + dCofele(3,1,1,2,x,yc) - dCofele(1,1,3,2,x,yc);
		dCof.dCof{9} = @(x,yc) dCofele(2,2,1,1,x,yc)-dCofele(1,2,2,1,x,yc) + dCofele(1,1,2,2,x,yc) - dCofele(2,1,1,2,x,yc);
       
    end
    
end