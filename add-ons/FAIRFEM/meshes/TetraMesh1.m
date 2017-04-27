%==============================================================================
% This code is part of the Finite Element Method app for the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR/FAIRFEM 
%==============================================================================
%
% classdef TetraMesh1 < handle
%
% Finite Element Mesh based on tetrahedral subdivision of rectangular mesh.
%
% Each Cell is divided into 24 tetrahedra:
%
%
%  To construct an instance of this class type:
%
%  >> Mesh = TriMesh1(omega,m)
%
% Input:
% 	omega - description of spatial domain
%   m     - number of cells
%
% Properties:
%   xn     - node list
%   tri    - triangle list
%   dim    - space dimension
%   omega  - description of spatial domain
%   m      - number of cells
%   type   - type of partition
%   vol    - volume of triangles
%   nnodes - number of nodes
%   ntri   - number of triangles
%   dx1    - partial derivative operator
%   dx2    - partial derivative operator
%   dx3    - partial derivative operator
%   GRAD   - gradient operator
%   P1     - projection operator for node 1
%   P2     - projection operator for node 2
%   P3     - projection operator for node 3
%   P4     - projection operator for node 4
%   PC     - projection operator for Barycentrum
%   P      - prolongation operator
%   Pt     - restriction operator
%
%  Methods:
%   mfPu   - matrix free prolongation/restriction
%   mfPi   - matrix free edge projector
%   getP   - builds prolongation operator
%   tri2cc - averaging
%
%
% see also
% =========================================================================
classdef TetraMesh1 < handle
    
    properties
        % ===================================
        % node list
        % ===================================
        xn
        % ===================================
        % triangle list
        % ===================================
        tri
        % ===================================
        % space dimension
        % ===================================
        dim  = 3;
        % ===================================
        % description of computational domain
        % ===================================
        omega
        % ===================================
        % number of cells
        % ===================================
        m
        % ===================================
        % type of partition
        % ===================================
        type = 1;
        % ===================================
        % function handle to myself
        % ===================================
        me = @TetraMesh1;
        % ===================================
        % volume of triangles
        % ===================================
        vol
        % ===================================
        % number of nodes
        % ===================================
        nnodes
        % ===================================
        % number of triangles in mesh
        % ===================================
        ntri
    end
    
    properties (Access = public, Dependent) % These will be created when first callend and stores persistently
        % ===================================
        % dx1  - partial derivative operator
        % ===================================
        dx1
        % ===================================
        % dx2  - Partial derivative operator
        % ===================================
        dx2
        % ===================================
        % dx3  - Partial derivative operator
        % ===================================
        dx3
        % ===================================
        % GRAD - Gradient operator
        %
        %         | dx1 |
        %         |     |
        %  GRAD = | dx2 |
        %         |     |
        %         | dx3 |
        %
        % ===================================
        GRAD
        % ===================================
        % B - Vector gradient operator
        %
        %         | GRAD    0   0   |
        %         |                 |
        %   B =   |   0   GRAD  0   |
        %         |                 |
        %         |   0     0  GRAD |
        %
        % ===================================
        B
        % ===================================
        % P1 - Projection operator on Node 1
        % ===================================
        P1
        % ===================================
        % P2 - Projection operator on Node 2
        % ===================================
        P2
        % ===================================
        % P3 - Projection operator on Node 3
        % ===================================
        P3
        % ===================================
        % P4 - Projection operator on Node 4
        % ===================================
        P4
        % ===================================
        % PC - Projection operator on Barycenter
        % ===================================
        PC
        % ===================================
        % P  - Prolongation operator
        % ===================================
        P
        % ===================================
        % Pt - Restriction operator
        % ===================================
        Pt
        % ===================================
        % Boundary indices
        % ===================================
        boundaryIdx
        % ===================================
        % Boundary projector
        % ===================================
        boundaryProj
        mfdx1
        mfdx2
        mfdx3
        mfGRAD
    end
    
    properties (Access = private)
        % These are where the dependent data is actually stored
        dx1_
        dx2_
        dx3_
        GRAD_
        B_
        P1_
        P2_
        P3_
        P4_
        PC_
        P_
        Pt_
        mfdx1_
        mfdx2_
        mfdx3_
        mfGRAD_
        boundaryIdx_
        boundaryProj_
    end
    
    methods
        function this = TetraMesh1(omega,m)
            if nargin==0,
                help(mfilename);
                this.runMinimalExample;
                return;
            end
            h          = (omega(2:2:end)-omega(1:2:end))./m;
            this.omega = omega;
            this.m     = m;
            
            % get nodes (combine nodal, stg-1,stg-2,stg-3,cc grids)
            xc = @(i) linspace(omega(2*i-1)+h(i)/2,omega(2*i)-h(i)/2,m(i))'; % cell centers
            xn = @(i) linspace(omega(2*i-1),omega(2*i),m(i)+1)'; % nodes
            % add nodal grid
            nn = reshape(1:prod(m+1),m+1);
            nn = reshape(nn(1:end-1,1:end-1,1:end-1),1,[]); % get indices of bottom left vertices
            this.xn = reshape(getNodalGrid(omega,m),[],3);
            % add stg-1 grid
            ns  = prod(ones(length(m),1)*m+eye(length(m)),2); %staggered
            ns1 = reshape(size(this.xn,1)+(1:ns(1)),[m(1)+1 m(2) m(3)]);
            ns1 = reshape(ns1(1:end-1,:,:),1,[]); % get indices of bottom left vertices
            [x1,x2,x3] = ndgrid(xn(1),xc(2),xc(3));
            this.xn = [this.xn; x1(:), x2(:), x3(:)];
            % add stg-2 grid
            ns2 = reshape(size(this.xn,1)+(1:ns(2)),[m(1) m(2)+1 m(3)]);
            ns2 = reshape(ns2(:,1:end-1,:),1,[]); % get indices of bottom left vertices
            [x1,x2,x3] = ndgrid(xc(1),xn(2),xc(3));
            this.xn = [this.xn; x1(:), x2(:), x3(:)];
            % add stg-3 grid
            ns3 = reshape(size(this.xn,1)+(1:ns(3)),[m(1) m(2) m(3)+1]);
            ns3 = reshape(ns3(:,:,1:end-1),1,[]); % get indices of bottom left vertices
            [x1,x2,x3] = ndgrid(xc(1),xc(2),xn(3));
            this.xn = [this.xn; x1(:), x2(:), x3(:)];
            % add cc grid
            nc = reshape(size(this.xn,1)+(1:prod(m)),m);
            nc = reshape(nc,1,[]); % get indices of bottom left vertices
            this.xn = [this.xn; reshape(getCellCenteredGrid(omega,m),[],3)];
            
            % specify triangles
            iyn = m(1)+1; izn = prod(m(1:2)+1); izc = prod(m(1:2));
            this.tri   = [
                nn;           nn+1;         ns2;      nc;
                nn+1;         nn+izn+1;     ns2;      nc;
                nn+izn+1;     nn+izn;       ns2;      nc;
                nn+izn;       nn;           ns2;      nc; % 4
                nn+1;         nn+iyn+1;     ns1+1;    nc;
                nn+iyn+1;     nn+iyn+1+izn; ns1+1;    nc;
                nn+iyn+1+izn; nn+1+izn;     ns1+1;    nc;
                nn+1+izn;     nn+1;         ns1+1;    nc; % 8
                nn+1;         nn;           ns3;      nc;
                nn;           nn+iyn;       ns3;      nc;
                nn+iyn;       nn+1+iyn;     ns3;      nc;
                nn+1+iyn;     nn+1;         ns3;      nc; % 12                
                
                nn+iyn;       nn;           ns1;      nc;
                nn;           nn+izn;       ns1;      nc;
                nn+izn;       nn+iyn+izn;   ns1;      nc;
                nn+iyn+izn;   nn+iyn;       ns1;      nc; % 16
                
                nn+iyn+1;     nn+iyn;       ns2+m(1); nc;
                nn+iyn;       nn+iyn+izn;   ns2+m(1); nc;
                nn+iyn+izn;   nn+iyn+1+izn; ns2+m(1); nc;
                nn+iyn+1+izn; nn+iyn+1;     ns2+m(1); nc; % 20
                nn+izn;       nn+izn+1;     ns3+izc;  nc;
                nn+izn+1;     nn+izn+1+iyn; ns3+izc;  nc;
                nn+iyn+1+izn; nn+iyn+izn;   ns3+izc;  nc;
                nn+iyn+izn;   nn+izn;       ns3+izc;  nc; % 24
                ];
            this.tri   = reshape(this.tri,4,[])';
            this.nnodes = size(this.xn,1);
            this.ntri  = size(this.tri,1);
            this.vol = prod((omega(2:2:end)-omega(1:2:end))./m)/24*ones(this.ntri,1);
        end
        
        function runMinimalExample(~)
            omega = [0 4 2 6 0 3]; m = [3 4 6];
            Mesh  = feval(mfilename,omega,m);
        end
        
        function x = mfPi(this,x,i)
            % =============================================================
            % function x = mfPi(this,x,i)
            %
            % matrix free edge projector
            % =============================================================
            switch i
                case 1
                    P = this.P1;
                case 2
                    P = this.P2;
                case 3
                    P = this.P3;
                case 4
                    P = this.P4;
                case 'C'
                    P = this.PC;
            end
            if size(x,1) == this.ntri,
                % ajoint
                x = P'*x;
            else
                x = P * x;
            end
        end
        
        function x = tri2cc(this,x)
            % =============================================================
            % function x = tri2cc(x)
            %
            % averaging or adjoint
            % =============================================================
            if numel(x)==this.ntri,
                x = mean(reshape(x,24,[]),1);
            else
                x = reshape(x,1,[]);
                x = (1/24)*repmat(x,[24 1]);
                x = x(:);
            end
            
        end
        
        
        
        function Pu = mfPuNodal(this,yn,m,flag)
            % =============================================================
            % function Pu = mfPuNodal(~,yn,m,flag)
            %
            % matrix free prolongation/restriction for nodal quantities
            % =============================================================
            yn = reshape(yn,[],this.dim);
            v = @(x) x(:);
            switch flag
                case 'Pu' % coarse --> fine
                    mf = 2*m;
                    
                    % indices of coarse grid nodes (split by nodal, stg-i, c)
                    ns   = prod(ones(length(m),1)*m+eye(length(m)),2); %staggered
                    inc  = reshape(1:prod(m+1),m+1);
                    is1c = reshape(prod(m+1)+(1:ns(1)),[m(1)+1 m(2) m(3)]);
                    is2c = reshape(prod(m+1)+ns(1)+ (1:ns(2)),[m(1) m(2)+1 m(3)]);
                    is3c = reshape(prod(m+1)+ns(1)+ns(2)+(1:ns(3)),[m(1) m(2) m(3)+1]);
                    icc  = reshape(prod(m+1)+ns(1)+ns(2)+ns(3)+(1:prod(m)),m);
                    % indices of fine grid nodes
                    ns   = prod(ones(length(mf),1)*mf+eye(length(mf)),2); %staggered
                    inf  = reshape(1:prod(mf+1),mf+1);
                    is1f = reshape(prod(mf+1)+(1:ns(1)),[mf(1)+1 mf(2) mf(3)]);
                    is2f = reshape(prod(mf+1)+ns(1)+ (1:ns(2)),[mf(1) mf(2)+1 mf(3)]);
                    is3f = reshape(prod(mf+1)+ns(1)+ns(2)+(1:ns(3)),[mf(1) mf(2) mf(3)+1]);
                    icf  = reshape(prod(mf+1)+ns(1)+ns(2)+ns(3)+(1:prod(mf)),mf);
                    
                    % allocate space
                    Pu = zeros(icf(end),this.dim);
                    
                    % include existing nodes  (nodal-->nodal)
                    Pu(v(inf(1:2:end,1:2:end,1:2:end)),:)  = yn(inc(:),:);
                    % include existing nodes  (stg-1-->nodal)
                    Pu(v(inf(1:2:end,2:2:end,2:2:end)),:) = yn(is1c(:),:);
                    % include existing nodes  (stg-2-->nodal)
                    Pu(v(inf(2:2:end,1:2:end,2:2:end)),:) = yn(is2c(:),:);
                    % include existing nodes  (stg-3-->nodal)
                    Pu(v(inf(2:2:end,2:2:end,1:2:end)),:) = yn(is3c(:),:);
                    % include existing nodes  (cc-->nodal)
                    Pu(v(inf(2:2:end,2:2:end,2:2:end)),:) = yn(icc(:),:);
                    
                    % average to get edge-stg-1
                    Pu(v(inf(2:2:end,1:2:end,1:2:end)),:) = ...
                        .5*(yn(v(inc(1:end-1,:,:)),:) + yn(v(inc(2:end,:,:)),:)) ;
                    % average to get edge-stg-2
                    Pu(v(inf(1:2:end,2:2:end,1:2:end)),:) = ...
                        .5*(yn(v(inc(:,1:end-1,:)),:) + yn(v(inc(:,2:end,:)),:)) ;
                    % average to get edge-stg-3
                    Pu(v(inf(1:2:end,1:2:end,2:2:end)),:) = ...
                        .5*(yn(v(inc(:,:,1:end-1)),:) + yn(v(inc(:,:,2:end)),:)) ;
                    % get face-stg-1
                    Pu(is1f(:),:) = .25*(   Pu(v(inf(1:end,1:end-1,1:end-1)),:) ...
                        + Pu(v(inf(1:end,2:end  ,1:end-1)),:) ...
                        + Pu(v(inf(1:end,2:end  ,2:end  )),:) ...
                        + Pu(v(inf(1:end,1:end-1,2:end  )),:));
                    
                    
                    % get face-stg-2
                    Pu(is2f(:),:) = .25*(   Pu(v(inf(1:end-1,1:end  ,1:end-1)),:) ...
                        + Pu(v(inf(2:end  ,1:end  ,1:end-1)),:) ...
                        + Pu(v(inf(2:end  ,1:end  ,2:end  )),:) ...
                        + Pu(v(inf(1:end-1,1:end  ,2:end  )),:));
                    % get face-stg-3
                    Pu(is3f(:),:) = .25*(   Pu(v(inf(1:end-1,1:end-1,1:end)),:) ...
                        + Pu(v(inf(2:end  ,1:end-1,1:end)),:) ...
                        + Pu(v(inf(2:end  ,2:end  ,1:end)),:) ...
                        + Pu(v(inf(1:end-1,2:end  ,1:end)),:));
                    
                    
                    % get new cell-centers
                    Pu(icf(:),:)  = .125*(   Pu(v(inf(1:end-1,1:end-1,1:end-1)),:) ...
                        + Pu(v(inf(2:end  ,1:end-1,1:end-1)),:) ...
                        + Pu(v(inf(2:end  ,2:end  ,1:end-1)),:) ...
                        + Pu(v(inf(1:end-1,2:end  ,1:end-1)),:) ...
                        + Pu(v(inf(1:end-1,1:end-1,2:end)),:) ...
                        + Pu(v(inf(2:end  ,1:end-1,2:end)),:) ...
                        + Pu(v(inf(2:end  ,2:end  ,2:end)),:) ...
                        + Pu(v(inf(1:end-1,2:end  ,2:end)),:));
                    
                    
                case 'PTu' % fine --> coarse
                    % include parent nodes
                    mf = m;
                    m  = m/2;
                    % indices of coarse grid nodes (split by nodal, stg-i, c)
                    ns   = prod(ones(length(m),1)*m+eye(length(m)),2); %staggered
                    inc  = reshape(1:prod(m+1),m+1);
                    is1c = reshape(prod(m+1)+(1:ns(1)),[m(1)+1 m(2) m(3)]);
                    is2c = reshape(prod(m+1)+ns(1)+ (1:ns(2)),[m(1) m(2)+1 m(3)]);
                    is3c = reshape(prod(m+1)+ns(1)+ns(2)+(1:ns(3)),[m(1) m(2) m(3)+1]);
                    icc  = reshape(prod(m+1)+ns(1)+ns(2)+ns(3)+(1:prod(m)),m);
                    % indices of fine grid nodes
                    ns   = prod(ones(length(mf),1)*mf+eye(length(mf)),2); %staggered
                    inf  = reshape(1:prod(mf+1),mf+1);
                    is1f = reshape(prod(mf+1)+(1:ns(1)),[mf(1)+1 mf(2) mf(3)]);
                    is2f = reshape(prod(mf+1)+ns(1)+ (1:ns(2)),[mf(1) mf(2)+1 mf(3)]);
                    is3f = reshape(prod(mf+1)+ns(1)+ns(2)+(1:ns(3)),[mf(1) mf(2) mf(3)+1]);
                    icf  = reshape(prod(mf+1)+ns(1)+ns(2)+ns(3)+(1:prod(mf)),mf);
                    
                    % allocate space
                    Pu = zeros(icc(end),this.dim);
                    
                    % push weights from cell-centers
                    yn(v(inf(1:end-1,1:end-1,1:end-1)),:)  = yn(v(inf(1:end-1,1:end-1,1:end-1)),:) + .125* yn(icf(:),:);
                    yn(v(inf(2:end  ,1:end-1,1:end-1)),:)  = yn(v(inf(2:end  ,1:end-1,1:end-1)),:) + .125* yn(icf(:),:);
                    yn(v(inf(2:end  ,2:end  ,1:end-1)),:)  = yn(v(inf(2:end  ,2:end  ,1:end-1)),:) + .125* yn(icf(:),:);
                    yn(v(inf(1:end-1,2:end  ,1:end-1)),:)  = yn(v(inf(1:end-1,2:end  ,1:end-1)),:) + .125* yn(icf(:),:);
                    yn(v(inf(1:end-1,1:end-1,2:end)),:)    = yn(v(inf(1:end-1,1:end-1,2:end)),:)   + .125* yn(icf(:),:);
                    yn(v(inf(2:end  ,1:end-1,2:end)),:)    = yn(v(inf(2:end  ,1:end-1,2:end)),:)   + .125* yn(icf(:),:);
                    yn(v(inf(2:end  ,2:end  ,2:end)),:)    = yn(v(inf(2:end  ,2:end  ,2:end)),:)   + .125* yn(icf(:),:);
                    yn(v(inf(1:end-1,2:end  ,2:end)),:)    = yn(v(inf(1:end-1,2:end  ,2:end)),:)   + .125* yn(icf(:),:);
                    % push weights from face-stg-1
                    yn(v(inf(1:end,1:end-1,1:end-1)),:)   = yn(v(inf(1:end,1:end-1,1:end-1)),:) + .25 * yn(is1f(:),:);
                    yn(v(inf(1:end,2:end  ,1:end-1)),:)   = yn(v(inf(1:end,2:end  ,1:end-1)),:) + .25 * yn(is1f(:),:);
                    yn(v(inf(1:end,2:end  ,2:end  )),:)   = yn(v(inf(1:end,2:end  ,2:end  )),:) + .25 * yn(is1f(:),:);
                    yn(v(inf(1:end,1:end-1,2:end  )),:)   = yn(v(inf(1:end,1:end-1,2:end  )),:) + .25 * yn(is1f(:),:);
                    % push weights from face-stg-2
                    yn(v(inf(1:end-1,1:end  ,1:end-1)),:) = yn(v(inf(1:end-1,1:end  ,1:end-1)),:) + .25 * yn(is2f(:),:);
                    yn(v(inf(2:end  ,1:end  ,1:end-1)),:) = yn(v(inf(2:end  ,1:end  ,1:end-1)),:) + .25 * yn(is2f(:),:);
                    yn(v(inf(2:end  ,1:end  ,2:end  )),:) = yn(v(inf(2:end  ,1:end  ,2:end  )),:) + .25 * yn(is2f(:),:);
                    yn(v(inf(1:end-1,1:end  ,2:end  )),:) = yn(v(inf(1:end-1,1:end  ,2:end  )),:) + .25 * yn(is2f(:),:);
                    % push weights from face-stg-3
                    yn(v(inf(1:end-1,1:end-1,1:end)),:)   = yn(v(inf(1:end-1,1:end-1,1:end)),:)  + .25 *  yn(is3f(:),:);
                    yn(v(inf(2:end  ,1:end-1,1:end)),:)   = yn(v(inf(2:end  ,1:end-1,1:end)),:)  + .25 *  yn(is3f(:),:);
                    yn(v(inf(2:end  ,2:end  ,1:end)),:)   = yn(v(inf(2:end  ,2:end  ,1:end)),:)  + .25 *  yn(is3f(:),:);
                    yn(v(inf(1:end-1,2:end  ,1:end)),:)   = yn(v(inf(1:end-1,2:end  ,1:end)),:)  + .25 *  yn(is3f(:),:);
                    
                    % include existing nodes  (nodal-->nodal)
                    Pu(inc(:),:)  = yn(v(inf(1:2:end,1:2:end,1:2:end)),:) ;
                    % include existing nodes  (stg-1-->nodal)
                    Pu(is1c(:),:) = yn(v(inf(1:2:end,2:2:end,2:2:end)),:);
                    % include existing nodes  (stg-2-->nodal)
                    Pu(is2c(:),:) = yn(v(inf(2:2:end,1:2:end,2:2:end)),:);
                    % include existing nodes  (stg-3-->nodal)
                    Pu(is3c(:),:) = yn(v(inf(2:2:end,2:2:end,1:2:end)),:);
                    % include existing nodes  (cc-->nodal)
                    Pu(icc(:),:)  = yn(v(inf(2:2:end,2:2:end,2:2:end)),:);
                    
                    % average to get edge-stg-1
                    Pu(v(inc(1:end-1,:,:)),:) = Pu(v(inc(1:end-1,:,:)),:) + .5 * yn(v(inf(2:2:end,1:2:end,1:2:end)),:);
                    Pu(v(inc(2:end,:,:)),:)   = Pu(v(inc(2:end,:,:)),:)   + .5 * yn(v(inf(2:2:end,1:2:end,1:2:end)),:);
                    
                    % average to get edge-stg-2
                    Pu(v(inc(:,1:end-1,:)),:) = Pu(v(inc(:,1:end-1,:)),:) + .5 * yn(v(inf(1:2:end,2:2:end,1:2:end)),:);
                    Pu(v(inc(:,2:end  ,:)),:) = Pu(v(inc(:,2:end  ,:)),:) + .5 * yn(v(inf(1:2:end,2:2:end,1:2:end)),:);
                    
                    % average to get edge-stg-3
                    Pu(v(inc(:,:,1:end-1)),:) = Pu(v(inc(:,:,1:end-1)),:) + .5 * yn(v(inf(1:2:end,1:2:end,2:2:end)),:);
                    Pu(v(inc(:,:,2:end  )),:) = Pu(v(inc(:,:,2:end  )),:) + .5 * yn(v(inf(1:2:end,1:2:end,2:2:end)),:);
                    
            end
            Pu = Pu(:);
            
        end
        function Pu = mfInterNodal(this,yn,m)
            % =============================================================
            % function Pu = mfInterNodal(~,yn,m,)
            %
            % interpolates to nodal grid of a finer resolution.
            % =============================================================
            yn = reshape(yn,[],this.dim);
            v = @(x) x(:);
            mf = 2*m;
            
            % indices of coarse grid nodes (split by nodal, stg-i, c)
            ns   = prod(ones(length(m),1)*m+eye(length(m)),2); %staggered
            inc  = reshape(1:prod(m+1),m+1);
            is1c = reshape(prod(m+1)+(1:ns(1)),[m(1)+1 m(2) m(3)]);
            is2c = reshape(prod(m+1)+ns(1)+ (1:ns(2)),[m(1) m(2)+1 m(3)]);
            is3c = reshape(prod(m+1)+ns(1)+ns(2)+(1:ns(3)),[m(1) m(2) m(3)+1]);
            icc  = reshape(prod(m+1)+ns(1)+ns(2)+ns(3)+(1:prod(m)),m);
            % indices of fine grid nodes
            inf  = reshape(1:prod(mf+1),mf+1);
            
            % allocate space
            Pu = zeros(inf(end),this.dim);
            
            % include existing nodes  (nodal-->nodal)
            Pu(v(inf(1:2:end,1:2:end,1:2:end)),:)  = yn(inc(:),:);
            % include existing nodes  (stg-1-->nodal)
            Pu(v(inf(1:2:end,2:2:end,2:2:end)),:) = yn(is1c(:),:);
            % include existing nodes  (stg-2-->nodal)
            Pu(v(inf(2:2:end,1:2:end,2:2:end)),:) = yn(is2c(:),:);
            % include existing nodes  (stg-3-->nodal)
            Pu(v(inf(2:2:end,2:2:end,1:2:end)),:) = yn(is3c(:),:);
            % include existing nodes  (cc-->nodal)
            Pu(v(inf(2:2:end,2:2:end,2:2:end)),:) = yn(icc(:),:);
            
            % average to get edge-stg-1
            Pu(v(inf(2:2:end,1:2:end,1:2:end)),:) = ...
                .5*(yn(v(inc(1:end-1,:,:)),:) + yn(v(inc(2:end,:,:)),:)) ;
            % average to get edge-stg-2
            Pu(v(inf(1:2:end,2:2:end,1:2:end)),:) = ...
                .5*(yn(v(inc(:,1:end-1,:)),:) + yn(v(inc(:,2:end,:)),:)) ;
            % average to get edge-stg-3
            Pu(v(inf(1:2:end,1:2:end,2:2:end)),:) = ...
                .5*(yn(v(inc(:,:,1:end-1)),:) + yn(v(inc(:,:,2:end)),:)) ;
            
            
            Pu = Pu(:);
            
        end
        
        
        function P = getPuNodal(~,m)
            % =============================================================
            % function P = getPuNodal(~,m)
            %
            % returns prolongation operator for input of cell-width m
            % =============================================================
            
            v = @(x) x(:);
            
            mf = 2*m;
            % indices of coarse grid nodes (split by nodal, stg-i, c)
            ns   = prod(ones(length(m),1)*m+eye(length(m)),2); %staggered
            inc  = reshape(1:prod(m+1),m+1);
            is1c = reshape(prod(m+1)+(1:ns(1)),[m(1)+1 m(2) m(3)]);
            is2c = reshape(prod(m+1)+ns(1)+ (1:ns(2)),[m(1) m(2)+1 m(3)]);
            is3c = reshape(prod(m+1)+ns(1)+ns(2)+(1:ns(3)),[m(1) m(2) m(3)+1]);
            icc  = reshape(prod(m+1)+ns(1)+ns(2)+ns(3)+(1:prod(m)),m);
            % indices of fine grid nodes
            ns   = prod(ones(length(mf),1)*mf+eye(length(mf)),2); %staggered
            inf  = reshape(1:prod(mf+1),mf+1);
            is1f = reshape(prod(mf+1)+(1:ns(1)),[mf(1)+1 mf(2) mf(3)]);
            is2f = reshape(prod(mf+1)+ns(1)+ (1:ns(2)),[mf(1) mf(2)+1 mf(3)]);
            is3f = reshape(prod(mf+1)+ns(1)+ns(2)+(1:ns(3)),[mf(1) mf(2) mf(3)+1]);
            icf  = reshape(prod(mf+1)+ns(1)+ns(2)+ns(3)+(1:prod(mf)),mf);
            
            % allocate space
            I = []; J = []; W = [];
            
            % include existing nodes  (nodal-->nodal)
            ii = inf(1:2:end,1:2:end,1:2:end);
            jj = inc;
            ww = ones(size(jj));
            I = [I; ii(:)]; J = [J; jj(:)]; W = [W; ww(:)];
            % include existing nodes  (stg-1-->nodal)
            ii = inf(1:2:end,2:2:end,2:2:end);
            jj = is1c;
            ww = ones(size(jj));
            I = [I; ii(:)]; J = [J; jj(:)]; W = [W; ww(:)];
            % include existing nodes  (stg-2-->nodal)
            ii = inf(2:2:end,1:2:end,2:2:end);
            jj = is2c;
            ww = ones(size(jj));
            I = [I; ii(:)]; J = [J; jj(:)]; W = [W; ww(:)];
            % include existing nodes  (stg-3-->nodal)
            ii = inf(2:2:end,2:2:end,1:2:end);
            jj = is3c;
            ww = ones(size(jj));
            I = [I; ii(:)]; J = [J; jj(:)]; W = [W; ww(:)];
            % include existing nodes  (cc-->nodal)
            ii = inf(2:2:end,2:2:end,2:2:end);
            jj = icc;
            ww = ones(size(jj));
            I = [I; ii(:)]; J = [J; jj(:)]; W = [W; ww(:)];
            
            % average to get edge-stg-1
            ii   = inf(2:2:end,1:2:end,1:2:end);
            jj   = [reshape(inc(1:end-1,:,:),[],1); reshape(inc(2:end,:,:),[],1)];
            ww   = .5*ones(size(jj));
            I    = [I; ii(:);ii(:)]; J = [J; jj(:)]; W = [W; ww(:)];
            % average to get edge-stg-2
            ii   = inf(1:2:end,2:2:end,1:2:end);
            jj   = [reshape(inc(:,1:end-1,:),[],1); reshape(inc(:,2:end,:),[],1)];
            ww   = .5*ones(size(jj));
            I    = [I; ii(:);ii(:)]; J = [J; jj(:)]; W = [W; ww(:)];
            % average to get edge-stg-3
            ii   = inf(1:2:end,1:2:end,2:2:end);
            jj   = [reshape(inc(:,:,1:end-1),[],1); reshape(inc(:,:,2:end),[],1)];
            ww   = .5*ones(size(jj));
            I    = [I; ii(:);ii(:)]; J = [J; jj(:)]; W = [W; ww(:)];
            
            % coarse --> fine(nodal only)
            P1  = sparse(I,J,W,prod(mf+1),icc(end));
            % allocate space
            I = []; J = []; W = [];
            
            % get face-stg-1
            ii  = is1f;
            jj  = [
                v(inf(1:end,1:end-1,1:end-1));
                v(inf(1:end,2:end  ,1:end-1));
                v(inf(1:end,2:end  ,2:end  ));
                v(inf(1:end,1:end-1,2:end  ));
                ];
            ww   = .25*ones(size(jj));
            I    = [I; ii(:);ii(:);ii(:);ii(:)]; J = [J; jj(:)]; W = [W; ww(:)];
            % get face-stg-2
            ii  = is2f;
            jj  = [
                v(inf(1:end-1,1:end  ,1:end-1));
                v(inf(2:end  ,1:end  ,1:end-1));
                v(inf(2:end  ,1:end  ,2:end  ));
                v(inf(1:end-1,1:end  ,2:end  ));
                ];
            ww   = .25*ones(size(jj));
            I    = [I; ii(:);ii(:);ii(:);ii(:)]; J = [J; jj(:)]; W = [W; ww(:)];
            % get face-stg-3
            ii  = is3f;
            jj  = [
                v(inf(1:end-1,1:end-1,1:end));
                v(inf(2:end  ,1:end-1,1:end));
                v(inf(2:end  ,2:end  ,1:end));
                v(inf(1:end-1,2:end  ,1:end));
                ];
            ww   = .25*ones(size(jj));
            I    = [I; ii(:);ii(:);ii(:);ii(:)]; J = [J; jj(:)]; W = [W; ww(:)];
            
            
            % get new cell-centers
            ii  = icf;
            jj  = [
                v(inf(1:end-1,1:end-1,1:end-1));
                v(inf(2:end  ,1:end-1,1:end-1));
                v(inf(2:end  ,2:end  ,1:end-1));
                v(inf(1:end-1,2:end  ,1:end-1));
                v(inf(1:end-1,1:end-1,2:end));
                v(inf(2:end  ,1:end-1,2:end));
                v(inf(2:end  ,2:end  ,2:end));
                v(inf(1:end-1,2:end  ,2:end));
                ];
            ww   = .125*ones(size(jj));
            I    = [I; repmat(ii(:),[8,1])]; J = [J; jj(:)]; W = [W; ww(:)];
            
            P2   = sparse(I-prod(mf+1),J,W,ns(1)+ns(2)+ns(3)+prod(mf),prod(mf+1));
            
            P    = [P1; P2*P1];
        end
        
        
        
        % ========== get methods ========================================
        function dx1 = get.dx1(this)
            if isempty(this.dx1_),
                [this.dx1_, this.dx2_, this.dx3_] = getGradientMatrixFEM(this,0);
            end
            dx1 = this.dx1_;
        end
        
        function dx2 = get.dx2(this)
            if isempty(this.dx2_),
                [this.dx1_, this.dx2_, this.dx3_] = getGradientMatrixFEM(this,0);
            end
            dx2 = this.dx2_;
        end
        
        function dx3 = get.dx3(this)
            if isempty(this.dx3_),
                [this.dx1_, this.dx2_, this.dx3_] = getGradientMatrixFEM(this,0);
            end
            dx3 = this.dx3_;
        end
        
        function GRAD = get.GRAD(this)
            if isempty(this.GRAD_),
                this.GRAD_ = [this.dx1;this.dx2;this.dx3];
            end
            GRAD = this.GRAD_;
        end
        function B = get.B(this)
            if isempty(this.B_),
                this.B_ = blkdiag(this.GRAD,this.GRAD,this.GRAD);
            end
            B = this.B_;
        end
        
        function P1 = get.P1(this)
            if isempty(this.P1_),
                A = speye(this.nnodes);
                this.P1_ = A(this.tri(:,1),:);
            end
            P1 = this.P1_;
        end
        function P2 = get.P2(this)
            if isempty(this.P2_),
                A = speye(this.nnodes);
                this.P2_ = A(this.tri(:,2),:);
            end
            P2 = this.P2_;
        end
        
        function P3 = get.P3(this)
            if isempty(this.P3_),
                A = speye(this.nnodes);
                this.P3_ = A(this.tri(:,3),:);
            end
            P3 = this.P3_;
        end
        
        function P4 = get.P4(this)
            if isempty(this.P4_),
                A = speye(this.nnodes);
                this.P4_ = A(this.tri(:,4),:);
            end
            P4 = this.P4_;
        end
        function PC = get.PC(this)
            if isempty(this.PC_),
                this.PC_ = (this.P1+this.P2+this.P3+this.P4)/4;
            end
            PC = this.PC_;
        end
        function P = get.P(this)
            if isempty(this.P_),
                this.P_ = this.getPuNodal(this.m);
            end
            P = this.P_;
        end
        function Pt = get.Pt(this)
            if isempty(this.Pt_),
                this.Pt_ = this.getPuNodal(this.m/2);
            end
            Pt = this.Pt_;
        end
        
        function mfdx1 = get.mfdx1(this)
            if isempty(this.mfdx1_),
                [this.mfdx1_, this.mfdx2_, this.mfdx3_] = getGradientMatrixFEM(this,1);
            end
            mfdx1 = this.mfdx1_;
        end
        
        function mfdx2 = get.mfdx2(this)
            if isempty(this.mfdx2_),
                [this.mfdx1_, this.mfdx2_, this.mfdx3_] = getGradientMatrixFEM(this,1);
            end
            mfdx2 = this.mfdx2_;
        end
        
        function mfdx3 = get.mfdx3(this)
            if isempty(this.mfdx3_),
                [this.mfdx1_, this.mfdx2_, this.mfdx3_] = getGradientMatrixFEM(this,1);
            end
            mfdx3 = this.mfdx3_;
        end
        
        function mfGRAD = get.mfGRAD(this)
            if isempty(this.mfGRAD_),
                this.mfGRAD_ = getGradientMatrixFEM(this,1);
            end
            mfGRAD = this.mfGRAD_;
        end
        
        function idx = get.boundaryIdx(this)
            if isempty(this.boundaryIdx_),
                mkvc = @(v) v(:);
                ns   = prod(ones(length(this.m),1)*this.m+eye(length(this.m)),2); %staggered
                           
                % nodal points
                id = reshape(1:prod(this.m+1),this.m+1);
                idx = [ ...
                            mkvc(id([1,end],:,:));
                            mkvc(id(:,[1,end],:));
                            mkvc(id(:,:,[1,end]));
                            ];
                % stg-1 points
                id = prod(this.m+1) + reshape(1:prod(this.m+[1,0,0]),this.m+[1,0,0]);
                idx = [idx; mkvc(id([1,end],:,:)) ];
                % stg-2 points
                id = prod(this.m+1) +ns(1) + reshape(1:prod(this.m+[0,1,0]),this.m+[0,1,0]);
                idx = [idx; mkvc(id(:,[1,end],:)) ];
                % stg-3 points
                id = prod(this.m+1) +ns(1) + ns(2) + reshape(1:prod(this.m+[0,0,1]),this.m+[0,0,1]);
                idx = [idx; mkvc(id(:,:,[1,end])) ];
                
                this.boundaryIdx_ = unique( idx);
            end
            idx = this.boundaryIdx_;
        end
        
        function idx = get.boundaryProj(this)
            if isempty(this.boundaryProj_),
                idx = this.boundaryIdx;
                
                P = speye(size(this.xn,1));
                P = P(idx,:);
                P = kron(speye(3),P);
                
                this.boundaryProj_ = P;
            end
            idx = this.boundaryProj_;
        end 
        
    end
end

