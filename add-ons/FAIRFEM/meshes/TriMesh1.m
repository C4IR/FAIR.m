%==============================================================================
% This code is part of the Finite Element Method app for the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR/FAIRFEM 
%==============================================================================
%
% classdef TriMesh1 < handle
%
% Finite Element Mesh based on triangular subdivision of rectangular mesh.
%
% Each Cell is divided into two triangles:
%
%  o---------o---------o---------o---------o
%  |        /|        /|        /|        /|
%  |      /  |      /  |      /  |      /  |
%  |    /    |    /    |    /    |    /    |
%  |  /      |  /      |  /      |  /      |
%  o---------o---------o---------o---------o
%  |        /|        /|        /|        /|
%  |      /  |      /  |      /  |      /  |
%  |    /    |    /    |    /    |    /    |
%  |  /      |  /      |  /      |  /      |
%  o---------o---------o---------o---------o
%  |       / |       / |       / |       / |
%  |     /   |     /   |     /   |     /   |
%  |   /     |   /     |   /     |   /     |
%  | /       | /       | /       | /       |
%  o---------o---------o---------o---------o
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
%   GRAD   - gradient operator
%   P1     - projection operator for node 1
%   P2     - projection operator for node 2
%   P3     - projection operator for node 3
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
classdef TriMesh1 < handle
    
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
        dim  = 2;
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
        me = @TriMesh1
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
        % GRAD - Gradient operator
        %
        %         | dx1 |
        %  GRAD = |     |
        %         | dx2 |
        %
        % ===================================
        GRAD
        % ===================================
        % B - Vector gradient operator
        %
        %         | GRAD  0    |
        %  B =    |            |
        %         | 0     GRAD |
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
        mfGRAD
    end
    
    properties (Access = private)
        % These are where the dependent data is actually stored
        dx1_
        dx2_
        GRAD_
        B_
        P1_
        P2_
        P3_
        PC_
        P_
        Pt_
        mfdx1_
        mfdx2_
        mfGRAD_
        boundaryIdx_
        boundaryProj_
    end
    
    methods
        function this = TriMesh1(omega,m)
            if nargin==0,
                help(mfilename);
                this.runMinimalExample;
                return;
            end
            this.omega = omega;
            this.m     = m;
            
            this.xn = reshape(getNodalGrid(omega,m),[],2);
            % get indices of bottom left vertices
            nodes = reshape(1:prod(m+1),m+1);
            nodes = reshape(nodes(1:end-1,1:end-1),1,[]);
            % specify triangles
            this.tri   = [nodes+1; nodes+m(1)+2; nodes; nodes+m(1)+1; nodes; nodes+m(1)+2];
            this.tri   = reshape(this.tri,3,[])';
            this.nnodes = prod(m+1);
            this.ntri  = size(this.tri,1);
            this.vol = prod((omega(2:2:end)-omega(1:2:end))./m)/2*ones(this.ntri,1);
        end
        
        function runMinimalExample(~)
            omega = [0 4 2 6]; m = [8 16];
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
                x = mean(reshape(x,2,[]),1);
            else
                x = reshape(x,1,[]);
                x = .5*[x;x];
                x = x(:);
            end
            
        end
        

        
        function Pu = mfPuNodal(~,yn,m,flag)
            % =============================================================
            % function Pu = mfPuNodal(~,yn,m,flag)
            %
            % matrix free prolongation/restriction for nodal quantities
            % =============================================================
            d  = numel(yn)/prod(m+1);
            yn = reshape(yn,[m+1, d]);
            switch flag
                case 'Pu' % coarse --> fine
                    Pu = zeros([2*m+1 d]);
                    % include existing nodes
                    Pu(1:2:end,1:2:end,:) = yn;
                    % prolongate in x direction
                    Pu(2:2:end-1,:,:) =  .5*(Pu(1:2:end-2,:,:) + Pu(3:2:end,:,:));
                    % prolongate in y direction
                    Pu(1:2:end,2:2:end-1,:) =  .5*(Pu(1:2:end,1:2:end-2,:) + Pu(1:2:end,3:2:end,:));
                    % prolongate diagonal
                    Pu(2:2:end,2:2:end,:) = .5* (Pu(1:2:end-2,1:2:end-2,:) + Pu(3:2:end,3:2:end,:));
                    
                    Pu = reshape(Pu,[],1);
                case 'PTu'
                    % include parent nodes
                    Pu = yn(1:2:end,1:2:end,:);
                    
                    % distribute in x direction
                    t = .5* yn(2:2:end-1,1:2:end,:);
                    Pu(1:end-1,:,:) = Pu(1:end-1,:,:) + t;
                    Pu(2:end,:,:)   = Pu(2:end,:,:) + t;
                    
                    % distribute in y direction
                    t = .5* yn(1:2:end,2:2:end-1,:);
                    Pu(:,1:end-1,:) = Pu(:,1:end-1,:) + t;
                    Pu(:,2:end,:) = Pu(:,2:end,:) + t;
                    
                    % distribute diagonal
                    t = .5* yn(2:2:end-1,2:2:end-1,:);
                    Pu(1:end-1,1:end-1,:) = Pu(1:end-1,1:end-1,:) + t;
                    Pu(2:end,2:end,:) = Pu(2:end,2:end,:) + t;
                    
                    Pu = reshape(Pu,[],1);
            end
        end
        
        function P = getPuNodal(~,m)
            % =============================================================
            % function P = getPuNodal(~,m)
            %
            % returns prolongation operator for input of cell-width m
            % =============================================================
            
            mf = 2*m;
            % indices of fine and coarse grid nodes
            indf = reshape(1:prod(mf+1),mf+1);
            indc = reshape(1:prod(m+1),m+1);
            
            % allocate space
            I = []; J = []; W = [];
            
            % include existing nodes
            ii = indf(1:2:end,1:2:end);
            jj = indc(:);
            ww = ones(size(jj));
            I = [I; ii(:)]; J = [J; jj(:)]; W = [W; ww(:)];
            
            % prolongate in x direction
            ii   = indf(2:2:end-1,1:2:end);
            jj   = [reshape(indc(1:end-1,:),[],1); reshape(indc(2:end,:),[],1)];
            ww = .5*ones(size(jj));
            I = [I; ii(:);ii(:)]; J = [J; jj(:)]; W = [W; ww(:)];
            
            % prolongate in y direction
            ii   = indf(1:2:end,2:2:end-1);
            jj   = [reshape(indc(:,1:end-1),[],1); reshape(indc(:,2:end),[],1)];
            ww   = .5*ones(size(jj));
            I = [I; ii(:); ii(:)]; J = [J; jj(:)]; W = [W; ww(:)];
            
            % prolongate diagonally
            ii   = indf(2:2:end-1,2:2:end-1);
            jj   = [reshape(indc(1:end-1,1:end-1),[],1); reshape(indc(2:end,2:end),[],1)];
            ww   = .5*ones(size(jj));
            I = [I; ii(:);ii(:)]; J = [J; jj(:)]; W = [W; ww(:)];
            
            P = sparse(I,J,W,prod(mf+1),prod(m+1));
        end
        
        
        
        % ========== get methods ========================================
        function dx1 = get.dx1(this)
            if isempty(this.dx1_),
                [this.dx1_, this.dx2_] = getGradientMatrixFEM(this,0);
            end
            dx1 = this.dx1_;
        end
        
        function dx2 = get.dx2(this)
            if isempty(this.dx2_),
                [this.dx1_, this.dx2_] = getGradientMatrixFEM(this,0);
            end
            dx2 = this.dx2_;
        end
        
        function GRAD = get.GRAD(this)
            if isempty(this.GRAD_),
                this.GRAD_ = [this.dx1;this.dx2];
            end
            GRAD = this.GRAD_;
        end
        function B = get.B(this)
            if isempty(this.B_),
                this.B_ = blkdiag(this.GRAD,this.GRAD);
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
        function PC = get.PC(this)
            if isempty(this.PC_),
                this.PC_ = (this.P1+this.P2+this.P3)/3;
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
                [this.mfdx1_, this.mfdx2_] = getGradientMatrixFEM(this,1);
            end
            mfdx1 = this.mfdx1_;
        end
        
        function mfdx2 = get.mfdx2(this)
            if isempty(this.mfdx2_),
                 [this.mfdx1_, this.mfdx2_] = getGradientMatrixFEM(this,1);
            end
            mfdx2 = this.mfdx2_;
        end
        
        function mfGRAD = get.mfGRAD(this)
            if isempty(this.mfGRAD_),
                this.mfGRAD_ = getGradientMatrixFEM(this,1);
            end
            mfGRAD = this.mfGRAD_;
        end
        function idx = get.boundaryIdx(this)
            if isempty(this.boundaryIdx_),
                id = reshape(1:prod(this.m+1),this.m+1);
                idx = [reshape(id([1,end],:),[],1);reshape(id(2:end-1,[1,end]),[],1)] ;
                this.boundaryIdx_ = idx;
            end
            idx = this.boundaryIdx_;
        end
        function idx = get.boundaryProj(this)
            if isempty(this.boundaryProj_),
                idx = this.boundaryIdx;
                
                P = speye(prod(this.m+1));
                P = P(idx,:);
                P = kron(speye(2),P);
                
                this.boundaryProj_ = P;
            end
            idx = this.boundaryProj_;
        end
      
       
        % ========== set methods ========================================
        function set.xn(this,xn)
            if numel(xn)~=(2*prod(this.m+1)),
                error('Invalid number of nodes');
            end
            if isempty(this.xn),
                this.xn = xn;
            else
                
                this.xn = xn;
                % delete operators that are sensitive to nodes
                this.dx1_  = [];
                this.dx2_  = [];
                this.GRAD_ = [];
                this.B_    = [];
                
                this.vol = volTetraGrid(this,this.xn,'matrixFree',1);
            end
        end
    end
end

