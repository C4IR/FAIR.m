%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% function [G,dG] = geometry(yc,m,flag,varargin)
%
% area, volume and det computation for tetrahedral partition of yc
%
% Input:
% ------
%    yc         nodal grid
%    m             number of discretization points
%    flag        can be {'A','V','Jac','Vrange'}
%    varargin   optional parameter (see below)
%
% Output:
% ------
%    G           value
%    dG          derivative 
%   if ~matrixFree, dG is sparse matrix; else, dG is struct endif        
%
% See also 
%  @article{2013-BMR,
%        Author    = {Burger M., Modersitzki J., Ruthotto L.},
%        Title     = {A hyperelastic regularization energy for image registration},
%     Journal   = {SIAM SISC},
%     Volume    = {35},
%        Year      = {2013},
%  }
% 
% =============================================================================

function [G,dG] = geometry(yc,m,flag,varargin)


if nargin==0,
    runMinimalExample;
    G = 'endOfMinimalExample';
    return;
end

matrixFree   = regularizer('get','matrixFree');
doDerivative = (nargout>1);
dG           = [];

for k=1:2:length(varargin) % overwrites defaults
    eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end

dim = length(m);        % image dimension
% n   = prod(m);          % number gid-points
% nn  = prod(m+1);        % number of nodal points
if exist('omega','var') && not(isempty(omega)),
    hd = prod((omega(2:2:end)-omega(1:2:end)) ./m);
end


flags = {'-mb-','-mf-'};
flag  = [flag, flags{matrixFree + 1},num2str(dim)];


switch flag
    %{
    
    Volumes for 2D and 3D
    
    Compute volume for 2D
    We assume a yc on a nodal grid like



            PA --- PB ----- * ----- * ----- *
            |       |       |       |       |
            |   PM  |   o   |   o   |   o   |
            |       |       |       |       |
            PC --- PD ----- * ----- * ----- *
            |       |       |       |       |
            |   o   |   o   |   o   |   o   |
            |       |       |       |       |
            * ----- * ----- * ----- * ----- *
            |       |       |       |       |
            |   o   |   o   |   o   |   o   |
            |       |       |       |       |
            * ----- * ----- * ----- * ----- *

    Let S = speye(prod(m+1),prod(m+1)) by an indicator matrix for indices
    The indices are indicated by
    PA : top    left  locations, S(bottom,:) = [] and S(right,:) = []
    PB : top    right locations, S(bottom,:) = [] and S(left,:)  = []
    PC : bottom left  locations, S(top,:)    = [] and S(right,:) = []
    PD : bottom right locations, S(top,:)    = [] and S(left,:)  = []
    PM : average to cell centers : cc(:,1:2) = PM*yc(:,1:2)

        Compute the area of the four triangles 

            (PA,PB,PM),  (PB,PD,PM), (PD,PC,PM), and (PC,PA,PM)
             
            PA --- PB         PB                       PA   
              \   /         / |                        |  \
               PM         PM  |           PM           |   PM
                            \ |         /   \          |  /
                              PD       PC --- PD        PC
    
    
    
    Same idea for 3D
                 PE ----- PF
                /|      / |
               / |     /  |
              /  PG --/-- PH
            PA ----- PB  /
            |  /     |  /
            | /      | /
            PC ----- PD      
    
    
    
    %}
    
    case {'V-mb-2','Jac-mb-2'} % matrixFree = 0
        P = getNodeProjection(m);
        
        [volABM dvolABM] = areaOfTriangle2D(P.A,P.B,P.M,yc,doDerivative);
        [volBDM dvolBDM] = areaOfTriangle2D(P.B,P.D,P.M,yc,doDerivative);
        [volDCM dvolDCM] = areaOfTriangle2D(P.D,P.C,P.M,yc,doDerivative);
        [volCAM dvolCAM] = areaOfTriangle2D(P.C,P.A,P.M,yc,doDerivative);
        
        
        if strcmp(flag,'V-mb-2')    
            % volume
            G = [volABM ; volBDM; volDCM; volCAM];
            if doDerivative,
                dG =[dvolABM ; dvolBDM; dvolDCM; dvolCAM];
            end
        else
            % Jacobian
            G = (volABM + volBDM + volDCM + volCAM) /hd;
            if doDerivative,
                dG =(dvolABM + dvolBDM + dvolDCM + dvolCAM) / hd;
            end
        end;
        

    case {'V-mb-3','Jac-mb-3'} % matrixFree = 0
        P = getNodeProjection(m);
        
        % compute volumes of 24 tetras per voxel
        [V1,  dV1] = volTetra3D(yc,P.A,P.B,P.ABDC,P.M,doDerivative);
        [V2,  dV2] = volTetra3D(yc,P.B,P.D,P.ABDC,P.M,doDerivative);
        [V3,  dV3] = volTetra3D(yc,P.D,P.C,P.ABDC,P.M,doDerivative);
        [V4,  dV4] = volTetra3D(yc,P.C,P.A,P.ABDC,P.M,doDerivative);
        
        [V5,  dV5] = volTetra3D(yc,P.B,P.D,P.BDHF,P.M,doDerivative);
        [V6,  dV6] = volTetra3D(yc,P.D,P.H,P.BDHF,P.M,doDerivative);
        [V7,  dV7] = volTetra3D(yc,P.H,P.F,P.BDHF,P.M,doDerivative);
        [V8,  dV8] = volTetra3D(yc,P.F,P.B,P.BDHF,P.M,doDerivative);
        
        [V9,  dV9] = volTetra3D(yc,P.A,P.B,P.ABFE,P.M,doDerivative);
        [V10,dV10] = volTetra3D(yc,P.B,P.F,P.ABFE,P.M,doDerivative);
        [V11,dV11] = volTetra3D(yc,P.F,P.E,P.ABFE,P.M,doDerivative);
        [V12,dV12] = volTetra3D(yc,P.E,P.A,P.ABFE,P.M,doDerivative);
        
        [V13,dV13] = volTetra3D(yc,P.C,P.A,P.CAEG,P.M,doDerivative);
        [V14,dV14] = volTetra3D(yc,P.A,P.E,P.CAEG,P.M,doDerivative);
        [V15,dV15] = volTetra3D(yc,P.E,P.G,P.CAEG,P.M,doDerivative);
        [V16,dV16] = volTetra3D(yc,P.G,P.C,P.CAEG,P.M,doDerivative);
        
        [V17,dV17] = volTetra3D(yc,P.D,P.C,P.DCGH,P.M,doDerivative);
        [V18,dV18] = volTetra3D(yc,P.C,P.G,P.DCGH,P.M,doDerivative);
        [V19,dV19] = volTetra3D(yc,P.G,P.H,P.DCGH,P.M,doDerivative);
        [V20,dV20] = volTetra3D(yc,P.H,P.D,P.DCGH,P.M,doDerivative);
        
        [V21,dV21] = volTetra3D(yc,P.E,P.F,P.EFHG,P.M,doDerivative);
        [V22,dV22] = volTetra3D(yc,P.F,P.H,P.EFHG,P.M,doDerivative);
        [V23 dV23] = volTetra3D(yc,P.H,P.G,P.EFHG,P.M,doDerivative);
        [V24,dV24] = volTetra3D(yc,P.G,P.E,P.EFHG,P.M,doDerivative);
        
        if strcmp(flag,'V-mb-3')
            % volume
            G = [-V1;-V2;-V3;-V4;V5;V6;V7;V8;V9;V10;V11;V12;...
                V13;V14;V15;V16;V17;V18;V19;V20;V21;V22;V23;V24];
            if doDerivative,
                dG = [-dV1;-dV2;-dV3;-dV4;dV5;dV6;dV7;dV8;dV9;dV10;dV11;dV12;...
                    dV13;dV14;dV15;dV16;dV17;dV18;dV19;dV20;dV21;dV22;dV23;dV24];
            end
        else
            % Jacobian
            G = (-V1-V2-V3-V4+V5+V6+V7+V8...
                +V9+V10+V11+V12+V13+V14+V15+V16...
                +V17+V18+V19+V20+V21+V22+V23+V24) /hd;
            if doDerivative,
                dG = (-dV1-dV2-dV3-dV4+dV5+dV6+dV7+dV8...
                    +dV9+dV10+dV11+dV12+dV13+dV14+dV15+dV16...
                    +dV17+dV18+dV19+dV20+dV21+dV22+dV23+dV24) /hd;
            end
        end;
        
        
    case {'VRange-mb-2','VRange-mb-3'}      % matrixFree = 0
        vol = geometry(yc,m,'V','matrixFree',matrixFree);
        G   = [ min(vol(:)); max(vol(:)) ];

    case 'A-mb-3'
        P = getNodeProjection(m);
        [A1,  dA1] = areaOfTriangle3D(yc,P.A,P.B,P.ABDC,doDerivative);
        [A2,  dA2] = areaOfTriangle3D(yc,P.B,P.D,P.ABDC,doDerivative);
        [A3,  dA3] = areaOfTriangle3D(yc,P.D,P.C,P.ABDC,doDerivative);
        [A4,  dA4] = areaOfTriangle3D(yc,P.C,P.A,P.ABDC,doDerivative);
        
        [A5,  dA5] = areaOfTriangle3D(yc,P.B,P.D,P.BDHF,doDerivative);
        [A6,  dA6] = areaOfTriangle3D(yc,P.D,P.H,P.BDHF,doDerivative);
        [A7,  dA7] = areaOfTriangle3D(yc,P.H,P.F,P.BDHF,doDerivative);
        [A8,  dA8] = areaOfTriangle3D(yc,P.F,P.B,P.BDHF,doDerivative);
        
        [A9,  dA9]  = areaOfTriangle3D(yc,P.A,P.B,P.ABFE,doDerivative);
        [A10, dA10] = areaOfTriangle3D(yc,P.B,P.F,P.ABFE,doDerivative);
        [A11, dA11] = areaOfTriangle3D(yc,P.F,P.E,P.ABFE,doDerivative);
        [A12, dA12] = areaOfTriangle3D(yc,P.E,P.A,P.ABFE,doDerivative);
        
        [A13, dA13] = areaOfTriangle3D(yc,P.C,P.A,P.CAEG,doDerivative);
        [A14, dA14] = areaOfTriangle3D(yc,P.A,P.E,P.CAEG,doDerivative);
        [A15, dA15] = areaOfTriangle3D(yc,P.E,P.G,P.CAEG,doDerivative);
        [A16, dA16] = areaOfTriangle3D(yc,P.G,P.C,P.CAEG,doDerivative);
        
        [A17, dA17] = areaOfTriangle3D(yc,P.D,P.C,P.DCGH,doDerivative);
        [A18, dA18] = areaOfTriangle3D(yc,P.C,P.G,P.DCGH,doDerivative);
        [A19, dA19] = areaOfTriangle3D(yc,P.G,P.H,P.DCGH,doDerivative);
        [A20, dA20] = areaOfTriangle3D(yc,P.H,P.D,P.DCGH,doDerivative);
        
        [A21, dA21] = areaOfTriangle3D(yc,P.E,P.F,P.EFHG,doDerivative);
        [A22, dA22] = areaOfTriangle3D(yc,P.F,P.H,P.EFHG,doDerivative);
        [A23, dA23] = areaOfTriangle3D(yc,P.H,P.G,P.EFHG,doDerivative);
        [A24, dA24] = areaOfTriangle3D(yc,P.G,P.E,P.EFHG,doDerivative);                

        G = [A1;A2;A3;A4;A5;A6;A7;A8;A9;A10;A11;A12;....
            A13;A14;A15;A16;A17;A18;A19;A20;A21;A22;A23;A24];
        if doDerivative,
            dG = [dA1;dA2;dA3;dA4;dA5;dA6;dA7;dA8;dA9;dA10;dA11;dA12;...
                dA13;dA14;dA15;dA16;dA17;dA18;dA19;dA20;dA21;dA22;dA23;dA24];
        end
        

    case {'V-mf-2','V-mf-3'},               % matrixFree = 1
        try
            G = geometryMexC(yc(:),m,'V');
        catch err
            throw(err);
        end


    case {'VRange-mf-2','VRange-mf-3'}      % matrixFree = 1
        try
            G = geometryMexC(yc(:),m,'VRange');
        catch err
            throw(err);
        end


    case {'Jac-mf-2','Jac-mf-3'}            % matrixFree = 1
        try
            G = geometryMexC(yc(:),m,'Jac')/hd;
        catch err
            throw(err);
        end
        if doDerivative,
            dG.dJac    = @(y,m,x) geometryMexC(y(:),m,'dJacx',   x(:))/hd;
            dG.dJacadj = @(y,m,x) geometryMexC(y(:),m,'dJacadjx',x(:))/hd;
        end
        
        
    case 'A-mf-3'
        try
            G = geometryMexC(yc(:),m,'A');
        catch err
            throw(err);
        end

        %{
    case {'dJacx-mf-2','dJacx-mf-3'}
        keyboard
        try
            G = geometryMexC(yc(:),m,'dJacx',x(:));
        catch err
            throw(err);
        end
        
    case {'dJacadjx-mf-2','dJacadjx-mf-3'}
        keyboard
        try
            G = geometryMexC(yc(:),m,'dJacadjx',x(:));
        catch err
            throw(err);
        end
    case {'V-mf-2','V-mf-3'}
        keyboard
        try
            G = geometryMexC(yc(:),m,'V');
        catch err
            throw(err);
        end
    %}
        
        
    otherwise
        error('Error in file geometry.m : flag '' %s '' nyi!',flag);

end
%------------------------------------------------------------------------------
% end geometry.m
%------------------------------------------------------------------------------



%------------------------------------------------------------------------------
% supprt functions for the matric based case
%------------------------------------------------------------------------------


%------------------------------------------------------------------------------
function Q = getNodeProjection(m)
%------------------------------------------------------------------------------
%{
    setup of projection matrices to the node
            A ----- B ----- * ----- * ----- *
            |       |       |       |       |
            |   M   |   o   |   o   |   o   |
            |       |       |       |       |
            C ----- D ----- * ----- * ----- *
            |       |       |       |       |
            |   o   |   o   |   o   |   o   |
            |       |       |       |       |
            * ----- * ----- * ----- * ----- *
            |       |       |       |       |
            |   o   |   o   |   o   |   o   |
            |       |       |       |       |
            * ----- * ----- * ----- * ----- *

    Let S = speye(prod(m+1),prod(m+1)) by an indicator matrix for indices
    The indices are indicated by
    A : top    left  locations, S(bottom,:) = [] and S(right,:) = []
    B : top    right locations, S(bottom,:) = [] and S(left,:)  = []
    C : bottom left  locations, S(top,:)    = [] and S(right,:) = []
    D : bottom right locations, S(top,:)    = [] and S(left,:)  = []
    M : average to cell centers : cc(:,1:2) = PM*yc(:,1:2)

    
    Same idea for 3D
                 E ----- F
                /|      /|
               / |     / |
              /  G ---/--H
            A ------ B  /
            |  /     | /
            | /      |/
            C ------ D      

%}

persistent S P
if isempty(S) || size(S,1)~=prod(m+1),
    S   = speye(prod(m+1),prod(m+1));
    
    if length(m) == 2,
        P.A = S(1:end-(m(1)+1),:);  P.A((m(1)+1)*(1:m(2)),:) = [];
        P.B = S(1:end-(m(1)+1),:);  P.B((m(1)+1)*(0:m(2)-1)+1,:) = [];
        P.C = S(m(1)+2:end,:);      P.C((m(1)+1)*(1:m(2)),:) = [];
        P.D = S(m(1)+2:end,:);      P.D((m(1)+1)*(0:m(2)-1)+1,:) = [];
        P.M = (1/4) * (P.A + P.B + P.C + P.D);
    elseif length(m) == 3,
        
        left  = @(k) ((1:m(k)+1) < (m(k)+1));
        right = @(k) ((1:m(k)+1) > 1    );
        
        P.A = S(kron(kron(left(3),left(2)),left(1)) == 1,:);
        P.B = S(kron(kron(left(3),left(2)),right(1)) == 1,:);
        P.C = S(kron(kron(left(3),right(2)),left(1)) == 1,:);
        P.D = S(kron(kron(left(3),right(2)),right(1)) == 1,:);
        P.E = S(kron(kron(right(3),left(2)),left(1)) == 1,:);
        P.F = S(kron(kron(right(3),left(2)),right(1)) == 1,:);
        P.G = S(kron(kron(right(3),right(2)),left(1)) == 1,:);
        P.H = S(kron(kron(right(3),right(2)),right(1)) == 1,:);
        
        P.ABDC = ( P.A + P.B + P.C + P.D) /4;
        P.BDHF = ( P.B + P.D + P.H + P.F) /4;
        P.ABFE = ( P.A + P.B + P.F + P.E) /4;
        P.CAEG = ( P.C + P.A + P.E + P.G) /4;
        P.DCGH = ( P.D + P.C + P.G + P.H) /4;
        P.EFHG = ( P.E + P.F + P.H + P.G) /4;
        P.M    = ( P.A + P.B + P.C + P.D + P.E + P.F + P.G + P.H ) /8;
    else
        error('dimesion nyi')
    end;
end;
Q = P;
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
function [V,dV] =  areaOfTriangle2D(P1,P2,PM,yc,doDerivative)
%------------------------------------------------------------------------------
%{
  Given a list of points yc, and filters P1 P2, PM for triangle points 

  Compute the area of the triangles by cross-product or determinant
  area = 0.5 * (P1-PM)x(P2-PM) = 0.5 * det([P1,PM,P2-PM])
%}

  dV = [];
  yc = reshape(yc,[],2);
  V  = ((P1*yc(:,1)-PM*yc(:,1)) .* (P2*yc(:,2)-PM*yc(:,2)) ...
      - (P2*yc(:,1)-PM*yc(:,1)) .* (P1*yc(:,2)-PM*yc(:,2)))/2;

  if doDerivative,
    dV =[
        sdiag(P2*yc(:,2)-PM*yc(:,2))*(P1-PM) ...
        - sdiag(P1*yc(:,2)-PM*yc(:,2))*(P2-PM),...
        sdiag(P1*yc(:,1)-PM*yc(:,1))*(P2-PM) ...
        - sdiag(P2*yc(:,1)-PM*yc(:,1))*(P1-PM)
        ]/2;
end
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
function [vol dVol] = volTetra3D(yc,P1,P2,PF,PM,doDerivative)
%------------------------------------------------------------------------------
% computes volume of tetrahedra with edges P1 and P2, face-staggered point
% PF and cell-center PM
%
% volume is given by the determinant of the matrix
%
%     |       |       |       |
% M = | P1-PM | P2-PM | PF-PM |
%     |       |       |       |
%       col1  col2  col3
%
dVol = [];
yc = reshape(yc,[],3);
P1 = (P1-PM);
P2 = (P2-PM);
P3 = (PF-PM);
col1 = P1*yc;
col2 = P2*yc;
col3 = P3*yc;

vol  = (1/6)*(col1(:,1) .* col2(:,2) .* col3(:,3) ...
    + col2(:,1) .* col3(:,2) .* col1(:,3) ...
    + col3(:,1) .* col1(:,2) .* col2(:,3) ...
    - col1(:,3) .* col2(:,2) .* col3(:,1) ...
    - col2(:,3) .* col3(:,2) .* col1(:,1) ...
    - col3(:,3) .* col1(:,2) .* col2(:,1));
if ~doDerivative,
    return;
end

dVol = (1/6)*[...
    sdiag(col2(:,2) .* col3(:,3) - col2(:,3) .* col3(:,2)) * P1 ...
    + sdiag(col3(:,2) .* col1(:,3) - col3(:,3) .* col1(:,2)) * P2 ...
    + sdiag(col1(:,2) .* col2(:,3) - col1(:,3) .* col2(:,2)) * P3,...
    sdiag(col1(:,1) .* col3(:,3) - col1(:,3) .* col3(:,1)) * P2...
    + sdiag(col2(:,1) .* col1(:,3) - col2(:,3) .* col1(:,1)) * P3...
    + sdiag(col3(:,1) .* col2(:,3) - col3(:,3) .* col2(:,1)) * P1,...
    + sdiag(col1(:,1) .* col2(:,2) - col1(:,2) .* col2(:,1)) * P3...
    + sdiag(col2(:,1) .* col3(:,2) - col2(:,2) .* col3(:,1)) * P1...
    + sdiag(col3(:,1) .* col1(:,2) - col3(:,2) .* col1(:,1)) * P2...
    ];
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
function [A dA] = areaOfTriangle3D(yc,P1,P2,PF,doDerivative)
%------------------------------------------------------------------------------
%{
    computes (squared) volume of triangle with edges P1 and P2 and the
    face-staggered point PF

    area is given by the euclidean length of the cross product

    ( B-A ) x ( C-A )
    =: v      =: w
%}

dA = [];
yc = reshape(yc,[],3);
Pv = (P1-PF);
Pw = (P2-PF);
v  = Pv*yc;
w  = Pw*yc;

vec = [...
    (v(:,2) .* w(:,3) - v(:,3) .* w(:,2)), ...
    (v(:,1) .* w(:,3) - v(:,3) .* w(:,1)), ...
    (v(:,1) .* w(:,2) - v(:,2) .* w(:,1)) ];

A   = sum(vec.^2,2);

if ~doDerivative,
    return;
end

dA = 2* [
    sdiag( (v(:,1) .* w(:,3) - v(:,3) .* w(:,1)).* w(:,3) ...
    +(v(:,1) .* w(:,2) - v(:,2) .* w(:,1)).* w(:,2)) * Pv ...
    - sdiag((v(:,1) .* w(:,3) - v(:,3) .* w(:,1)).* v(:,3) ...
    +(v(:,1) .* w(:,2) - v(:,2) .* w(:,1)).* v(:,2))* Pw,...
    sdiag( (v(:,2) .* w(:,3) - v(:,3) .* w(:,2)).* w(:,3) ...
    -(v(:,1) .* w(:,2) - v(:,2) .* w(:,1)).* w(:,1)) * Pv ...
    - sdiag((v(:,2) .* w(:,3) - v(:,3) .* w(:,2)).* v(:,3)...
    -((v(:,1) .* w(:,2) - v(:,2) .* w(:,1)).* v(:,1)))* Pw,...
    - sdiag((v(:,2) .* w(:,3) - v(:,3) .* w(:,2)).* w(:,2) ...
    +(v(:,1) .* w(:,3) - v(:,3) .* w(:,1)).* w(:,1) )* Pv...
    + sdiag( (v(:,2) .* w(:,3) - v(:,3) .* w(:,2)).* v(:,2)...
    +(v(:,1) .* w(:,3) - v(:,3) .* w(:,1)).* v(:,1)) * Pw ...
    ];
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
function a = sdiag(a)
% shortcut for sparse diagonal matrices
a = spdiags(reshape(a,[],1),0,length(a),length(a));
%------------------------------------------------------------------------------

function runMinimalExample

help(mfilename)

%------------------------------------------------------------------------------
fprintf('Testing 2D \n');
omega = [0 2 1 3];
m     = [4 9];
tol   = 1e-14; %% test against condition*machine precision
xc    = getNodalGrid(omega,m);
xc    = xc + 1e-2 * randn(size(xc));


show = @(a,e) fprintf('%-30s, rel. error = %8.2e, passed = %1d\n',a,e,e<tol);

% Volume
[mbV, mbdV] = geometry(xc,m,'V','x',[],'matrixFree',0);
[mfV, mfdV] = geometry(xc,m,'V','x',[],'matrixFree',1);
RE = norm(mbV-mfV)/norm(mbV);
show('Test Volume (Option [V])',RE)

% Volume Range
[mbVRange] = geometry(xc,m,'VRange','x',[],'matrixFree',0);
[mfVRange] = geometry(xc,m,'VRange','x',[],'matrixFree',1);
RE = norm(mbVRange-mfVRange)/norm(mbVRange);
show('Test Range (Option [Vrange])',RE)

% Jacobian
xt = randn(size(xc));
[mbJac, mbdJac] = geometry(xc,m,'Jac','x',xt,'matrixFree',0,'omega',omega);
[mfJac, mfdJac] = geometry(xc,m,'Jac','x',xt,'matrixFree',1,'omega',omega);
RE =norm(mbJac-mfJac)/norm(mbJac);
show('Test Jacobian (Option [Jac])',RE)
  
% Derivative of Jacobian
xt = 1e2*randn(size(xc));
RE = norm(mbdJac*xt-mfdJac.dJac(xc,m,xt)) / norm(mbdJac*xt);
show('Test derivative Jacobian',RE)

% Adjoint of derivative of Jacobian
xt = 1e2*randn(prod(m),1);
RE = norm(mbdJac'*xt-mfdJac.dJacadj(xc,m,xt)) /norm(mbdJac'*xt);
show('Test adjoint of Jacobian',RE)
  


fprintf('\nTesting 3D \n');
omega = [0 2 1 3 4 6];
m     = [2 3 5];
xc    = getNodalGrid(omega,m);
xt    = xc + 1e-2 * randn(size(xc));

% Volume
[mbV, mbdV] = geometry(xc,m,'V','x',[],'matrixFree',0);
[mfV, mfdV] = geometry(xc,m,'V','x',[],'matrixFree',1);
RE = norm(mbV-mfV)/norm(mbV);
show('Test Volume (Option [V])',RE)

% Volume Range
[mbVRange] = geometry(xc,m,'VRange','x',[],'matrixFree',0);
[mfVRange] = geometry(xc,m,'VRange','x',[],'matrixFree',1);
RE = norm(mbVRange-mfVRange)/norm(mbVRange);
show('Test Range (Option [Vrange])',RE)

% Jacobian
[mbJac, mbdJac] = geometry(xc,m,'Jac','x',xt,'matrixFree',0,'omega',omega);
[mfJac, mfdJac] = geometry(xc,m,'Jac','x',xt,'matrixFree',1,'omega',omega);
RE = norm(mbJac-mfJac)/norm(mbJac);
show('Test Jacobian (Option [Jac])',RE)
  
% Derivative of Jacobian
xt = randn(size(xc));
RE = norm(mbdJac*xt-mfdJac.dJac(xc,m,xt)) / norm(mbdJac*xt);
show('Test derivative Jacobian',RE)

% Adjoint of derivative of Jacobian
xt = randn(prod(m),1);
RE = norm(mbdJac'*xt-mfdJac.dJacadj(xc,m,xt)) /norm(mbdJac'*xt);
show('Test adjoint of Jacobian',RE)

% Test area
mbA = geometry(xc,m,'A','matrixFree',0);
mfA = geometry(xc,m,'A','matrixFree',1);
RE  = norm(mbA-mfA) /norm(mbA);
show('Test area (option [A])',RE)
% ==================================================================================

