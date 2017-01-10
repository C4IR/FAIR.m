%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% function [rho,drho] = rhoSpline(Tc,Rc,minT,maxT,nT,minR,maxR,nR,doDerivative)
%
% Parzen-Window Based density estimator for gray values of Tc and Rc, spline based.
%
% Input: 
%  Tc, Rc          template and reference
%  minT, maxT, nT  discreted points in gray value range where splines are located
%  minR, maxR, nR  discreted points in gray value range where splines are located
%  doSerivative         FLAG for derivative compuatation
%
% Output:
%  rho          joint density estimator
%  drho         derivative of joint density estimator
% 
% see also MIspline for an example.
%==============================================================================

function [rho,drho] = rhoSpline(Tc,Rc,minT,maxT,nT,minR,maxR,nR,doDerivative)

if nargin == 0,
  help(mfilename);
  return;
end;

if ~exist('doDerivative','var'), doDerivative = 1;       end;
doDerivative = doDerivative & (nargout>1);

% prepare the output and organize input
rho  = []; drho = []; Tc = reshape(Tc,[],1); Rc = reshape(Rc,[],1);

widthT = (maxT-minT)/nT;
widthR = (maxR-minR)/nR;

Tt   = linspace(minT,maxT,nT);
Rt   = linspace(minR,maxR,nR);
rho  = zeros(length(Tt),length(Rt));
drho = spalloc(numel(rho),length(Tc),5*length(Tc));

for j=1:length(Tc), % run over all samples
  IT = find(abs(Tt-Tc(j))<widthT); % find locations of interest in Tc
  IR = find(abs(Rt-Rc(j))<widthR); % find locations of interest in Rc
  if ~isempty(IT) && ~isempty(IR),
    [KT,dKT] = spline1D(Tc(j) - Tt(IT),widthT,doDerivative);
    KR       = spline1D(Rc(j) - Rt(IR),widthR,0);
    rho(IT,IR) = rho(IT,IR) +KT*KR'; 
    if doDerivative,
      drhoj = zeros(size(rho));
      drhoj(IT,IR) = dKT*KR';
      drho(:,j) = drho(:,j) + sparse(drhoj(:));
    end;
  end;
end;

% normalize rho
fac = (Tt(2)-Tt(1))*(Rt(2)-Rt(1))/length(Tc);
rho  = fac*rho(:);
drho = fac*drho;

%------------------------------------------------------------------------------

% The Parzen Window Kernel is a spline function
%  function [s,ds] = spline1D(x,sigma,doDerivative);
% (c) Jan Modersitzki 2008/02/12, see FAIR.
% This evaluates a spline function s
% s(x,sigma) = 0 for x\notin[-sigma,sigma], \int s(x,sigma) dx = 1

function [s,ds] = spline1D(x,sigma,doDerivative);

if ~exist('doDerivative','var'), doDerivative = 1;       end;
doDerivative = doDerivative & (nargout>1);

% map [-sigm2,sigma] \to [-2,2]
x = 4/sigma*x;

if ~doDerivative,
  s  = 0*reshape(x,[],1);
  ds = [];
  J = find(x>=-2 & x<-1);   s(J)   = (x(J)+2).^3;
  J = find(x>=-1 & x< 0);   s(J)   = -x(J).^3 - 2*(x(J)+1).^3+6*(x(J)+1);
  J = find(x>= 0 & x< 1);   s(J)   = -2*(1-x(J)).^3 + x(J).^3 - 6*x(J) + 6;
  J = find(x>= 0 & x< 1);   s(J)   = -2*(1-x(J)).^3 + x(J).^3 - 6*x(J) + 6;
  J = find(x>= 1 & x< 2);   s(J)   = (2-x(J)).^3;
else
  s  = 0*reshape(x,[],1);
  ds = s;
  J = find(x>=-2 & x<-1);   s(J)   = (x(J)+2).^3;
                            ds(J)  = 3*(x(J)+2).^2;        
  J = find(x>=-1 & x< 0);   s(J)   = -x(J).^3 - 2*(x(J)+1).^3+6*(x(J)+1);
                            ds(J)  = -3*x(J).^2 - 6*(x(J)+1).^2+6;  
  J = find(x>= 0 & x< 1);   s(J)   = -2*(1-x(J)).^3 + x(J).^3 - 6*x(J) + 6;
                            ds(J)  = 6*(1-x(J)).^2 + 3*x(J).^2 - 6; 
  J = find(x>= 1 & x< 2);   s(J)   = (2-x(J)).^3;
                            ds(J)  = -3*(2-x(J)).^2;
end;
% scale to make integral 1
s  = 2/(3*sigma)*s;
ds = 8/(3*sigma^2)*ds;
%==============================================================================

