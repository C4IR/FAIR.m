function [Dc,rho,dD,drho,d2psi] = MI(Tc,Rc,omega,m,varargin);
% This is a link to MIcc, the most efficient implementation of Mutual Information

if nargin == 0,  
  help('MImex');  
  MImex;  
  return; 
end;

[Dc,rho,dD,drho,d2psi] = MImex(Tc,Rc,omega,m,varargin{:});
