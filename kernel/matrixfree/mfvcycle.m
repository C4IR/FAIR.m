%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% function [u,res,r] = mfvcycle(para,u,rhs,tol,level,out)
%
% !! Multigrid, handle with care
%==============================================================================

function [u,res,r] = mfvcycle(para,u,rhs,tol,level,out)

if nargin == 0,
  help(mfilename);
  MGsolver;
  Pu = 'endOfMinimalExample';
  return
end;;

omega   = para.omega; 
m       = para.m; 
h       = (omega(2:2:end)-omega(1:2:end))./m;
MGcycle = para.MGcycle;
para.MGcycle   = 1;
dimstr  = @(m) sprintf('[%s]',sprintf(' %d',m));

if ~isfield(para,'M'), para.M = 0;       end;

if level == para.MGlevel & out==1,
  fprintf('MG: ');
end;
if (level>1) & (min(m)==1)
  %if out, fprintf('grid %s --> coarse grid',dimstr(m)); end;
  level = 1;
end;
%out = 2;
if out==1,
  fprintf('%d>',level);
end;

% -----------------------------------------------------------------------------
if level > 1,
  % ---------------------------------------------------------------------------
  % NOT on COARSE grid
  str = sprintf('  - grid %s',dimstr(m));
  
  % initialize
  r = rhs - mfAy(u,para);

  % -------------------------------------------------------------
  for i = 1:MGcycle, % loop over the V/W cycles
  % -------------------------------------------------------------
    n = length(r);
    if out>1,
      norm_r_in = norm(r)/n;
      fprintf('%20s |r-in|   = %e\n',str,norm_r_in);
    end;

    % pre-smoothing, exp: [u,r] Richardson(x,b,para,MGcycle)
    ur = feval(para.MGsmoother,0*r,r,para,para.MGpresmooth);
    u = u + ur; r = rhs - mfAy(u,para);

    if out>1,
      norm_r_pre = norm(r)/n;
      fprintf('%20s |r-pre|  = %e\n','',norm_r_pre);
    end;
    
    % PREPARE FOR COARSER GRID ------------------------------------------------
    M = para.d2D.M;
    if ~(isempty(m) || all(size(M)==[1,1]) && M==0 )
      para.d2D.M = mfPu(M,para.omega,para.m,'PTu')/4;
    end;
    rc = mfPu(r,para.omega,para.m,'PTu')/4;
    para.m = m/2;

    % RECURSIVE CALL -------------------------------------------------
    % Solve the coarse grid system
    uc = mfvcycle(para,0*rc,rc,1e-16,level-1,out);

    % prolongate, back to fine grid
    ur = mfPu(uc,para.omega,para.m,'Pu');
    
    
    para.m = m;
    para.h = h;        
    para.d2D.M = M;
    
    % ----------------------------------------------------------------
    
    % update
    u = u + ur;  r = rhs - mfAy(u,para);

    if out>0,
      norm_r_app = norm(r)/n;
    end;
    
    % post-smoothing 
    ur = feval(para.MGsmoother,0*r,r,para,para.MGpostsmooth);
    u = u + ur;  r = rhs - mfAy(u,para);
    
    if out==1,
      fprintf('<%d',level);
    elseif out>1,
      norm_r_post = norm(r)/n;
      
      fprintf('%20s |r-in|   = %e\n',str,norm_r_in);
      fprintf('%20s |r-pre|  = %e\n','',norm_r_pre);
      fprintf('%20s |r-app|  = %e\n','',norm_r_app);
      fprintf('%20s |r-post| = %e\n','',norm_r_post);
    end;
    
    if para.MGcycle > 2,
      fprintf('iter = %d,  |res| = %e\n',i,norm(r))
    end;
    res(i) = norm(r);
    if res(i) < tol, 
      resi = res(i);
      return; 
    end;
    
    if level == para.MGlevel && out  ~= 0,
      fprintf('>');
    end;
  end;

  % FINE GRID DONE
  % ---------------------------------------------------------------------------
else
  % ---------------------------------------------------------------------------  
  % COARSE grid
  str = sprintf('  - coarse grid %s',dimstr(m));
  norm_r_in = norm(rhs);

  if length(rhs)>100, error('rhs too big');  end;
  
  B = para.d2S.B(omega,m);

  A = spdiags(para.d2D.M,0,length(rhs),length(rhs)) + para.d2S.alpha*prod(h)*B'*B;  
  u = pinv(full(A))*rhs;
  r = rhs - A*u;

  norm_r_out = norm(r);
  res = norm_r_out;
  if out==1,
    fprintf('coarse(%s)',dimstr(m));
  elseif out>1,
    fprintf('%20s |r-in|   = %e\n',str,norm_r_in);
    fprintf('%20s |r-out|  = %e\n',str,norm_r_out);
  end;
end;
% -----------------------------------------------------------------------------

if level == para.MGlevel & out>0,
  if out==1,
    fprintf('.\n');
  end;
 
  %testMG = norm(rhs - mfAy(u,para));
  %fprintf('MG:|rhs-Au(u)|=%12.4e\n',testMG)
end;

%==============================================================================
