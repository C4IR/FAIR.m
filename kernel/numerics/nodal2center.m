%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
% function P = nodal2center(y,m,mex)
%
% transfers a nodal grid to a cell-centered grid
% if nargin==1, builds P explicitely, else return results of P(y), matrix free; endif
% depending on numel(Y), the matrix free version also handles P'(Y)
%
% Input:
%   y        input points,  nodal or cell-centered
%   m        number of discretization points
%   mex      bool (use mex code for matrix-free implementation)
%
% Output:
%   P        the projection matrix P, if nargin == 2
%            P*y,                     if y is nodal
%            P'*y,                    if y is cell-centered
%==============================================================================

function P = nodal2center(y,m,mex)

if nargin == 0 % help and minimal example
  help(mfilename);
  runMinimalExample;
  P = 'endMinimalExample';
  return
end

Js = {[1 2 3], [2 1 3], [3 1 2]}; % different permutations of dimensions

if nargin==1, m = y; end
dim = length(m);

% Here starts the matrix-based code
if nargin==1
  % ----------------------------------
  % build the matrix P and return it
  % ----------------------------------
  av = @(i) spdiags(ones(m(i),1)*[ 1,1],[0,1],m(i),m(i)+1)/2;
  zero = sparse(prod(m),prod(m+1));
  switch dim
    case 2
      A = kron(av(2), av(1));
      P = sparse([A zero; zero A]);
    case 3
      A = kron(av(3), kron(av(2), av(1)));
      P = sparse([A zero zero; zero A zero; zero zero A]);
    otherwise
      error('Dimension must be either 2 or 3.')
  end
  return; % end of story
end

% Here starts the MEX based code
if (dim > 1) && exist(fullfile(FAIRpath,'kernel',[mfilename,'C.',mexext])) == 3,
  if numel(y) == length(m)*prod(m)
    % cell-centered ->  nodal
    status = true;
  else
    % nodal -> cell-centered
    status = false;
  end;
  P = nodal2centerC(y, m, dim, status);
  return;
end;


% Here starts the matrix-free code

if numel(y) == length(m)*prod(m)
  % -----------------------
  % cell-centered ->  nodal
  % -----------------------
  y = reshape(y,prod(m),[]);
  Z = zeros(prod(m+1),size(y,2));         % allocate memory
  for i=1:size(y,2)                       % run over all components of Y
    yi = reshape(y(:,i),m);
    for j=1:dim                           % run over all dimensions in yi
      % J = [j,setdiff(1:dim,j)];         % make j-th dimension first
      J = Js{j};                          % make j-th dimension first
      yi = permute(yi,J);
      
      zi = yi([1 1:end],:,:);
      zi(2:end-1,:,:) = zi(2:end-1,:,:) + yi(2:end,:,:);
      
      yi = .5* ipermute(zi,J);            % undo permutation
    end
    Z(:,i) = yi(:);                       % store i-th component
  end
else
  % -----------------------
  % nodal -> cell-centered
  % -----------------------
  % MATLAB implementation
  y = reshape(y,prod(m+1),[]);
  Z = zeros(prod(m),size(y,2));           % allocate memory
  for i=1:size(y,2)                       % run over all components of Y
    if dim==1                             % 1D case is simplr and does not admit oermutations
      P = 0.5*(y(1:end-1)+y(2:end));
      return
    else
      yi = reshape(y(:,i),m+1);           % reorganize Y
    end
    for j=1:dim                           % run over all dimensions in yi
      % J = [j,setdiff(1:dim,j)];         % make j-th dimension first
      J = Js{j};                          % make j-th dimension first
      yi = permute(yi,J);
      % average nodal to center
      yi = 0.5*(yi(1:end-1,:,:)+yi(2:end,:,:));
      yi = ipermute(yi,J);                % undo permutation
    end
    Z(:,i) = reshape(yi,[],1);            % store i-th component
  end
  
end
P = reshape(Z,[],1);                      % reshape
%------------------------------------------------------------------------------

function runMinimalExample

omega = [0 2 0 1]; 
m     = [8,7];
yn    = getNodalGrid(omega,m);
xc    = nodal2center(yn,m);
xc    = reshape(xc,[],2);

FAIRfigure(1); clf;
subplot(2,1,1); spy(nodal2center([4,5,6])); title('spy(nodal2gcenter operator)');
subplot(2,1,2); plotGrid(yn,omega,m,'color','b'); hold on; plot(xc(:,1),xc(:,2),'rx');
title(mfilename);

%==============================================================================
