
function varargout = viewImage2D(T,omega,m,varargin)

if nargin==0
  help(mfilename);
    runMinimalExample; 
    return;
end

if ~exist('m',    'var'), 
    m     = size(T); 
end;
if ~exist('omega','var'), 
    omega = reshape([ones(1,length(m));m],1,[]); 
end;

h  = (omega(2:2:end)-omega(1:2:end))./m;
xi = @(i) (omega(2*i-1)+h(i)/2:h(i):omega(2*i)-h(i)/2)';

ih  = image(xi(1),xi(2),reshape(T,m)'); axis xy image

% the following lines add some nice stuff to the code.
% if varargin = {'title','FAIR','xlabel','x'}
% the code evaluates "title('FAIR');xlabel('x');"
for k=1:2:length(varargin), 
  if ~isempty(varargin{k}), feval(varargin{k},varargin{k+1}); end;
end;

if nargout > 0, varargout = {ih};  end;

%------------------------------------------------------------------------------

function runMinimalExample

setup2DhandData;
FAIRfigure(1); clf;
viewImage2D(dataT,omega,m,'colormap','bone(256)');

%==============================================================================

