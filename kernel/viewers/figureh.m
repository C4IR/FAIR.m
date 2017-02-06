% amenity: open figure but do not shift window focus to the figure
function varargout = figureh(h)

if nargin == 0,
  h = figure;
  set(h,'position',position(h));
elseif not(ishandle(h))
  h = figure(h);
  set(h,'position',position(h));
else
  set(0, 'CurrentFigure', h);
end
if nargout > 0,
  varargout{1} = h;
end;


function pos = position(h)
if not(isnumeric(h))
  h = h.Number;
end;
pos = get(0,'MonitorPositions');
dx = mod((h-1)*25,floor(pos(1,4)/2));
dy = mod((h-1)*25,floor(pos(1,3)/2));
pos = round([pos(1,1)+dx,0.6*pos(1,4)-dy,0.33*pos(1,3),0.33*pos(1,4)]);
