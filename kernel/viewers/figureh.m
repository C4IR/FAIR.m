% amenity: open figure but do not shift window focus to the figure
function varargout = figureh(h)
if nargin>=1 
    if ishandle(h)
        set(0, 'CurrentFigure', h);
    else
        h = figure(h);
    end
else
    h = figure;
end
if nargout>0, varargout = {h}; end;
