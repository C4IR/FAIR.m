%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
% function [ph,th] = plotLM(LM,varargin);
%
% a convenient tool for plotting landmarks
%                   
% Input
%   LM            landmark positions
%   varargin    optional parameter, see below
%
% Output:
%  ph           handle to plot
%  th           handle to text
%
% see also E5_Hands_TPS for an example
%==============================================================================

function [ph,th] = plotLM(LM,varargin)

if nargin==0
    help(mfilename)
    runMinimalExample; 
    return;
end

% setup default parameter
dim       = size(LM,2);
numbering = 'off';
dx        = 0.1;
fontsize  = 20;
color     = 'y';
for k=1:2:length(varargin), % overwrites default parameter
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

ph = []; th = []; % create empty handles

% handle options 'numbering', 'fontsize' and 'dx'
J = strcmp(varargin,'numbering') ...
  | strcmp(varargin,'dx') | strcmp(varargin,'fontsize');
if any(J),
  I = find(J);
  K = setdiff(1:length(varargin),[I,I+1]);
  varargin = {varargin{K}};
end;

% use plot in 2D or 3D
varargin = {varargin{:},'color',color};
switch dim,
  case 2, ph = plot(LM(:,1),LM(:,2),'x',varargin{:});
  case 3, ph = plot3(LM(:,1),LM(:,2),LM(:,3),'x',varargin{:});
end;

if strcmp(numbering,'off'), return;  end;

% add labels to the landmarks
dx = [dx,zeros(1,dim-1)];
pos = @(j) LM(j,:) + dx;
for j=1:size(LM,1),
  th(j) = text('position',pos(j),'string',sprintf('%d',j));
end;
set(th,'fontsize',fontsize,'color',color);
%------------------------------------------------------------------------------
function runMinimalExample
help(mfilename);
clear
close all
setup2DhandData
figure;
subplot(1,2,1);
viewImage2D(dataT,omega,m,'colormap','bone(265)');
hold on; ph =plotLM(LM(:,1:2),'color','g'); hold off;
set(ph,'marker','s')
subplot(1,2,2);
viewImage2D(dataR,omega,m,'colormap','bone(265)');
hold on; ph = plotLM(LM(:,3:4),'color','r'); hold off;
set(ph,'marker','o')
%------------------------------------------------------------------------------

