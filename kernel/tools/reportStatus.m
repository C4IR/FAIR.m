%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
% reports FAIRs's current setting
%==============================================================================

caller = dbstack; % identify the name of the calling function
caller = caller(min(length(caller),2)).name;
FAIRmessage(caller,'-')
viewImage('disp');
imgModel('disp');
distance('disp');
trafo('disp');
regularizer('disp');
FAIRmessage('-');
%==============================================================================
