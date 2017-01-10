%==============================================================================
% ##1
% ##2
%==============================================================================
%
% KERNEL/VIEWERS
% 
% Contents of FAIR VIEWER
%
% contents                  this file
% testViewer                tests the view components
% showResults               convenient way of visualizing results
% FAIRplots                 visualization for MLIR and such
%                           +----------------------------------+
%                           | [R]       [T(y0)]    [T(yc)]     |
%                           | [T+grid]  [T(y0)-R]  [T(yc)-R]   |
%                           +----------------------------------+
% FAIRfigure                opens a non-gray figure
% FAIRposition              computes positions for a figure
% figureh                   open figure but do not shift window focus to the figure
%
% viewImage                 standard (configurable) viewer for data
% viewImage2D               2D viewer for data (essentially image.m)
% viewImage2Dsc             2D viewer for data (essentially imagesc.m)
% viewSlices                3D viewer
% volView                   3D renderer
% imgmontage                montage 3D data
% overlayImage2D            view overlay of two 2D datasets
% viewIP 					visualizes the intensity projections of a 3D image
%                           (average/minimun/maximum)
% plotGrid                  plot a grid
% plotStaggeredGrid         plot a staggered grid
% plotIterationHistory      plot the iteration history of an optimization scheme
% plotMLIterationHistory    plot the iteration history of a multilevel scheme
%==============================================================================

function debit = contents

if nargout == 0, 
    help(mfilename); 
    return; 
end;

debit = {
  'contents.m'
  'testViewer.m'
  'showResults.m'
  'FAIRplots.m'

  'FAIRfigure.m'
  'FAIRposition.m'
  'figureh.m'

  'viewImage.m'
  'viewImage2D.m'
  'viewImage2Dsc.m'
  'viewSlices.m'
  'volView.m'
  'imgmontage.m'
  'overlayImage2D.m'
  'viewIP.m'
   
  'plotGrid.m'
  'plotIterationHistory.m'
  'plotMLIterationHistory.m'
  'plotStaggeredGrid.m'
  
};
%==============================================================================