%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% KERNEL/DATA
%
% This folder provides the basis data (examples) that come with the toolbox 
% and some convenient tools for handling the data;
% see also references at the end of the mfiles
%
% To get the data "X" (e.g. 2DHandData), conveniently run setup"X";
% e.g. setup2DHandData.
% For faster access, the generated data is saved in an automatically 
% generated file in the temp folder and later simply re-used.
%
%------------------------------------------------------------------------------
%
% The following examples are available:
%
% 2D data:
%   2Dhands     -  hand xrays,                  see ref. [1] and setup2DHandData
%   2DHNSP      -  brain serial sectioning,     see ref. [2] and setup2DHNSPData
%   2DMRIs      -  MR slices from human brain,  see ref. [3] and setup2DMRIData
%   2DPETCT     -  PET-CT study,                see ref. [4] and setup2DPETCTdata
%   2DUS        -  US image of Vibe Heldmann,
%                  courtesy by Dr. Stefan Heldmnn, Fraunhofer MEVIS, Lübeck
%   2Ddisc2C    -  binary images of a disc and a C a la Gary Christensen,
%                  test case for large deformations, academic testcase
%   2DGaussian  -  Gaussian blobs, test-case for mass-preserving schemes,
%                  academic testcase
%
% 3D data:
%   3Dbrains    -  MRI data of a human brain,
%                  data courtesy Ron Kikinis, Brigham & Women's Hospital, Boston, USA
%   3Dknees     -  MRI data from a human knee,
%                  data courtesy Thomas Netsch, Philips Medical Solutions, Hamburg
%   3Dmice      -  PET data of a mouse heart,
%                  data courtesy European Institute for Molecular Imaging (EIMI),
%                  University of Muenster, Germany
%   3Dboxes     -  academical example with deformed binary boxes
%
% Tools (M-files)
%
%   setup2DhandData     - requires hands-T.jpg,     hands-R.jpg
%   setup2DHNSPData     - requires HNSP-T.jpg,      HNSP-R.jpg
%   setup2DMRIData      - requires MRIhead-T.jpg,   MRIhead-R.jpg
%   setup2DPETCTData    - requires PET-CT-CT.jpg,   PET-CT-PET.jpg
%   setup2DUSData       - requires US.jpg
%   setup2Ddisc2CData   - requires disc.jpg and c.jpg
%   setup2DGaussianData - requires Gauss-T.jpg, Gauss-r.jpg
%   setup3DbrainData    - requires brain3D.mat
%   setup3DkneeData     - requires createPhilipsKnees.mat
%   setup3DmouseData    - requires mice3D.mat
%   setup3DboxData      - academic, a box 
%
%------------------------------------------------------------------------------
% More tools:
%   contents            -  this file
%   checkDataFile       -  in the first run, the data is saved and later simply
%                          loaded for faster access,
%                          this file handles the administration
%   jpgs2data           -  convenient tool to convert two jpgs files into the
%                          FAIR data structure
%   testData            -  tests this part of the toolbox
%
%------------------------------------------------------------------------------
% JPG files:
%   Gauss-R.jpg,           Gauss-T.jpg
%   HNSP-R.jpg             HNSP-T.jpg
%   MRIhead-R.jpg          MRIhead-T.jpg
%   PET-CT-CT.jpg          PET-CT-PET.jpg
%   c.jpg                  disc.jpg
%   hands-R.jpg            hands-T.jpg
%   US.jpg
%   EPIslice-R.jpg         EPIslice-T.jpg [5]
%
% MAT files:
%   brain3D.mat
%   knees3D.mat
%   mice3D.mat
%
%
% [5] 2D slices of MR images, Institute of Clinical Radiology, Univesity
%     Hospital Muenster, Germany.
%   phantom3D.mat
%     MR data of a hardware phantom,
% Institute for Clinical Radiology, University Hospital Muenster, Germany
% for data sources see the corresponding setup file
%
%==============================================================================

function debit = contents
if nargout == 0, 
  help(mfilename); 
  return; 
end;

debit = {
  'contents.m'
  'checkDataFile.m'
  'jpgs2data.m'

  'setup2DGaussianData.m'
  'setup2DHNSPData.m'
  'setup2DMRIData.m'
  'setup2DPETCTData.m'
  'setup2DUSData.m'
  'setup2Ddisc2CData.m'
  'setup2DhandData.m'
  'setup3DboxData.m'
  'setup3DbrainData.m'
  'setup3DkneeData.m'
  'setup3DmouseData.m'

  'testData.m'

  'Gauss-R.jpg'
  'Gauss-T.jpg'
  'HNSP-R.jpg'
  'HNSP-T.jpg'
  'MRIhead-R.jpg'
  'MRIhead-T.jpg'
  'PET-CT-CT.jpg'
  'PET-CT-PET.jpg'
  'US.jpg'
  'c.jpg'
  'disc.jpg'
  'hands-R.jpg'
  'hands-T.jpg'
  'brain3D.mat'
  'EPIslice-R.jpg'
  'EPIslice-T.jpg'
  'knee3D.mat'
  'mice3D.mat'
  'phantom3D.mat'
  };
%==============================================================================
