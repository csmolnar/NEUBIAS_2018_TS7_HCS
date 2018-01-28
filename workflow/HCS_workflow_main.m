%% main file 
%
% Run this file
% Please set "options.data.rawDataDir" variable in file
% step_0_loadOptions.m file.
% The raw data of plates required in the following structure:
% 
% rawDatadir/
%       plateName1/
%       plateName2/
%       ...
%
% All the other folders will be generated automatically.

global options;

p = mfilename('fullpath');
[HCSWorkflowPath,~,~] = fileparts(p);
addpath(genpath(HCSWorkflowPath));

options.HCSWorkflowPath = HCSWorkflowPath;

%% set options
% setting all the necessary variables to personalize each component in the
% workflow

step_0_loadOptions();

%!!! the following part does not work without running all the code above !!!

%% refactoring image names
% convert image names to the SSS structure if necessary

% under construction and bugfixing
% step_0_1_refactorImages();


%% flat field correction
% CIDRE illumination correction is implemented

step_1_flatFieldCorrection();

%% focus detection
% creating one image from aquired images of different focus planes for each
% field of view

% a) adaptive focus: 'adaptive'
% b) focus plane selection: 'bestplane'

step_2_focusDetection();

%% field of view image quality control

% under construction and bugfixing
% step_3_fovQualityControl();

%% segmentation & feature extraction
% 
% For segmentation of trivial cases, which can be solved by threshold,
% we use the base version of CellProfiler. In case of more complex problems
% preprocessing steps (e.g. ilastik) or special modules are included in the
% workflow.
%
% To use you own pipeline please run CellProfiler and save the pipeline and
% set "options.segmentation.pipelineName" variable in step_0_loadOptions.m.
% 
% 1. nuclei segmentation: problems, low to high level algorithms
% 2. cell body segmentation: distance and intensity based methods

% 3. feature extraction
%   - area/shape features
%   - intensity features
%   - texture features

% ignore warning messages

step_4_segmentation_featureExtraction();

%% cell-level quality control

% step_5_cellLevelQualityControl();

%% annotation
% 
% To discover and analysis of extracted cell image data we use Advanced
% Cell Classifier.
% Source code is available here: http://www.cellclassifier.org/download/
% Start the software by running startup.m in ACC source folder.
%
% You can find tutorial videos about getting started and all
% functionalities here: http://www.cellclassifier.org/about-acc/.
% 
% In case of any troubles with the preceding steps sample dataset is
% available on the download page.
%
% Getting started with ACC
%
% The ACC file structure is the following:
%
% segmentedDir/
%       plateName1/
%           anal1/
%               plateName1_wA01_s01*.ext
%               plateName1_wA02_s01*.ext
%               ...
%           anal2/
%               plateName1_wA01_s01*.txt
%               plateName1_wA03_s01*.txt
%               ...
%           anal3/
%               plateName1_wA01_s01*.ext
%               plateName1_wA03_s01*.ext
%               ...
%
% where anal1 folder contains the images for 'contour' view, anal3 folder
% contains images for 'colour' (='not contour' view), and anal2 folder
% contains text files for the extracted features for each cell line by
% line (rows: cells, columns: features, column number 1 and 2 is for
% centroid of the detected cells). anal2 folder also contains the file that
% contains the names of features each in separate lines: featureNames.acc.
% The content of anal2 folder can be exported by the exporttoacc modules of 
% CellProfiler.
%
% In this workflow the 4_segmented folder should contain the data.
%  
% 6_annotations folder is created to store the already labelled cells by
% Save project option.
% 

%% classification
% 
% In ACC you can press "Predict selected plates" button for classification
% of each cell in the selected plates with some base statistics.


%% statistics
%
% For predicted plates it is possible to report several statistics from
% extracted features under Classification->Feature based statistics
