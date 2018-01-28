%% main file 
% Run this pipeline
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

step_0_loadOptions();

%% refactoring image names
% convert image names to the SSS structure if necessary

% not double checked
% step_0_1_refactorImages();


%% flat field correction

% a) CIDRE
% b) mean image estimation

step_1_flatFieldCorrection();

%% focus detection

% a) adaptive focus
% b) focus plane selection
% c) extended depth of field?

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

step_4_segmentation_featureExtraction();

%% cell-level quality control

% step_5_cellLevelQualityControl();

%% annotation
% 
% To discover and analysis of extracted cell image data we use Advanced
% Cell Classifier
% Source code is available here: http://www.cellclassifier.org/download/
% Start the software by running startup.m in ACC source folder
%
% You can find tutorial videos about getting started and all
% functionalities here: http://www.cellclassifier.org/about-acc/
% 
% In case of any troubles with the preceding steps sample dataset is
% available on the download page
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
