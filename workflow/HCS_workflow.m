%% main file 
% prepare a GUI?
global options;

p = mfilename('fullpath');
[HCSWorkflowPath,~,~] = fileparts(p);
addpath(genpath(HCSWorkflowPath));

options.HCSWorkflowPath = HCSWorkflowPath;

%% set options

loadOptions();

%% refactoring image names
% convert image names to the SSS structure if necessary

refactorImages();


%% flat field correction

% a) CIDRE
% b) mean image estimation

flatFieldCorrection();

%% focus detection

% a) adaptive focus
% b) focus plane selection
% c) extended depth of field?

focusDetection();

%% segmentation & feature extraction
% showing CellProfiler GUI, run batch version

% 1. nuclei segmentation: problems, low to high level algorithms
% 2. cell body segmentation: distance and intensity based methods

% - area/shape features
% - intensity features
% - texture features

segmentation();

if ~options.segmentation.isFeatureExtractionIncluded
    featureExtraction();
end

%% field of view image quality control

filterOnImageQuality();

%% cell-level quality control

filterOnSegmentationQuality(options);