%% main file 
% prepare a GUI?

p = mfilename('fullpath');
[HCSWorkflowPath,~,~] = fileparts(p);
addpath(genpath(HCSWorkflowPath))

%% set options

loadOptions();

%% field of view image quality control

filterOnImageQuality(options);

%% flat field correction

% a) CIDRE
% b) mean image estimation

flatFieldCorrection(options);

%% focus detection

% a) adaptive focus
% b) focus plane selection
% c) extended depth of field?

focusDetection(options);

%% segmentation & feature extraction
% showing CellProfiler GUI, run batch version

% 1. nuclei segmentation: problems, low to high level algorithms
% 2. cell body segmentation: distance and intensity based methods

% - area/shape features
% - intensity features
% - texture features

segmentation(options);


