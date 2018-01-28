%% this file contains the settings for the HCS workflow

fprintf('Loading options... ');

global options;

options.data.rawDataChannels = {'*ch1*','*ch2*','*ch4*'};

% location of folder where the plates are saved
options.data.rawDataDir = 'd:\NEUBIAS\test_data\1_rawdata';

existingFolder = isdir(options.data.rawDataDir);

if ~existingFolder
    
    choice = questdlg(sprintf('The input data was not found in %s', options.data.rawDataDir), ...
        'Select raw data folder', ...
        'Browse','Cancel','Browse');
    
    switch choice
        case 'Browse'
            options.data.rawDataDir = uigetdir(options.HCSWorkflowPath,'Please select input folder where raw data is');
            if ~options.data.rawDataDir
                error('HCS:input','Input data was not found. Please set "options.data.rawDataFolder" in this file.' );
            end
            
            existingFolder = isdir(options.data.rawDataDir);
            
            if existingFolder
                plateList = dir(options.data.rawDataDir);
                plateList(~cat(1,plateList.isdir)) = [];
                for i=length(plateList):-1:1
                    if strcmp(plateList(i).name,'.') || strcmp(plateList(i).name,'..')
                        plateList(i) = [];
                    end
                end
                emptyInputFolder = isempty(plateList);
                if emptyInputFolder
                    error('HCS:loadsetting','Input folder was empty. Please set "options.data.rawDataFolder" in step_0_loadOptions.m file.' );
                end
            end
        otherwise
            error('HCS:loadsettings','Input data was not found. Please set "options.data.rawDataFolder" in step_0_loadOptions.m file.' );
    end
end

[mainpath,name,~] = fileparts(options.data.rawDataDir);
if isempty(name)
    [mainpath,name,~] = fileparts(mainpath);
end

options.data.correctedDataDir = fullfile(mainpath, '2_corrected');
if ~exist(options.data.correctedDataDir,'dir')
    mkdir(options.data.correctedDataDir);
end
options.data.focusedDataDir = fullfile(mainpath, '3_focused');
if ~exist(options.data.focusedDataDir,'dir')
    mkdir(options.data.focusedDataDir);
end
options.data.segmentedDataDir = fullfile(mainpath, '4_segmented');
if ~exist(options.data.segmentedDataDir,'dir')
    mkdir(options.data.segmentedDataDir);
end
options.data.qcDiscardedDir = fullfile(mainpath, '5_qcDiscarded');
if ~exist(options.data.qcDiscardedDir,'dir')
    mkdir(options.data.qcDiscardedDir);
end
options.data.machineLearningDir = fullfile(mainpath, '6_annotations');
if ~exist(options.data.machineLearningDir,'dir')
    mkdir(options.data.machineLearningDir);
end
options.data.statisticsDir = fullfile(mainpath, '7_reports');
if ~exist(options.data.statisticsDir,'dir')
    mkdir(options.data.statisticsDir);
end

options.data.allowedFileTypes = {'.bmp', '.gif', '.jpg', '.jpeg', '.tif', '.tiff', '.png',...
    '.BMP', '.GIF', '.JPG', '.JPEG', '.TIF', '.TIFF', '.PNG'};

% options.refactor.doRefactor = 0;
% options.refactor.inputFormat = [];
% options.refactor.outputFormat = '%s_w%c%02d_s%02d_p%02'; % SSS Structure

% available flat field correction method: 'cidre'
options.ffc.name = 'cidre';
options.ffc.vizualize = 1;

% available focus detection methods: 'bestplane' and 'adaptive'
options.focusDetection.name = 'bestplane';

options.segmentation.name = 'cp1_identify_objects';
% you can change the path of the pipeline
if ismac % I tried to work with a different segmentation method to avoid Mac error
    options.segmentation.pipelineName = fullfile(options.HCSWorkflowPath,'segmentation','pipelines','pipe_batch_Carraher_nucleiExtensionForMac_20180128.mat.mat');
else % the basic method that worked on windows and linux
    options.segmentation.pipelineName = fullfile(options.HCSWorkflowPath,'segmentation','pipelines','pipe_batch_Carraher_20180127.mat');
end
options.segmentation.poolSize = 4;
options.segmentation.batchSize = 4;
options.features = {};

fprintf('DONE.\n');
