clear all;

warning off all;

% CellProfiler pipeline name
pipeLineName = 'S8_SCREEN_PIPE_20171117_spot_detection_mlgoc.mat';

% number of CPUs working parallel
poolSize = 6;

% number of images in each batch
batchSize = 2;

% number of new attempts if a processes fails
numAttempts = 4;

mainFolder = 'd:\Projects\S8 SCREEN\plates\';

plateList = dir(fullfile(mainFolder, '1001*'));

% filter only on directories
plateList(~cat(1,plateList.isdir)) = [];

for i=1:numel(plateList)

    folderName = fullfile(mainFolder, plateList(i).name);
    failedJobs = evalBatch(folderName, pipeLineName, poolSize, batchSize, numAttempts);

end
