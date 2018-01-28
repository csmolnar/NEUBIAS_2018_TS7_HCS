function CP1segmentation()
% Batch segmentation via CellProfiler 1.0

global options;

inputMainDir = options.data.focusedDataDir;
outputMainDir = options.data.segmentedDataDir;

plateList = dir(fullfile(inputMainDir));

pipelineName = options.segmentation.pipelineName;

% filter only on directories
plateList(~cat(1,plateList.isdir)) = [];

for i=length(plateList):-1:1
    if strcmp(plateList(i).name,'.') || strcmp(plateList(i).name,'..')
        plateList(i) = [];
    end
end

% runc CP pipeline in batch mode
for i=1:numel(plateList)

    inputFolderName = fullfile(inputMainDir, plateList(i).name);
    outputFolderName = fullfile(outputMainDir, plateList(i).name);
    failedJobs = evalBatch(inputFolderName, outputFolderName, pipelineName, options.segmentation.poolSize, options.segmentation.batchSize, 3);

end