function step_3_fovQualityControl()
% field of view quality controls

global options;

mainInputDir = options.data.focusedDataDir;
discardedFilesDir = options.data.qcDiscardedDir;
if ~exist(discardedFilesDir,'dir')
    mkdir(discardedFilesDir);
end

plateList = dir(fullfile(mainInputDir));

% clear not plate folders
plateList(~cat(1,plateList.isdir)) = [];

for i=length(plateList):-1:1
    if strcmp(plateList(i).name,'.') || strcmp(plateList(i).name,'..')
        plateList(i) = [];
    end
end

regexpr='*ch1.tif';
report=false;
numChannels=3;
thres=1.45;
blueThres=300000;
blueSize=Inf;
blueThresAdj=Inf;
adjust=true;
fprintf('running focus + stain error detection...\n');

runQCFcn(mainInputDir,plateList,discardedFilesDir,discardedFilesDir,regexpr,report,numChannels,thres,adjust,blueThres,blueSize,blueThresAdj);
