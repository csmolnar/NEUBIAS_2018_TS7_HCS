function fovQC()
% field of view quality controls

global options;

mainInputDir = options.data.focusedDataDir;
discardedFilesDir = options.data.qcDiscardedDir;

plateList = dir(fullfile(mainDir));

% clear not plate folders
plateList(~cat(1,plateList.isdir)) = [];

for i=length(plateList):-1:1
    if strcmp(plateList(i).name,'.') || strcmp(plateList(i).name,'..')
        plateList(i) = [];
    end
end

for i=1:length(plateList)
    
    
    
end




