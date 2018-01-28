function prepareBatchAnalysis(pipeLineName, imageFolder, outFolder)

p = mfilename('fullpath');
[scriptPath,~,~] = fileparts(p);

load(fullfile(scriptPath,'temphandles.mat'),'handles');

% pipeLineName
% imageFolder

subFolderList = {'anal1', 'anal2', 'anal3', 'anal4', 'anal5'};

for i=1:length(subFolderList)
    if exist(fullfile(outFolder, subFolderList{i}), 'dir')
        delete(fullfile(outFolder, subFolderList{i},'*.png'));
    else
        mkdir(fullfile(outFolder, subFolderList{i}));
    end
end

% create folders, remove analysis content if exist
% if exist(fullfile(imageFolder, 'anal1'), 'dir')    
%     delete(fullfile(imageFolder, 'anal1','*.png'));
% else
%     mkdir(fullfile(imageFolder, 'anal1'));
% end
% 
% if exist(fullfile(imageFolder, 'anal2'), 'dir')    
%     delete(fullfile(imageFolder, 'anal2','*.png'));
% else
%     mkdir(fullfile(imageFolder, 'anal2'));
% end
% 
% if exist(fullfile(imageFolder, 'anal3'), 'dir')    
%     delete(fullfile(imageFolder, 'anal3','*.png'));
% else
%     mkdir(fullfile(imageFolder, 'anal3'));
% end
% 
% if exist(fullfile(imageFolder, 'anal4'), 'dir')    
%     delete(fullfile(imageFolder, 'anal4','*.png'));
% else
%     mkdir(fullfile(imageFolder, 'anal4'));
% end
% 
% if exist(fullfile(imageFolder, 'anal5'), 'dir')    
%     delete(fullfile(imageFolder, 'anal5','*.png'));
% else
%     mkdir(fullfile(imageFolder, 'anal5'));
% end

handles.Current.DefaultImageDirectory = imageFolder;

% fake the figure handles
handles.Current.FigureNumberForModule02 = 111;
handles.Current.FigureNumberForModule03 = 111;
handles.Current.FigureNumberForModule04 = 111;
handles.Current.FigureNumberForModule05 = 111;
handles.Current.FigureNumberForModule06 = 111;
handles.Current.FigureNumberForModule07 = 111;
handles.Current.FigureNumberForModule08 = 111;
handles.Current.FigureNumberForModule09 = 111;
handles.Current.FigureNumberForModule10 = 111;
handles.Current.FigureNumberForModule11 = 111;
handles.Current.FigureNumberForModule12 = 111;
handles.Current.FigureNumberForModule13 = 111;
handles.Current.FigureNumberForModule14 = 111;
handles.Current.FigureNumberForModule15 = 111;
handles.Current.FigureNumberForModule16 = 111;
handles.Current.FigureNumberForModule17 = 111;
handles.Current.FigureNumberForModule18 = 111;
handles.Current.FigureNumberForModule19 = 111;
handles.Current.FigureNumberForModule20 = 111;
handles.Current.FigureNumberForModule21 = 111;
handles.Current.FigureNumberForModule22 = 111;
handles.Current.FigureNumberForModule23 = 111;
handles.Current.FigureNumberForModule24 = 111;
handles.Current.FigureNumberForModule25 = 111;
handles.Current.FigureNumberForModule26 = 111;
handles.Current.FigureNumberForModule27 = 111;
handles.Current.FigureNumberForModule28 = 111;
handles.Current.FigureNumberForModule29 = 111;
handles.Current.FigureNumberForModule30 = 111;
handles.Current.FigureNumberForModule31 = 111;
handles.Current.FigureNumberForModule32 = 111;
handles.Current.StartingImageSet = 1;

% set I/O folders
handles.Preferences.DefaultImageDirectory = imageFolder;
handles.Preferences.DefaultOutputDirectory = outFolder;
handles.Current.DefaultImageDirectory = imageFolder;
handles.Current.DefaultOutputDirectory = outFolder;

handles.Settings = load(pipeLineName);

handles.Settings = handles.Settings.Settings;

handles = LoadImages(handles);

%handles.Pipeline.FilenameOrigImage = {};
%handles.Pipeline = rmfield(handles.Pipeline, {'OrigImage', });

handles.Current.NumberOfModules = length(handles.Settings.NumbersOfVariables);

if length(handles.Settings.NumbersOfVariables) < 10
    handles.Current.CurrentModuleNumber = ['0' num2str(length(handles.Settings.NumbersOfVariables))];
else
    handles.Current.CurrentModuleNumber = num2str(length(handles.Settings.NumbersOfVariables));
end

save([outFolder 'CPBatchInfo.mat'], 'handles');
