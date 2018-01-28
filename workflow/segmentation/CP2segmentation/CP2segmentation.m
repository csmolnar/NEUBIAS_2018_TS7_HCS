function CP2segmentation()
% Batch segmentation via CellProfiler 2.2

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

% runc CP2 pipeline in batch mode
for i=1:numel(plateList)
    % run the CellProfiler from command line
    % system(sprintf('CellProfiler.exe -c -r -p %s -i %s -o %s --plugins-directory=%s --ij-plugins-directory=%s',...
%                 pipePath, inputDir, outputDir, cp2PluginDir, cp2IJPluginDir));
end

% 
% inputDir = 'c:\Users\mcsaba\Downloads\ExampleHuman\images\';
% outputDir = 'c:\Users\mcsaba\Downloads\ExampleHuman\images\';
% 
% cp2PluginDir = 'd:\Projects\CellProfiler2Plugins\cpplugins\';
% cp2IJPluginDir = 'd:\Projects\CellProfiler2Plugins\ijplugins\';
% 
% pipePath = 'd:\Projects\NEUBIAS\workflow\tempCode\Example_pipeline_project.cppipe';
% 
% system(sprintf('CellProfiler.exe -c -r -p %s -i %s -o %s --plugins-directory=%s --ij-plugins-directory=%s',...
%                 pipePath, inputDir, outputDir, cp2PluginDir, cp2IJPluginDir));