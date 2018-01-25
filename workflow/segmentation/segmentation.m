function segmentation()
% batch segmentation of images

global options;

fprintf('\n====================\nSegmentation started..\n');

if ~isfield(options,'segmentation') || isempty(options.segmentation)
    error('HCS:segmentation','No information about segmentation exists');
end

switch options.segmentation.name
    case 'cp_identify_objects'
        for i=1:length(options.segmentation.commands)
            system(options.segmentation.commands{i});
        end
    otherwise
        error('HCS:segmentation', 'No segmentation method %s is implemented', options.segmentation.name);
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

fprintf('\nSegmentation DONE.\n');