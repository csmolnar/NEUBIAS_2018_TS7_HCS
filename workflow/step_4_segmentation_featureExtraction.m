function step_4_segmentation_featureExtraction()
% batch segmentation of images

global options;

fprintf('\n====================\nSegmentation started..\n');

if ~isfield(options,'segmentation') || isempty(options.segmentation)
    error('HCS:segmentation','No information about segmentation exists');
end

switch options.segmentation.name
    case 'cp2_identify_objects' % segmentation via CellProfiler 2.2 (Python version)
        CP2segmentation();
    case 'cp1_identify_objects' % segmentation via CellProfiler 1.0 (Matlab version)
        CP1segmentation();
    otherwise
        error('HCS:segmentation', 'No segmentation method %s is implemented', options.segmentation.name);
end

fprintf('\nSegmentation DONE.\n');