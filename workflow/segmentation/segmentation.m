function segmentation(options)
% batch segmentation of images

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

fprintf('\nSegmentation DONE.\n');