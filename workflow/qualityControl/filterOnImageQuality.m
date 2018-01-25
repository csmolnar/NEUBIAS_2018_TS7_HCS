function filterOnImageQuality(filterType, inputDir, outputDir)
% quality control

global options;

switch filterType
    case 'oversaturation'
    case 'segmentationerror'
    otherwise
        error('HCS:qualitycontrol', 'No quality filter type %s is implemented.', filterType);        
end

end