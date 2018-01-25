function focusDetection()
% focus detection

global options;

fprintf('\n====================\nFocus detection started...\n');

% setting default focus detection method if it does not exist
if ~isfield(options, 'focusDetection') || isempty(options.focusDetection)
    options.focusDetection = 'adaptive';
end

switch options.focusDetection
    case 'adaptive'
        focusMethod = @adaptiveFocus;
    case 'bestplane'
        focusMethod = @bestPlaneFocus;
    otherwise
        error('HCS:focus','No focus method %s is implemented.',options.focusDetecion);
end

for i=1:length(options.data.rawDataChannels)
%     imageList = dir(fullfile(options.data.correctedDataDir, options.data.rawDataChannels{i}));
    for r=1:16
        for c=1:24
            for f=1:1000
                imageRegexp = sprintf('%s_w%c%02d_s%02d_p*%s',plateName,r+'A'-1,c,f, options.data.rawDataChannels{i});
                stackList = dir(fullfile(options.data.refactoredDataDir,imageRegexp));
                for j = 1:length(stackList)
                    planeImagesStack(:,:,j) = imread(fullfile(options.data.refactoredDataDir,stackList(i).name));
                end
                focusedImage = focusMethod(planeImagesStack);
                outImageName = sprintf('%s_w%c%02_s%02d_%s.tiff',plateName,r+'A'-1,c,f,channelRegexp);
                imwrite(focusedImage, fullfile(options.data.focusedDataDir, outImageName));
            end
        end
    end
end

fprintf('\nFocus detection DONE.\n');

