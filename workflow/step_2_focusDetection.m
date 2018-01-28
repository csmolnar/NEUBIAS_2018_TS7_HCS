function step_2_focusDetection()
% focus detection

global options;

fprintf('\n====================\nFocus detection started...\n');

% setting default focus detection method if it does not exist
if ~isfield(options, 'focusDetection') || isempty(options.focusDetection)
    options.focusDetection = 'bestplane';
end

switch options.focusDetection.name
    case 'adaptive'
        focusMethod = @adaptiveFocus;
    case 'bestplane'
        focusMethod = @selectFocusPlane;
    otherwise
        error('HCS:focus','No focus method %s is implemented.',options.focusDetecion);
end

plateList = dir( fullfile(options.data.rawDataDir));

plateList(~cat(1,plateList.isdir)) = [];

for i=length(plateList):-1:1
    if strcmp(plateList(i).name,'.') || strcmp(plateList(i).name,'..')
        plateList(i) = [];
    else
        mkdir(fullfile(options.data.focusedDataDir, plateList(i).name) );
    end
end

for ii=1:length(plateList)
    fprintf('Focus detection of plate %s is in progress...\n', plateList(ii).name);
    for i=1:length(options.data.rawDataChannels)
        baseChannelName = options.data.rawDataChannels{i};
        baseChannelName(baseChannelName=='*') = '';
        % r is the row number in plate
        for r=1:16
            % c is column number in plate
            for c=1:24
                % f is the number of field of view
                for f=1:4
                    imageRegexp = sprintf('%s_w%c%02d_s%02d%s_p*',plateList(ii).name,r+'A'-1,c,f, options.data.rawDataChannels{i});
                    stackList = dir(fullfile(options.data.correctedDataDir,plateList(ii).name,imageRegexp));
                    if ~isempty(stackList)
                        for j = 1:length(stackList)
                            planeImagesStack(:,:,j) = imread(fullfile(options.data.correctedDataDir,plateList(ii).name,stackList(j).name));
                        end
                        focusedImage = focusMethod(planeImagesStack);
                        outImageName = sprintf('%s_w%c%02d_s%02d_%s.tiff',plateList(ii).name,r+'A'-1,c,f,baseChannelName);
                        imwrite(focusedImage, fullfile(options.data.focusedDataDir, plateList(ii).name, outImageName));
                    end
                end
            end
        end
    end
    
end

fprintf('\nFocus detection DONE.\n');

