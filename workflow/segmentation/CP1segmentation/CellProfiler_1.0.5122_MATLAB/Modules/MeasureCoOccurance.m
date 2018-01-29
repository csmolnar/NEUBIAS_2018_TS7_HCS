function handles = MeasureCoOccurance(handles)

% Help for the Measure Correlation module:
% Category: Measurement
%
% SHORT DESCRIPTION:
% Measures the cooccurances on objects between two intensity channels.
% Given the upper and lower intensity limit, bins the rannge into an nxn
% matrix and calculates the co-occurances such that the sum of the matrix
% is normalized to 1.
% *************************************************************************
%

% CellProfiler is distributed under the GNU General Public License.
% See the accompanying file LICENSE for details.
%
% Developed by Peter Horvath
% ETH Zurich 2011
%
%
% $Revision: 5025 $

%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%
drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = Choose two images to measure cooccurance between (first image):
%infotypeVAR01 = imagegroup
ImageName{1} = char(handles.Settings.VariableValues{CurrentModuleNum,1});
%inputtypeVAR01 = popupmenu

%textVAR02 = Lower and upper bounds on the measurement, in the range [0,1]
%defaultVAR02 = 0,0.0625
Image1LowUpRange = char(handles.Settings.VariableValues{CurrentModuleNum,2});

%textVAR03 = Choose two images to measure cooccurance between (second image):
%infotypeVAR03 = imagegroup
ImageName{2} = char(handles.Settings.VariableValues{CurrentModuleNum,3});
%inputtypeVAR03 = popupmenu

%textVAR04 = Lower and upper bounds on the measurement, in the range [0,1]
%defaultVAR04 = 0,0.0625
Image2LowUpRange = char(handles.Settings.VariableValues{CurrentModuleNum,4});

%textVAR05 = Number of bins
%defaultVAR05 = 5
NumberOfBins = str2double(char(handles.Settings.VariableValues{CurrentModuleNum,5}));

%textVAR06 = Which objects do you want to measure?
%infotypeVAR06 = objectgroup
ObjectName = char(handles.Settings.VariableValues{CurrentModuleNum,6});
%inputtypeVAR06 = popupmenu

%%%VariableRevisionNumber = 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRELIMINARY CALCULATIONS & FILE HANDLING %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

%%% Get the images
if strcmp(ImageName{1},'') || strcmp(ImageName{2},'')
    error(['Image processing was canceled in the ', ModuleName, ' module because two image types must be chosen.'])
else
    try
        %%% Checks whether image has been loaded.
        Image1 = CPretrieveimage(handles,ImageName{1},ModuleName,'MustBeGray','DontCheckScale');
        Image2 = CPretrieveimage(handles,ImageName{2},ModuleName,'MustBeGray','DontCheckScale');
    catch
        error(['Image processing was canceled in the ', ModuleName, ' module because there was a problem loading the image.']);
    end
end


%%% Get the masks of segmented objects
if strcmp(ObjectName,'')
    error(['Image processing was canceled in the ', ModuleName, ' module because object type must be chosen.'])
else
    if ~strcmp(ObjectName,'Image')
        %%% Retrieves the label matrix image that contains the
        %%% segmented objects which will be used as a mask.
        LabelMatrixImage = CPretrieveimage(handles,['Segmented', ObjectName],ModuleName,'MustBeGray','DontCheckScale');
    else
        LabelMatrixImage = ones(size(Image1));
    end
    
end;

% Minimum and maximum bin ranges
% image1
index = strfind(Image1LowUpRange,',');
if isempty(index)
    error(['Image processing was canceled in the ', ModuleName, ' module because the Min and Max imag bounds are invalid.'])
end
Image1LowRange = str2double(Image1LowUpRange(1:index-1));
Image1UpRange  = str2double(Image1LowUpRange(index+1:end));

% image2
index = strfind(Image2LowUpRange,',');
if isempty(index)
    error(['Image processing was canceled in the ', ModuleName, ' module because the Min and Max imag bounds are invalid.'])
end
Image2LowRange = str2double(Image2LowUpRange(1:index-1));
Image2UpRange  = str2double(Image2LowUpRange(index+1:end));

%%%%%%%%%%%%%%%%%%%%%%
%%% IMAGE ANALYSIS %%%
%%%%%%%%%%%%%%%%%%%%%%
drawnow


% check if the images and label matrices are in the same size
if any(size(Image1) < size(LabelMatrixImage)) || any(size(Image2) < size(LabelMatrixImage))
    error(['Image processing was canceled in the ', ModuleName, ' module. The size of the image you want to measure is not the same as the size of the image from which the ',ObjectName{ObjectNameNbr},' objects were identified.']);
end;




%%% Calculate the cooccurance matrices for all objects
binImage1 = linspace(Image1LowRange, Image1UpRange, NumberOfBins+1);
binImage2 = linspace(Image2LowRange, Image2UpRange, NumberOfBins+1);
NbrOfObjects = max(LabelMatrixImage(:));          % Get number of segmented objects
cooccurance = zeros(NumberOfBins);
cooccuranceFeatureSet = zeros(NbrOfObjects, NumberOfBins^2);
for ObjectNbr = 1:NbrOfObjects
    index = find(LabelMatrixImage == ObjectNbr);   % Get the indexes for the this object number
    intensitiesImage1 = Image1(index);
    intensitiesImage2 = Image2(index);
    % find indices in image1 where the intensity is in bin_n
    
    for binIterator1=1:length(binImage1)-1
        indicesImage1 = intensitiesImage1>=binImage1(binIterator1) & intensitiesImage1<binImage1(binIterator1+1);
        % get a histogram on image2
        currentIntensitiesImage2 = intensitiesImage2(indicesImage1);
        for binIterator2=1:length(binImage2)-1
            cooccurance(binIterator1, binIterator2) = sum(currentIntensitiesImage2>=binImage2(binIterator2) & currentIntensitiesImage2<binImage2(binIterator2+1));
        end;
    end;
    % normalization
    if sum(cooccurance(:)) > 0
        cooccurance = cooccurance / sum(cooccurance(:));
    end;
    cooccuranceFeatureSet(ObjectNbr, :) = cooccurance(:);
end;
if strcmp (ObjectName,'Image')
    handles.Measurements.(ObjectName).CooccuranceFeatures = ['Cooccurence_bin_' num2str(NumberOfBins)];
    handles.Measurements.(ObjectName).Cooccurance(handles.Current.SetBeingAnalyzed) = {cooccuranceFeatureSet};
else
    handles.Measurements.(ObjectName).CooccuranceFeatures = ['Cooccurence_bin_' num2str(NumberOfBins)];
    handles.Measurements.(ObjectName).Cooccurance(handles.Current.SetBeingAnalyzed) = {cooccuranceFeatureSet};
end


%%%%%%%%%%%%%%%%%%%%%%%
%%% DISPLAY RESULTS %%%
%%%%%%%%%%%%%%%%%%%%%%%
drawnow

