%%% LamdbaTilde a scaling parameter for active contour model
LambdaTilde = 1;

%%% Rhatstar energy normalization parameter, its values have to be between
%%% 0.69 and 0.78 for positive alpha and beta values
Rhatstar = 0.75;

%%% ContourParameters contains the active contour, phase field and MRF
%%% parameters of the GOC inflection point model
MRGOCIPMParameters = computeMRGOCIPMparameters(LambdaTilde, Radius, Rhatstar);

%%% PriorPhasefieldParameters contains the phase field GOC parameters
PriorPhasefieldParameters = MRGOCIPMParameters(2);
PriorPhasefieldParameters.discrete = 1;
PriorPhasefieldParameters.marie = 1;

SegmentedGrayScaleObjects = CPretrieveimage(handles,['Segmented', SeedName],ModuleName,'DontCheckColor','DontCheckScale',size(OrigImage));

muin = mean2(OrigImage(SegmentedGrayScaleObjects>0));
sigmain = std2(OrigImage(SegmentedGrayScaleObjects>0));

muout = mean2(OrigImage(SegmentedGrayScaleObjects==0));
sigmaout = std2(OrigImage(SegmentedGrayScaleObjects==0));

MaskSegmentedObject = SegmentedGrayScaleObjects>0;
SE = strel('disk', round(Radius*0.2), 0);
MaskSegmentedObject = imerode(MaskSegmentedObject, SE);

%%% GradientWeight is the positive weight of image gradient term. It is not
%%% used in adaptive model
GradientWeight = 0.0;
% DataParameters = struct('muin', 0.15, 'sigmain', 0.01, 'muout', 0.03, 'sigmaout', 0.01, 'win', 1, 'wout', 1, 'gamma1', GradientWeight, 'gamma2', DataWeight);
DataParameters = struct('muin', muin, 'sigmain', sigmain, 'muout', muout, 'sigmaout', sigmaout, 'win', 1, 'wout', 1, 'gamma1', GradientWeight, 'gamma2', DataWeight);

Maxd = int32(max([PriorPhasefieldParameters.d]));

ExtendedImage = extendImage(OrigImage, Maxd, DataParameters(1).muout);
[hExtended, wExtended] = size(ExtendedImage);

% rng('shuffle');

if ~isnumeric(LayerNumString)
    % sort objects to layers
%     InitialPhi = sortGrayscaleObjects2Layers(SegmentedGrayScaleObjects.*MaskSegmentedObject, Radius);
    
    InitialPhi = sortGrayscaleObjects2LayersColoring(SegmentedGrayScaleObjects.*MaskSegmentedObject, 4);
    
end

OptimizationParameters = struct('maxIts', MaxIterations, 'saveFreq', -1, 'tolerance', Tolerance);

%%% Kappa is the weight of overlap penalty - not used with the adaptive
%%% model
Kappa = 0.0;

%%% Initialization can be Seeded, Neutral or Squared
if strcmp(Initialization, 'Seeds (manual)')
%     Initialization = 'manual';
    LayerNumber = size(InitialPhi, 3);
    ExtendedInitialPhi = zeros(hExtended, wExtended, LayerNumber) - 1;
    ExtendedInitialPhi(Maxd+1:Maxd+hOriginal, Maxd+1:Maxd+wOriginal, :) = (InitialPhi*2-1)*0.1;
elseif strcmp(Initialization, 'Neutral')
%     Initialization = 'neutral';
    LayerNumber = str2num(LayerNumString);
    ExtendedInitialPhi = randn(hExtended,wExtended,LayerNumber) + ones(hExtended,wExtended,LayerNumber)*PriorPhasefieldParameters(1).alpha/PriorPhasefieldParameters(1).lambda;
elseif strcmp(Initialization, 'Squares')
%     Initialization = 'squares';
    LayerNumber = str2num(LayerNumString);
    %change to squares!!!
    [InitialPhi] = createInitMLPhiSquares(hOriginal, wOriginal, PriorPhasefieldParameters, LayerNumber);
    ExtendedInitialPhi = zeros(hExtended, wExtended, LayerNumber) - 1;
    ExtendedInitialPhi(Maxd+1:Maxd+hOriginal, Maxd+1:Maxd+wOriginal, :) = (InitialPhi*2-1)*0.1;
else
%     Initialization = 'neutral';
    LayerNumber = str2num(LayerNumString);
    ExtendedInitialPhi = randn(hExtended,wExtended,LayerNumber) + ones(hExtended,wExtended,LayerNumber)*PriorPhasefieldParameters(1).alpha/PriorPhasefieldParameters(1).lambda;
end

drawnow

ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);
if any(findobj == ThisModuleFigureNumber)
    %%% Activates the appropriate figure window.
    CPfigure(handles,'Image',ThisModuleFigureNumber);
    if handles.Current.SetBeingAnalyzed == handles.Current.StartingImageSet
        CPresizefigure(OrigImage,'TwoByTwo',ThisModuleFigureNumber)
    end
    %%% A subplot of the figure window is set to display the original image.
    h11 = subplot(2,2,1);
    CPimagesc(OrigImage,handles);
    title(['Input Image, cycle # ', num2str(handles.Current.SetBeingAnalyzed)]);
    
    %%% A subplot of the figure window is set to display the input objects
    %%% for correction
    h12 = subplot(2,2,2);
    im = CPlabel2rgb(handles,SegmentedGrayScaleObjects.*MaskSegmentedObject);
    CPimagesc(im,handles);
    title({'Masked Segmented Objects'; ['#num of generated layers: ' num2str(LayerNumber)]});
end

drawnow

PriorPhasefieldParameters = repmat(PriorPhasefieldParameters, LayerNumber);

finalPhi = MLGOCSegmentationGM(PriorPhasefieldParameters, ExtendedImage, DataParameters, Kappa, ExtendedInitialPhi, OptimizationParameters);

Threshold = zeros(LayerNumber,1);
for ll=1:LayerNumber
    Threshold(ll) = PriorPhasefieldParameters(ll).alpha/PriorPhasefieldParameters(ll).lambda;
end

%MISSING: filter fully overlapping objects out 
% throw out objects that is fully covered by a bigger object, or the
% intercesction of 2 objects is more than a fixed amount of their union

[ContourImage, Contours] = createContourImage(OrigImage, finalPhi, Threshold);

%%%%%%%%%%%%%%%%%%%%%%%
%%% DISPLAY RESULTS %%%
%%%%%%%%%%%%%%%%%%%%%%%
% drawnow
% 
% ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);
%%% Check whether that figure is open. This checks all the figure handles
%%% for one whose handle is equal to the figure number for this module.
if any(findobj == ThisModuleFigureNumber)
    %%% Activates the appropriate figure window.
%     CPfigure(handles,'Image',ThisModuleFigureNumber);
%     if handles.Current.SetBeingAnalyzed == handles.Current.StartingImageSet
%         CPresizefigure(OrigImage,'TwoByTwo',ThisModuleFigureNumber)
%     end
%     %%% A subplot of the figure window is set to display the original image.
%     h11 = subplot(2,2,1);
%     CPimagesc(OrigImage,handles);
%     title(['Input Image, cycle # ',num2str(handles.Current.SetBeingAnalyzed)]);
%     
%     %%% A subplot of the figure window is set to display the input objects
%     %%% for correction
%     h12 = subplot(2,2,2);
%     im = CPlabel2rgb(handles,SegmentedGrayScaleObjects.*MaskSegmentedObject);
%     CPimagesc(im,handles);
%     title('Masked Segmented Objects');
    
    %%% A subplot of the figure window is set to display the sum of final
    %%% phase field layers
    h21 = subplot(2,2,3);
    CPimagesc(sum( (finalPhi+1)/2, 3),handles); 
    title('Sum of final layers');
    
    %%% A subplot of the figure window is set to display the outlined image
    h22 = subplot(2,2,4);
    CPimshow(ContourImage,handles); 
    title('Segmentation');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SAVE DATA TO HANDLES STRUCTURE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

%%% The Rescaled image is saved to the handles structure so it can be
%%% used by subsequent modules.
handles.Pipeline.(['MLGOCFinalPhi' ObjectName]) = finalPhi;
handles.Pipeline.(['MLGOCThreshold' ObjectName]) = Threshold;

SegmentedLayers = zeros(size(finalPhi));

for ll=1:LayerNumber
    
    SegmentedLayer = bwlabel( finalPhi(:,:,ll) > Threshold(ll) );
    if ll==1
%         MaxId = max( reshape(SegmentedLayer,1,[]) );
        SegmentedLayers(:,:,ll) = SegmentedLayer;
    else
        MaxId = max( reshape(SegmentedLayers(:,:,ll-1), 1, []));
        SegmentedLayer(SegmentedLayer>0) = SegmentedLayer(SegmentedLayer>0) + MaxId;
        SegmentedLayers(:,:,ll) = SegmentedLayer;
    end
end
handles.Pipeline.(['Segmented' ObjectName]) = SegmentedLayers;
handles.Pipeline.(['UneditedSegmented' ObjectName]) = SegmentedLayers;
handles.Pipeline.(['SmallRemovedSegmented' ObjectName]) = SegmentedLayers;