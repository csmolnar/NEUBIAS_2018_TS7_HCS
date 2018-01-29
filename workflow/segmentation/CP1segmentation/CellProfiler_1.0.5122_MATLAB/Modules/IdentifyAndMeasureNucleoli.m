function handles = IdentifyAndMeasureNucleoli(handles)

% Help for the Identify Nucleoli module:
% Category: Object Processing
%
% SHORT DESCRIPTION:
% Identifies objects (e.g. nucleoli) inside "seed" objects identified by
% an Identify Primary module (e.g. nuclei).
% *************************************************************************
%
% This module identifies secondary objects (e.g. nucleoli) based on two
% inputs: (1) a previous module's identification of primary objects (e.g.
% nuclei) and (2) an image stained for the secondary objects (not required
% for the Distance - N option). Each secondary object is assumed to be completely
% within a primary object (e.g. nucleoli are completely within nuclei
% stained for actin).
%
% It accomplishes two tasks:
% (a) finding the dividing lines between secondary objects which touch each
% other. Three methods are available: Propagation, Watershed (an older
% version of Propagation), and Distance (but the method is not using it).
% (b) finding the dividing lines between the secondary objects and the
% background of the image. 
% (c) measures basic shape and intensity metric are implemented in
% MeasureObjectAreaShape and MeasureObjectIntensity modules.
% (d) relates each secondary object to exactly one primary object like
% 'parent' similarly to the Relate module.
%
% Settings:
%
% Ratio of overlapping where objects should be merged: a real number
% between 0.0 and 1.0 which controls what degree of overlapping between 2
% objects implicates the merging (union) of those objects.
%
% TECHNICAL DESCRIPTION OF THE PROPAGATION OPTION:

% See also Identify primary, MeasureObjectAreaShape, MeasureObjectIntensity and
% Relate modules.

% CellProfiler is distributed under the GNU General Public License.
% See the accompanying file LICENSE for details.
%
% Developed by Krisztian Koos.
% Copyright 2015.
%
% Please see the AUTHORS file for credits.
%
% Website: http://www.cellprofiler.org
%
% $Revision: 5025 $

%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%
%drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%%% Sets up loop for test mode.
% if strcmp(char(handles.Settings.VariableValues{CurrentModuleNum,12}),'Yes') %%% TestMode
%     IdentChoiceList = {'Distance - N' 'Distance - B' 'Watershed' 'Propagation'};
% else
%     IdentChoiceList = {char(handles.Settings.VariableValues{CurrentModuleNum,3})};
% end

%textVAR01 = What did you call the primary objects you want to create secondary objects around?
%infotypeVAR01 = objectgroup
PrimaryObjectName = char(handles.Settings.VariableValues{CurrentModuleNum,1});
%inputtypeVAR01 = popupmenu

%textVAR02 = What do you want to call the objects identified by this module?
%defaultVAR02 = Nucleoli
%infotypeVAR02 = objectgroup indep
SecondaryObjectName = char(handles.Settings.VariableValues{CurrentModuleNum,2});

%textVAR03 = Select the method to identify the secondary objects (Distance - B uses background; Distance - N does not):
%choiceVAR03 = Hill-descending
%inputtypeVAR03 = popupmenu
OriginalIdentChoice = char(handles.Settings.VariableValues{CurrentModuleNum,3});

%textVAR04 = What did you call the images to be used to find the edges of the secondary objects? For DISTANCE - N, this will not affect object identification, only the final display.
%infotypeVAR04 = imagegroup
ImageName = char(handles.Settings.VariableValues{CurrentModuleNum,4});
%inputtypeVAR04 = popupmenu

%textVAR05 = What do you want to call the outlines of the identified nucleoli (optional)?
%defaultVAR05 = Do not save
%infotypeVAR05 = outlinegroup indep
SaveOutlines = char(handles.Settings.VariableValues{CurrentModuleNum,5});

%textVAR06 = Ratio of overlapping where objects should be merged
%defaultVAR06 = 0.95
MaxOverlapRatio = str2double(char(handles.Settings.VariableValues{CurrentModuleNum,6}));

%%%VariableRevisionNumber = 3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRELIMINARY CALCULATIONS & FILE HANDLING %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%drawnow

%%% Reads (opens) the image you want to analyze and assigns it to a
%%% variable.
OrigImage = CPretrieveimage(handles,ImageName,ModuleName,'MustBeGray','CheckScale');

%%% Retrieves the preliminary label matrix image that contains the primary
%%% segmented objects which have only been edited to discard objects
%%% that are smaller than a certain size.  This image
%%% will be used as markers to segment the secondary objects with this
%%% module.  Checks first to see whether the appropriate image exists.
PrelimPrimaryLabelMatrixImage = CPretrieveimage(handles,['SmallRemovedSegmented', PrimaryObjectName],ModuleName,'DontCheckColor','DontCheckScale',size(OrigImage));

%%% Retrieves the label matrix image that contains the edited primary
%%% segmented objects which will be used to weed out which objects are
%%% real - not on the edges and not below or above the specified size
%%% limits. Checks first to see whether the appropriate image exists.
EditedPrimaryLabelMatrixImage = CPretrieveimage(handles,['Segmented', PrimaryObjectName],ModuleName,'DontCheckColor','DontCheckScale',size(OrigImage));

%%% Converts the EditedPrimaryBinaryImage to binary.
EditedPrimaryBinaryImage = im2bw(EditedPrimaryLabelMatrixImage,.5);

%%%%%%%%%%%%%%%%%%%%%%
%%% IMAGE ANALYSIS %%%
%%%%%%%%%%%%%%%%%%%%%%
%drawnow

LogicalOutlines = false(size(OrigImage));

[L, num] = bwlabel(EditedPrimaryBinaryImage);
FinalLabelMatrixImage = false( size(EditedPrimaryBinaryImage) );
CentroidsImage = false( size(EditedPrimaryBinaryImage) );
CentroidsImage2 = uint16(zeros( size(EditedPrimaryBinaryImage) ));
ParentList = [];
ChildCounts = zeros(1,num);

numOfNucleoli = 0;

if num>0
    segstats = regionprops(L, 'Area', 'BoundingBox', 'Image', 'Orientation', 'Perimeter', ...
    'Solidity', 'Eccentricity', 'EquivDiameter', 'MinorAxisLength', 'MajorAxisLength', 'PixelIdxList');
    for i=1:num
        bb = segstats(i).BoundingBox;
        neoliBox = OrigImage(ceil(bb(2)):ceil(bb(2))+bb(4)-1, ceil(bb(1)):ceil(bb(1))+bb(3)-1);
        objects{i} = segmentNucleoli(neoliBox, segstats(i).Image, MaxOverlapRatio);
        
        ParentList = [ParentList; zeros(length(objects{i}),1)+i];
        ChildCounts(i) = length(objects{i});
        
        for j=1:length(objects{i})
            LogicalOutlines(ceil(bb(2)):ceil(bb(2))+bb(4)-1, ceil(bb(1)):ceil(bb(1))+bb(3)-1) = ...
                LogicalOutlines(ceil(bb(2)):ceil(bb(2))+bb(4)-1, ceil(bb(1)):ceil(bb(1))+bb(3)-1) | bwperim(objects{i}{j});
            FinalLabelMatrixImage(ceil(bb(2)):ceil(bb(2))+bb(4)-1, ceil(bb(1)):ceil(bb(1))+bb(3)-1) = ...
                FinalLabelMatrixImage(ceil(bb(2)):ceil(bb(2))+bb(4)-1, ceil(bb(1)):ceil(bb(1))+bb(3)-1) | objects{i}{j};
            props = regionprops(objects{i}{j},'Centroid');
            CentroidsImage(int16(bb(2)+props.Centroid(2)), int16(bb(1)+props.Centroid(1))) = 1;
            CentroidsImage2(int16(bb(2)+props.Centroid(2)), int16(bb(1)+props.Centroid(1))) = numOfNucleoli+j;
        end
        numOfNucleoli = numOfNucleoli+length(objects{i});
    end
end

ObjectCount = length(ParentList);

% [L2, num2] = bwlabel(CentroidsImage,4);

%%%%%%%%%%%%%%%%%%%%%%%
%%% DISPLAY RESULTS %%%
%%%%%%%%%%%%%%%%%%%%%%%

LineIntensity = max(OrigImage(:));

ObjectOutlinesOnOrigImage = OrigImage;
ObjectOutlinesOnOrigImage(LogicalOutlines) = LineIntensity;

BothOutlinesOnOrigImage = ObjectOutlinesOnOrigImage;
PrimaryObjectOutlines = bwperim(EditedPrimaryBinaryImage);
BothOutlinesOnOrigImage(PrimaryObjectOutlines) = LineIntensity;

ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);
if any(findobj == ThisModuleFigureNumber)
    %%% Activates the appropriate figure window.
    CPfigure(handles,'Image',ThisModuleFigureNumber);
    if handles.Current.SetBeingAnalyzed == handles.Current.StartingImageSet
        CPresizefigure(OrigImage,'TwoByTwo',ThisModuleFigureNumber);
    end
%     ObjectCoverage = 100*sum(sum(FinalLabelMatrixImage > 0))/numel(FinalLabelMatrixImage);
%     uicontrol(ThisModuleFigureNumber,'Style','Text','Units','Normalized','Position',[0.25 0.01 .6 0.04],...
%         'BackgroundColor',[.7 .7 .9],'HorizontalAlignment','Left','String',sprintf('Threshold:  %0.3f               %0.1f%% of image consists of objects',Threshold,ObjectCoverage),'FontSize',handles.Preferences.FontSize);
    %%% A subplot of the figure window is set to display the original image.
    subplot(2,2,1);
    CPimagesc(OrigImage,handles);
    title(['Input Image, cycle # ',num2str(handles.Current.SetBeingAnalyzed)]);
    %%% A subplot of the figure window is set to display the colored label
    %%% matrix image.
    subplot(2,2,2);
    %% Ray hp
    
    ColoredLabelMatrixImage = CPlabel2rgb(handles,FinalLabelMatrixImage);
    CPimagesc(ColoredLabelMatrixImage,handles);
    title(['Outlined ',SecondaryObjectName]);
    %%% A subplot of the figure window is set to display the original image
    %%% with secondary object outlines drawn on top.
    subplot(2,2,3);
    CPimagesc(ObjectOutlinesOnOrigImage,handles);
    title([SecondaryObjectName, ' Outlines on Input Image']);
    %%% A subplot of the figure window is set to display the original
    %%% image with outlines drawn for both the primary and secondary
    %%% objects.
    subplot(2,2,4);
    tmp = OrigImage/max(OrigImage(:));
    OutlinedObjectsR = tmp;
    OutlinedObjectsR(LogicalOutlines) = max(tmp(:));
    OutlinedObjectsR(PrimaryObjectOutlines) = min(tmp(:));
    OutlinedObjectsG = tmp;
    OutlinedObjectsG(LogicalOutlines) = min(tmp(:));
    OutlinedObjectsG(PrimaryObjectOutlines) = max(tmp(:));
    OutlinedObjectsB = tmp;
    OutlinedObjectsB(LogicalOutlines) = min(tmp(:));
    OutlinedObjectsB(PrimaryObjectOutlines) = min(tmp(:));
    OutlinedObjects = cat(3,OutlinedObjectsR,OutlinedObjectsG,OutlinedObjectsB);
    CPimagesc(OutlinedObjects,handles);
%     CPimagesc(BothOutlinesOnOrigImage,handles);
    title(['Outlines of ', PrimaryObjectName, ' and ', SecondaryObjectName, ' on Input Image']);
    drawnow
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     CALCULATE MEASUREMENTS     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Calculate and store basic features from MeasureObjectAreaShape module

%%% Retrieves the pixel size that the user entered (micrometers per pixel).
PixelSize = str2double(handles.Settings.PixelSize);

BasicObjectAreaShapeFeatures = {'Area',...
    'Eccentricity',...
    'Solidity',...
    'Extent',...
    'EulerNumber',...
    'Perimeter',...
    'FormFactor',...
    'MajorAxisLength',...
    'MinorAxisLength'...
    'Orientation'};

BasicObjectAreaShape = [];

for i=1:num
    for j=1:length(objects{i})
        props = regionprops(objects{i}{j},'Area','Eccentricity','Solidity','Extent','EulerNumber',...
            'MajorAxisLength','MinorAxisLength','Perimeter','Orientation');
        FormFactor = (4*pi*cat(1,props.Area)) ./ ((cat(1,props.Perimeter)+1).^2);       % Add 1 to perimeter to avoid divide by zero
        BasicObjectAreaShape = [BasicObjectAreaShape;...
            cat(1,props.Area)*PixelSize^2,...
            cat(1,props.Eccentricity),...
            cat(1,props.Solidity),...
            cat(1,props.Extent),...
            cat(1,props.EulerNumber),...
            cat(1,props.Perimeter)*PixelSize,...
            FormFactor,...
            cat(1,props.MajorAxisLength)*PixelSize,...
            cat(1,props.MinorAxisLength)*PixelSize,...
            cat(1,props.Orientation)];
        
    end
end

%%% Calculate and store basic features from MeasureObjectIntensity module

%%% Initialize measurement structure
% BasicObjectIntensity = zeros(ObjectCount,11);
BasicObjectIntensity = [];
BasicObjectIntensityFeatures    = {'IntegratedIntensity',...
    'MeanIntensity',...
    'StdIntensity',...
    'MinIntensity',...
    'MaxIntensity',...
    'IntegratedIntensityEdge',...
    'MeanIntensityEdge',...
    'StdIntensityEdge',...
    'MinIntensityEdge',...
    'MaxIntensityEdge',...
    'MassDisplacement'};

for i=1:num
    for j=1:length(objects{i})
        object = objects{i}{j};
        bb = segstats(i).BoundingBox;
        neoliBox = OrigImage(ceil(bb(2)):ceil(bb(2))+bb(4)-1, ceil(bb(1)):ceil(bb(1))+bb(3)-1);
        props = regionprops(object,'PixelIdxList');
        if isempty(props.PixelIdxList)
            BasicObjectIntensity = [BasicObjectIntensity; zeros(1,11)];
            continue;
        end
        
        %%% Measure basic set of Intensity features
        BasicObjectIntensityVector = [ sum(neoliBox(props.PixelIdxList)),...
            mean(neoliBox(props.PixelIdxList)),...
            std(neoliBox(props.PixelIdxList)),...
            min(neoliBox(props.PixelIdxList)),...
            max(neoliBox(props.PixelIdxList)) ];
        
        [sr, sc] = size(object);
        [r,c] = ind2sub([sr sc],props.PixelIdxList);
        rmax = min(sr,max(r));
        rmin = max(1,min(r));
        cmax = min(sc,max(c));
        cmin = max(1,min(c));
        BWim = object;
        Greyim = neoliBox(rmin:rmax,cmin:cmax);
        LabelBoundaryImage = CPlabelperim(object);
        Boundaryim = LabelBoundaryImage(rmin:rmax,cmin:cmax) == 1;
        perim = Greyim(Boundaryim);
        
        BasicObjectIntensityVector = [BasicObjectIntensityVector,...
            sum(perim),...
            mean(perim),...
            std(perim),...
            min(perim),...
            max(perim)];
        %%% Calculate the Mass displacment (taking the pixelsize into account), which is the distance between
        %%% the center of gravity in the gray level image and the binary
        %%% image.
        
        BWx = sum((1:size(BWim,2)).*sum(BWim,1))/sum(1:size(BWim,2));
        BWy = sum((1:size(BWim,1))'.*sum(BWim,2))/sum(1:size(BWim,1));
        Greyx = sum((1:size(Greyim,2)).*sum(Greyim,1))/sum(1:size(Greyim,2));
        Greyy = sum((1:size(Greyim,1))'.*sum(Greyim,2))/sum(1:size(Greyim,1));
        BasicObjectIntensityVector = [BasicObjectIntensityVector, sqrt((BWx-Greyx)^2+(BWy-Greyy)^2)*PixelSize];
        
        BasicObjectIntensity = [BasicObjectIntensity; BasicObjectIntensityVector];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SAVE DATA TO HANDLES STRUCTURE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%drawnow

fieldname = ['Segmented',SecondaryObjectName];
handles.Pipeline.(fieldname) = FinalLabelMatrixImage;

fieldname = ['SmallRemovedSegmented',SecondaryObjectName];
handles.Pipeline.(fieldname) = FinalLabelMatrixImage;

handles = CPsaveObjectCount(handles, SecondaryObjectName, CentroidsImage2);
handles = CPsaveObjectLocations(handles, SecondaryObjectName, CentroidsImage2);

handles = CPaddmeasurements(handles,SecondaryObjectName,'Parent',PrimaryObjectName,ParentList);
handles = CPaddmeasurements(handles,PrimaryObjectName,'Children',[SecondaryObjectName,'Count'],ChildCounts);

handles.Measurements.(SecondaryObjectName).AreaShapeFeatures = cat(2,BasicObjectAreaShapeFeatures);
handles.Measurements.(SecondaryObjectName).AreaShape{handles.Current.SetBeingAnalyzed} = BasicObjectAreaShape;

handles.Measurements.(SecondaryObjectName).(['Intensity_',ImageName,'Features']) = BasicObjectIntensityFeatures;
handles.Measurements.(SecondaryObjectName).(['Intensity_',ImageName]){handles.Current.SetBeingAnalyzed} = BasicObjectIntensity;

try
    if ~strcmpi(SaveOutlines,'Do not save')
        handles.Pipeline.(SaveOutlines) = LogicalOutlines;
    end
catch
    error(['The object outlines were not calculated by the ', ModuleName, ' module, so these images were not saved to the handles structure. The Save Images module will therefore not function on these images. This is just for your information - image processing is still in progress, but the Save Images module will fail if you attempted to save these images.'])
end

%{
for IdentChoiceNumber = 1:length(IdentChoiceList)

    IdentChoice = IdentChoiceList{IdentChoiceNumber};

    if strncmp(IdentChoice,'Distance',8)
        if strcmp(IdentChoice(12),'N')
            %%% Creates the structuring element using the user-specified size.
%             StructuringElement = strel('disk', DistanceToDilate);
%             %%% Dilates the preliminary label matrix image (edited for small only).
%             DilatedPrelimSecObjectLabelMatrixImage = imdilate(PrelimPrimaryLabelMatrixImage, StructuringElement);
%             %%% Converts to binary.
%             DilatedPrelimSecObjectBinaryImage = im2bw(DilatedPrelimSecObjectLabelMatrixImage,.5);
%             %%% Computes nearest neighbor image of nuclei centers so that the dividing
%             %%% line between secondary objects is halfway between them rather than
%             %%% favoring the primary object with the greater label number.

            [dist, Labels] = bwdist( full( PrelimPrimaryLabelMatrixImage>0) ); %#ok We want to ignore MLint error checking for this line.
            DilatedPrelimSecObjectBinaryImage = dist < DistanceToDilate;            
            %drawnow
            %%% Remaps labels in Labels to labels in PrelimPrimaryLabelMatrixImage.
            if max(Labels(:)) == 0,
                Labels = ones(size(Labels));
            end
            ExpandedRelabeledDilatedPrelimSecObjectImage = PrelimPrimaryLabelMatrixImage(Labels);
            %%% Removes the background pixels (those not labeled as foreground in the
            %%% DilatedPrelimSecObjectBinaryImage). This is necessary because the
            %%% nearest neighbor function assigns *every* pixel to a nucleus, not just
            %%% the pixels that are part of a secondary object.
            RelabeledDilatedPrelimSecObjectImage = zeros(size(ExpandedRelabeledDilatedPrelimSecObjectImage));
            RelabeledDilatedPrelimSecObjectImage(DilatedPrelimSecObjectBinaryImage) = ExpandedRelabeledDilatedPrelimSecObjectImage(DilatedPrelimSecObjectBinaryImage);
            %drawnow
        elseif strcmp(IdentChoice(12),'B')
            [labels_out,d]=IdentifySecPropagateSubfunction(PrelimPrimaryLabelMatrixImage,OrigImage,ThresholdedOrigImage,1.0);
            labels_out(d>DistanceToDilate) = 0;
            labels_out((PrelimPrimaryLabelMatrixImage > 0)) = PrelimPrimaryLabelMatrixImage((PrelimPrimaryLabelMatrixImage > 0));
            RelabeledDilatedPrelimSecObjectImage = labels_out;
        end


        %%% Removes objects that are not in the edited EditedPrimaryLabelMatrixImage.
                
        Map = sparse(1:prod(size(PrelimPrimaryLabelMatrixImage)), PrelimPrimaryLabelMatrixImage(:)+1, EditedPrimaryLabelMatrixImage(:));
        LookUpColumn = full(max(Map,[], 1));
        LookUpColumn(1)=0;
        FinalLabelMatrixImage = LookUpColumn(RelabeledDilatedPrelimSecObjectImage+1);

    elseif strcmp(IdentChoice,'Propagation')
        %%% STEP 2: Starting from the identified primary objects, the secondary
        %%% objects are identified using the propagate function, written by Thouis
        %%% R. Jones. Calls the function
        %%% "IdentifySecPropagateSubfunction.mexmac" (or whichever version is
        %%% appropriate for the computer platform being used), which consists of C
        %%% code that has been compiled to run quickly within Matlab.
        
        % 2007-Jul-16 Kyungnam: If you want to get additional outputs, then
        % add more output arguments as follows:
        %%% [PropagatedImage, dist, diff_count, pop_count] = IdentifySecPropagateSubfunction(PrelimPrimaryLabelMatrixImage,OrigImage,ThresholdedOrigImage,RegularizationFactor);
        PropagatedImage = IdentifySecPropagateSubfunction(PrelimPrimaryLabelMatrixImage,OrigImage,ThresholdedOrigImage,RegularizationFactor);
        %drawnow

        %%% STEP 3: We used the PrelimPrimaryLabelMatrixImage as the
        %%% source for primary objects, but that label-matrix is built
        %%% before small/large objects and objects touching the
        %%% boundary are removed.  We need to filter the label matrix
        %%% from propagate to make the labels match, and remove any
        %%% secondary objects that correspnd to size- or
        %%% boundary-filtered primaries.
        %%%
        %%% Map preliminary labels to edited labels based on maximum
        %%% overlap from prelim to edited.  We can probably assume
        %%% that no borders are adjusted during editing (i.e., any
        %%% changes from Prelim to Edited only involves removing
        %%% entire objects), but this is safer.
        %%% 
        %%% (add one so that zeros are remapped correctly.)
        PrelimToEditedHist = sparse(EditedPrimaryLabelMatrixImage(:) + 1, PrelimPrimaryLabelMatrixImage(:) + 1, 1);
        [ignore, PrelimToEditedRemap] = sort(PrelimToEditedHist, 1);
        PrelimToEditedRemap = PrelimToEditedRemap(end, :) - 1;
        %%% make sure zeros map to zeros (note the off-by-one for the
        %%% index because Matlab doesn't do 0-indexing).
        PrelimToEditedRemap(1) = 0;
        EditedLabelMatrixImage = PrelimToEditedRemap(PropagatedImage + 1);

        %%% STEP 4:
        %%%
        %%% Fill holes (any contiguous, all-0 regions that are
        %%% surrounded by a single value).
        FinalLabelMatrixImage = CPfill_holes(EditedLabelMatrixImage);
        
    elseif strcmp(IdentChoice,'Watershed')
        %%% In order to use the watershed transform to find dividing lines between
        %%% the secondary objects, it is necessary to identify the foreground
        %%% objects and to identify a portion of the background.  The foreground
        %%% objects are retrieved as the binary image of primary objects from the
        %%% previously run image analysis module.   This forces the secondary
        %%% object's outline to extend at least as far as the edge of the primary
        %%% objects.

        %%% Inverts the image.
        InvertedThresholdedOrigImage = imcomplement(ThresholdedOrigImage);

        %%% NOTE: There are two other ways to mark the background prior to
        %%% watershedding; I think the method used above is best, but I have
        %%% included the ideas for two alternate methods.
        %%% METHOD (2): Threshold the original image (or a smoothed image)
        %%% so that background pixels are black.  This is overly strong, so instead
        %%% of weakly thresholding the image as is done in METHOD (1),  you can then "thin"
        %%% the background pixels by computing the SKIZ
        %%% (skeleton of influence zones), which is done by watershedding the
        %%% distance transform of the thresholded image.  These watershed lines are
        %%% then superimposed on the marked image that will be watershedded to
        %%% segment the objects.  I think this would not produce results different
        %%% from METHOD 1 (the one used above), since METHOD 1 overlays the
        %%% outlines of the primary objects anyway.
        %%% This method is based on the Mathworks Image Processing Toolbox demo
        %%% "Marker-Controlled Watershed Segmentation".  I found it online; I don't
        %%% think it is in the Matlab Demos that are found through help.  It uses
        %%% an image of a box of oranges.
        %%%
        %%% METHOD (3):  (I think this method does not work well for clustered
        %%% objects.)  The distance transformed image containing the marked objects
        %%% is watershedded, which produces lines midway between the marked
        %%% objects.  These lines are superimposed on the marked image that will be
        %%% watershedded to segment the objects. But if marked objects are
        %%% clustered and not a uniform distance from each other, this will produce
        %%% background lines on top of actual objects.
        %%% This method is based on Gonzalez, et al. Digital Image Processing using
        %%% Matlab, page 422-425.

        %%% STEP 2: Identify the outlines of each primary object, so that each
        %%% primary object can be definitely separated from the background.  This
        %%% solves the problem of some primary objects running
        %%% right up against the background pixels and therefore getting skipped.
        %%% Note: it is less accurate and less fast to use edge detection (sobel)
        %%% to identify the edges of the primary objects.
        %drawnow
        %%% Converts the PrelimPrimaryLabelMatrixImage to binary.
        PrelimPrimaryBinaryImage = im2bw(PrelimPrimaryLabelMatrixImage,.5);
        %%% Creates the structuring element that will be used for dilation.
        StructuringElement = strel('square',3);
        %%% Dilates the Primary Binary Image by one pixel (8 neighborhood).
        DilatedPrimaryBinaryImage = imdilate(PrelimPrimaryBinaryImage, StructuringElement);
        %%% Subtracts the PrelimPrimaryBinaryImage from the DilatedPrimaryBinaryImage,
        %%% which leaves the PrimaryObjectOutlines.
        PrimaryObjectOutlines = DilatedPrimaryBinaryImage - PrelimPrimaryBinaryImage;

        %%% STEP 3: Produce the marker image which will be used for the first
        %%% watershed.
        %drawnow
        %%% Combines the foreground markers and the background markers.
        BinaryMarkerImagePre = PrelimPrimaryBinaryImage | InvertedThresholdedOrigImage;
        %%% Overlays the PrimaryObjectOutlines to maintain distinctions between each
        %%% primary object and the background.
        BinaryMarkerImage = BinaryMarkerImagePre;
        BinaryMarkerImage(PrimaryObjectOutlines == 1) = 0;

        %%% STEP 4: Calculate the Sobel image, which reflects gradients, which will
        %%% be used for the watershedding function.
        %drawnow
        %%% Calculates the 2 sobel filters.  The sobel filter is directional, so it
        %%% is used in both the horizontal & vertical directions and then the
        %%% results are combined.
        filter1 = fspecial('sobel');
        filter2 = filter1';
        %%% Applies each of the sobel filters to the original image.
        I1 = imfilter(OrigImage, filter1);
        I2 = imfilter(OrigImage, filter2);
        %%% Adds the two images.
        %%% The Sobel operator results in negative values, so the absolute values
        %%% are calculated to prevent errors in future steps.
        AbsSobeledImage = abs(I1) + abs(I2);

        %%% STEP 5: Perform the first watershed.
        %drawnow

        %%% Overlays the foreground and background markers onto the
        %%% absolute value of the Sobel Image, so there are black nuclei on top of
        %%% each dark object, with black background.
        Overlaid = imimposemin(AbsSobeledImage, BinaryMarkerImage);
        %%% Perform the watershed on the marked absolute-value Sobel Image.
        BlackWatershedLinesPre = watershed(Overlaid);
        %%% Bug workaround (see step 9).
        BlackWatershedLinesPre2 = im2bw(BlackWatershedLinesPre,.5);
        BlackWatershedLines = bwlabel(BlackWatershedLinesPre2);

        %%% STEP 6: Identify and extract the secondary objects, using the watershed
        %%% lines.
        %drawnow
        %%% The BlackWatershedLines image is a label matrix where the watershed
        %%% lines = 0 and each distinct object is assigned a number starting at 1.
        %%% This image is converted to a binary image where all the objects = 1.
        SecondaryObjects1 = im2bw(BlackWatershedLines,.5);
        %%% Identifies objects in the binary image using bwlabel.
        %%% Note: Matlab suggests that in some circumstances bwlabeln is faster
        %%% than bwlabel, even for 2D images.  I found that in this case it is
        %%% about 10 times slower.
        LabelMatrixImage1 = bwlabel(SecondaryObjects1,4);
        %drawnow

        %%% STEP 7: Discarding background "objects".  The first watershed function
        %%% simply divides up the image into regions.  Most of these regions
        %%% correspond to actual objects, but there are big blocks of background
        %%% that are recognized as objects. These can be distinguished from actual
        %%% objects because they do not overlap a primary object.

        %%% The following changes all the labels in LabelMatrixImage1 to match the
        %%% centers they enclose (from PrelimPrimaryBinaryImage), and marks as background
        %%% any labeled regions that don't overlap a center. This function assumes
        %%% that every center is entirely contained in one labeled area.  The
        %%% results if otherwise may not be well-defined. The non-background labels
        %%% will be renumbered according to the center they enclose.

        %%% Finds the locations and labels for different regions.
        area_locations = find(LabelMatrixImage1);
        area_labels = LabelMatrixImage1(area_locations);
        %%% Creates a sparse matrix with column as label and row as location,
        %%% with the value of the center at (I,J) if location I has label J.
        %%% Taking the maximum of this matrix gives the largest valued center
        %%% overlapping a particular label.  Tacking on a zero and pushing
        %%% labels through the resulting map removes any background regions.
        map = [0 full(max(sparse(area_locations, area_labels, PrelimPrimaryBinaryImage(area_locations))))];
        ActualObjectsBinaryImage = map(LabelMatrixImage1 + 1);

        %%% STEP 8: Produce the marker image which will be used for the second
        %%% watershed.
        %drawnow
        %%% The module has now produced a binary image of actual secondary
        %%% objects.  The gradient (Sobel) image was used for watershedding, which
        %%% produces very nice divisions between objects that are clumped, but it
        %%% is too stringent at the edges of objects that are isolated, and at the
        %%% edges of clumps of objects. Therefore, the stringently identified
        %%% secondary objects are used as markers for a second round of
        %%% watershedding, this time based on the original (intensity) image rather
        %%% than the gradient image.

        %%% Creates the structuring element that will be used for dilation.
        StructuringElement = strel('square',3);
        %%% Dilates the Primary Binary Image by one pixel (8 neighborhood).
        DilatedActualObjectsBinaryImage = imdilate(ActualObjectsBinaryImage, StructuringElement);
        %%% Subtracts the PrelimPrimaryBinaryImage from the DilatedPrimaryBinaryImage,
        %%% which leaves the PrimaryObjectOutlines.
        ActualObjectOutlines = DilatedActualObjectsBinaryImage - ActualObjectsBinaryImage;
        %%% Produces the marker image which will be used for the watershed. The
        %%% foreground markers are taken from the ActualObjectsBinaryImage; the
        %%% background markers are taken from the same image as used in the first
        %%% round of watershedding: InvertedThresholdedOrigImage.
        BinaryMarkerImagePre2 = ActualObjectsBinaryImage | InvertedThresholdedOrigImage;
        %%% Overlays the ActualObjectOutlines to maintain distinctions between each
        %%% secondary object and the background.
        BinaryMarkerImage2 = BinaryMarkerImagePre2;
        BinaryMarkerImage2(ActualObjectOutlines == 1) = 0;

        %%% STEP 9: Perform the second watershed.
        %%% As described above, the second watershed is performed on the original
        %%% intensity image rather than on a gradient (Sobel) image.
        %drawnow
        %%% Inverts the original image.
        InvertedOrigImage = imcomplement(OrigImage);
        %%% Overlays the foreground and background markers onto the
        %%% InvertedOrigImage, so there are black secondary object markers on top
        %%% of each dark secondary object, with black background.
        MarkedInvertedOrigImage = imimposemin(InvertedOrigImage, BinaryMarkerImage2);
        %%% Performs the watershed on the MarkedInvertedOrigImage.
        SecondWatershedPre = watershed(MarkedInvertedOrigImage);
        %%% BUG WORKAROUND:
        %%% There is a bug in the watershed function of Matlab that often results in
        %%% the label matrix result having two objects labeled with the same label.
        %%% I am not sure whether it is a bug in how the watershed image is
        %%% produced (it seems so: the resulting objects often are nowhere near the
        %%% regional minima) or whether it is simply a problem in the final label
        %%% matrix calculation. Matlab has been informed of this issue and has
        %%% confirmed that it is a bug (February 2004). I think that it is a
        %%% reasonable fix to convert the result of the watershed to binary and
        %%% remake the label matrix so that each label is used only once. In later
        %%% steps, inappropriate regions are weeded out anyway.
        SecondWatershedPre2 = im2bw(SecondWatershedPre,.5);
        SecondWatershed = bwlabel(SecondWatershedPre2);
        %drawnow

        %%% STEP 10: As in step 7, remove objects that are actually background
        %%% objects.  See step 7 for description. This time, the edited primary object image is
        %%% used rather than the preliminary one, so that objects whose nuclei are
        %%% on the edge of the image and who are larger or smaller than the
        %%% specified size are discarded.

        %%% Finds the locations and labels for different regions.
        area_locations2 = find(SecondWatershed);
        area_labels2 = SecondWatershed(area_locations2);
        %%% Creates a sparse matrix with column as label and row as location,
        %%% with the value of the center at (I,J) if location I has label J.
        %%% Taking the maximum of this matrix gives the largest valued center
        %%% overlapping a particular label.  Tacking on a zero and pushing
        %%% labels through the resulting map removes any background regions.
        map2 = [0 full(max(sparse(area_locations2, area_labels2, EditedPrimaryBinaryImage(area_locations2))))];
        FinalBinaryImagePre = map2(SecondWatershed + 1);
        %%% Fills holes in the FinalBinaryPre image.
        FinalBinaryImage = imfill(FinalBinaryImagePre, 'holes');
        %%% Converts the image to label matrix format. Even if the above step
        %%% is excluded (filling holes), it is still necessary to do this in order
        %%% to "compact" the label matrix: this way, each number corresponds to an
        %%% object, with no numbers skipped.
        ActualObjectsLabelMatrixImage3 = bwlabel(FinalBinaryImage);
        %%% The final objects are relabeled so that their numbers
        %%% correspond to the numbers used for nuclei.
        %%% For each object, one label and one label location is acquired and
        %%% stored.
        [LabelsUsed,LabelLocations] = unique(EditedPrimaryLabelMatrixImage);
        %%% The +1 increment accounts for the fact that there are zeros in the
        %%% image, while the LabelsUsed starts at 1.
        LabelsUsed(ActualObjectsLabelMatrixImage3(LabelLocations(2:end))+1) = EditedPrimaryLabelMatrixImage(LabelLocations(2:end));
        FinalLabelMatrixImagePre = LabelsUsed(ActualObjectsLabelMatrixImage3+1);
        %%% The following is a workaround for what seems to be a bug in the
        %%% watershed function: very very rarely two nuclei end up sharing one
        %%% "cell" object, so that one of the nuclei ends up without a
        %%% corresponding cell.  I am trying to determine why this happens exactly.
        %%% When the cell is measured, the area (and other
        %%% measurements) are recorded as [], which causes problems when dependent
        %%% measurements (e.g. perimeter/area) are attempted.  It results in divide
        %%% by zero errors and the mean area = NaN and so on.  So, the Primary
        %%% label matrix image (where it is nonzero) is written onto the Final cell
        %%% label matrix image pre so that every primary object has at least some
        %%% pixels of secondary object.
        FinalLabelMatrixImage = FinalLabelMatrixImagePre;
        FinalLabelMatrixImage(EditedPrimaryLabelMatrixImage ~= 0) = EditedPrimaryLabelMatrixImage(EditedPrimaryLabelMatrixImage ~= 0);
    end

    %%% Calculates the ColoredLabelMatrixImage for displaying in the figure
    %%% window in subplot(2,2,2).
 %   ColoredLabelMatrixImage = CPlabel2rgb(handles,FinalLabelMatrixImage);
    %%% Calculates OutlinesOnOrigImage for displaying in the figure
    %%% window in subplot(2,2,3).
    %%% Note: these outlines are not perfectly accurate; for some reason it
    %%% produces more objects than in the original image.  But it is OK for
    %%% display purposes.
    %%% Maximum filters the image with a 3x3 neighborhood.
    MaxFilteredImage = ordfilt2(FinalLabelMatrixImage,9,ones(3,3),'symmetric');
    %%% Determines the outlines.
    IntensityOutlines = FinalLabelMatrixImage - MaxFilteredImage;
    %%% Converts to logical.
    warning off MATLAB:conversionToLogical
    LogicalOutlines = logical(IntensityOutlines);
    warning on MATLAB:conversionToLogical
    %%% Determines the grayscale intensity to use for the cell outlines.
    LineIntensity = max(OrigImage(:));
    %%% Overlays the outlines on the original image.
    ObjectOutlinesOnOrigImage = OrigImage;
    ObjectOutlinesOnOrigImage(LogicalOutlines) = LineIntensity;
    %%% Calculates BothOutlinesOnOrigImage for displaying in the figure
    %%% window in subplot(2,2,4).
    %%% Creates the structuring element that will be used for dilation.
    StructuringElement = strel('square',3);
    %%% Dilates the Primary Binary Image by one pixel (8 neighborhood).
    DilatedPrimaryBinaryImage = imdilate(EditedPrimaryBinaryImage, StructuringElement);
    %%% Subtracts the PrelimPrimaryBinaryImage from the DilatedPrimaryBinaryImage,
    %%% which leaves the PrimaryObjectOutlines.
    PrimaryObjectOutlines = DilatedPrimaryBinaryImage - EditedPrimaryBinaryImage;
    BothOutlinesOnOrigImage = ObjectOutlinesOnOrigImage;
    BothOutlinesOnOrigImage(PrimaryObjectOutlines == 1) = LineIntensity;
    
    if strcmp(OriginalIdentChoice,IdentChoice)
        if ~isfield(handles.Measurements,SecondaryObjectName)
            handles.Measurements.(SecondaryObjectName) = {};
        end

        if ~isfield(handles.Measurements,PrimaryObjectName)
            handles.Measurements.(PrimaryObjectName) = {};
        end

        handles = CPrelateobjects(handles,SecondaryObjectName,PrimaryObjectName,FinalLabelMatrixImage,EditedPrimaryLabelMatrixImage,ModuleName);

        %%%%%%%%%%%%%%%%%%%%%%%
        %%% DISPLAY RESULTS %%%
        %%%%%%%%%%%%%%%%%%%%%%%


        ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);
        if any(findobj == ThisModuleFigureNumber)
                        %%% Activates the appropriate figure window.
            CPfigure(handles,'Image',ThisModuleFigureNumber);
            if handles.Current.SetBeingAnalyzed == handles.Current.StartingImageSet
                CPresizefigure(OrigImage,'TwoByTwo',ThisModuleFigureNumber);
            end
            ObjectCoverage = 100*sum(sum(FinalLabelMatrixImage > 0))/numel(FinalLabelMatrixImage);
            uicontrol(ThisModuleFigureNumber,'Style','Text','Units','Normalized','Position',[0.25 0.01 .6 0.04],...
                'BackgroundColor',[.7 .7 .9],'HorizontalAlignment','Left','String',sprintf('Threshold:  %0.3f               %0.1f%% of image consists of objects',Threshold,ObjectCoverage),'FontSize',handles.Preferences.FontSize);
            %%% A subplot of the figure window is set to display the original image.
            subplot(2,2,1);
            CPimagesc(OrigImage,handles);
            title(['Input Image, cycle # ',num2str(handles.Current.SetBeingAnalyzed)]);
            %%% A subplot of the figure window is set to display the colored label
            %%% matrix image.
            subplot(2,2,2);
            %% Ray hp
            ColoredLabelMatrixImage = CPlabel2rgb(handles,FinalLabelMatrixImage);                
            CPimagesc(ColoredLabelMatrixImage,handles);
            title(['Outlined ',SecondaryObjectName]);
            %%% A subplot of the figure window is set to display the original image
            %%% with secondary object outlines drawn on top.
            subplot(2,2,3);
            CPimagesc(ObjectOutlinesOnOrigImage,handles);
            title([SecondaryObjectName, ' Outlines on Input Image']);
            %%% A subplot of the figure window is set to display the original
            %%% image with outlines drawn for both the primary and secondary
            %%% objects.
            subplot(2,2,4);
            CPimagesc(BothOutlinesOnOrigImage,handles);
            title(['Outlines of ', PrimaryObjectName, ' and ', SecondaryObjectName, ' on Input Image']);
            drawnow
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% SAVE DATA TO HANDLES STRUCTURE %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %drawnow

        %%% Saves the final, segmented label matrix image of secondary objects to
        %%% the handles structure so it can be used by subsequent modules.
        fieldname = ['Segmented',SecondaryObjectName];
        handles.Pipeline.(fieldname) = FinalLabelMatrixImage;

        if strcmp(IdentChoice,'Propagation')
            %%% Saves the Threshold value to the handles structure.
            %%% Storing the threshold is a little more complicated than storing other measurements
            %%% because several different modules will write to the handles.Measurements.Image.Threshold
            %%% structure, and we should therefore probably append the current threshold to an existing structure.
            % First, if the Threshold fields don't exist, initialize them
            if ~isfield(handles.Measurements.Image,'ThresholdFeatures')
                handles.Measurements.Image.ThresholdFeatures = {};
                handles.Measurements.Image.Threshold = {};
            end
            %%% Search the ThresholdFeatures to find the column for this object type
            column = find(~cellfun('isempty',strfind(handles.Measurements.Image.ThresholdFeatures,SecondaryObjectName)));
            %%% If column is empty it means that this particular object has not been segmented before. This will
            %%% typically happen for the first cycle. Append the feature name in the
            %%% handles.Measurements.Image.ThresholdFeatures matrix
            if isempty(column)
                handles.Measurements.Image.ThresholdFeatures(end+1) = {SecondaryObjectName};
                column = length(handles.Measurements.Image.ThresholdFeatures);
            end
            handles.Measurements.Image.Threshold{handles.Current.SetBeingAnalyzed}(1,column) = Threshold;

            %%% Also add the thresholding quality metrics to the measurements
            if exist('WeightedVariance', 'var')
                FeatureName = [SecondaryObjectName '_WeightedVariance'];
                column = find(~cellfun('isempty',strfind(handles.Measurements.Image.ThresholdFeatures,FeatureName)));
                if isempty(column),
                    handles.Measurements.Image.ThresholdFeatures(end+1) = {FeatureName};
                    column = length(handles.Measurements.Image.ThresholdFeatures);
                end
                handles.Measurements.Image.Threshold{handles.Current.SetBeingAnalyzed}(1,column) = WeightedVariance;

                FeatureName = [SecondaryObjectName '_SumOfEntropies'];
                column = find(~cellfun('isempty',strfind(handles.Measurements.Image.ThresholdFeatures,FeatureName)));
                if isempty(column),
                    handles.Measurements.Image.ThresholdFeatures(end+1) = {FeatureName};
                    column = length(handles.Measurements.Image.ThresholdFeatures);
                end
                handles.Measurements.Image.Threshold{handles.Current.SetBeingAnalyzed}(1,column) = SumOfEntropies;
            end
        end

	handles = CPsaveObjectCount(handles, SecondaryObjectName, FinalLabelMatrixImage);
	handles = CPsaveObjectLocations(handles, SecondaryObjectName, FinalLabelMatrixImage);
	
        %%% Saves images to the handles structure so they can be saved to the hard
        %%% drive, if the user requested.
        try
            if ~strcmpi(SaveOutlines,'Do not save')
                handles.Pipeline.(SaveOutlines) = LogicalOutlines;
            end
        catch
            error(['The object outlines were not calculated by the ', ModuleName, ' module, so these images were not saved to the handles structure. The Save Images module will therefore not function on these images. This is just for your information - image processing is still in progress, but the Save Images module will fail if you attempted to save these images.'])
        end
    end
end

%}
end


%%%%%%%%%%%%%%%%%%%
%%% SUBFUNCTION %%%
%%%%%%%%%%%%%%%%%%%

function objects = segmentNucleoli( image, mask, varargin )
%SEGMENTNUCLEOLI Segments nucleoli in one nucleus

%% parse inputs
p = inputParser;
defaultMaxOverlapRatio = 0.95;

addOptional(p, 'maxOverlapRatio', defaultMaxOverlapRatio, @isnumeric);
parse(p, varargin{:});

maxOverlapRatio = p.Results.maxOverlapRatio;
%% init
avgimage = imfilter(image, fspecial('average'), 'replicate');
varimage = stdfilt(avgimage, ones(3));
stdev = std(avgimage(:));
meanVal = mean(avgimage(:));
se = ones(3);

%% Custom, 'hill-descending' algorithm
% Input: pre-segmented image (twomaskImg, optionally mask6) and segmented
% local maxima
% 1. Let the local maxima points be seed (and marked) points.
% 2. Let all other points be unmarked
% 3. If an unmarked pixel has a marked neighbour with greater or equal
%    intensity, then mark this pixel.
% 4. Repeat 3 until there are no new marked pixels.

area1 = mat2gray(varimage)>graythresh(mat2gray(varimage));
area1 = imfill(area1, 'holes');
area2 = mat2gray(avgimage)>=graythresh(mat2gray(avgimage));
maxarea = (area1 | area2) & mask;

localmax = imregionalmax(avgimage);
markStrict = avgimage.*localmax; % remove the probably unnecessary local maxima
markStrict(markStrict<meanVal+stdev) = 0;
markStrict(markStrict>0) = 1;
markPermissive = imdilate(markStrict, se) & imregionalmax(avgimage);
seedPoints = markStrict | markPermissive;
seedPoints = clearMarked2Borders(seedPoints);
objects = findMniPath(avgimage, seedPoints, maxarea);

%% merge objects that overlap at least maxOverlapPerc percent
numObj = numel(objects);
validIndices = true(numObj,1);
for i = 1:numObj-1
    for j = i+1:numObj
        roi = objects{i};
        roi = (roi + objects{j})>1;
        if any(roi(:))
            obj1 = objects{i};
            obj2 = objects{j};
            obj1size = sum(obj1(:));
            obj2size = sum(obj2(:));
            overlapsize = sum(roi(:));
            if overlapsize/obj1size>maxOverlapRatio % the area percentage should be parameter
                validIndices(i) = false;
                objects{j} = objects{j} | objects{i};
            elseif overlapsize/obj2size>maxOverlapRatio
                validIndices(j) = false;
                objects{i} = objects{i} | objects{j};
            end
        end
    end
end
objects = objects(validIndices);

%% Erode one step, because the results seem like to be one pixel bigger
for i = 1:numel(objects)
    objects{i} = imerode(objects{i}, se);
end

%% Dilate in roi until overlap
% There is a problem with multiple overlapping:
% - Which level to assign first
% - what to do when a single pixel line separates the two objects and it
% could be assigned to both
% CURRENT SOLUTION: dilated images are saved to a different variable and
% overlaps are calculated on the original
numObj = numel(objects);
newObjects = cell(numObj,1);
for i = 1:numObj
    newObjects{i} = zeros(size(objects{i}));
end
for i = 1:numObj-1
    for j = i+1:numObj
        roi = (objects{i} + objects{j})>1;
        if any(roi(:))
            overlap = roi;
            baseobj1 = objects{i};
            baseobj2 = objects{j};
            while true
                obj1 = objects{i} - overlap;
                obj2 = objects{j} - overlap;
                dilobj1 = imdilate(obj1,se) & baseobj1 & ~obj2;
                dilobj2 = imdilate(obj2,se) & baseobj2 & ~obj1;
                overlap2 = (dilobj1 + dilobj2)>1;
                newObjects{i} = newObjects{i} + (dilobj1 - overlap2);
                newObjects{j} = newObjects{j} + (dilobj2 - overlap2);
                overlap = overlap2 | (overlap & ~(dilobj1 | dilobj2));
                if isequal(objects{i}.*roi | objects{j}.*roi |  ((obj1 | obj2) & roi), roi)
                    break
                end
            end
        end
    end
end
for i = 1:numObj
    maxval = max(newObjects{i}(:));
    if maxval>0
        objects{i} = newObjects{i} == maxval;
    end
end

%% final morphological operations
objects = customMorphClean(objects);
for i = 1:numel(objects)
    objects{i} = imfill(objects{i}, 'holes');
    objects{i} = bwmorph(objects{i}, 'clean');
end

%%
numObj = numel(objects);
validIndices = true(numObj,1);
for i = 1:numObj
    areaprops = regionprops(objects{i}, 'Area');
    if isempty(areaprops)
        validIndices(i) = false;
        continue
    end
    if length(areaprops)>1
        dataprops = regionprops(objects{i}, 'BoundingBox', 'Image');
        [~, idx] = max([areaprops.Area]);
        objects{i}(:) = 0;
        bbox = round(dataprops(idx).BoundingBox);
        objects{i}(bbox(2):bbox(2)+bbox(4)-1, bbox(1):bbox(1)+bbox(3)-1) = dataprops(idx).Image;
    end
end
objects = objects(validIndices);

end

function A = clearMarked2Borders(A)
[m,n] = size(A);
B = get2BorderPixels(m,n);
A(B) = 0;
end

%% returns an image with 2 pixels tick border
function A = get2BorderPixels(m,n)
A = false(m,n);
A([1:2, end-1:end], :) = true;
A(3:end-2, [1:2, end-1:end]) = true;
end

%% custom morphological cleaning
function objects = customMorphClean(objects)
lut = makeCleanLut();
valid = true(numel(objects),1);
for i = 1:numel(objects)
    A = bwmorph(objects{i},lut);
    while any(A(:))
        objects{i}(A) = ~objects{i}(A);
        A = bwmorph(objects{i},lut);
    end
    if 0 == sum(objects{i}(:))
        valid(i) = false;
    end
end
objects = objects(valid);

end

function objects = findMniPath( img, init_marked, mask )
% FINDMNIPATH Finds monotonic non-increasing path.

%% Custom, 'hill-descending' algorithm
% 3. If an unmarked pixel has a marked neighbour with greater or equal
%    intensity, then mark this pixel.

rop = @le; % SHOULD NOT BE MODIFIED IN THIS FILE (because of the file name and description)
img = padarray(img, [1 1]);
init_marked = padarray(init_marked, [1 1]);
mask = padarray(mask, [1 1]);
[m, n] = size(img);
marked = init_marked;
marked_components = bwconncomp(marked, 8);
final = zeros(m,n);

%% ideas
% 1. mark only if it has at least X marked neighbours greater or equal
% 2. DONE if the dilated marked overlaps an initial marked, then mark that
%    component also
% 3. WONTFIX drop overlapped areas?
% 4. DONE dilate object with overlapped area as roi

%%
needToCheck = true(marked_components.NumObjects, 1);
objects = cell(0);
objCtr = 0;
for i = 1:marked_components.NumObjects
    if ~needToCheck(i)
        continue
    end
    one_marked = false(m, n);
    one_marked(marked_components.PixelIdxList{i}) = true;
    marked = one_marked;
    marked_new = false(m,n);
    while ~isequal(marked, marked_new)
        marked = marked | marked_new;
        marked_new = marked;
%         [marked_new, needToCheck] = findNearbyMarked(img, marked_new, init_marked, marked_components, ...
%             needToCheck);
        candidates = imdilate(marked,ones(3)) - marked;
        candidates = candidates & mask;
        if sum(candidates(:))==0
            break
        end
        marked_elements = find(marked_new);
        for j = 1:numel(marked_elements)
            idx = int32(marked_elements(j));
            val = img(idx);
            if candidates(idx-1) && rop(img(idx-1), val) % north
                marked_new(idx-1) = true;
            end
            if candidates(idx+1) && rop(img(idx+1), val) % south
                marked_new(idx+1) = true;
            end
            if candidates(idx-m) && rop(img(idx-m), val) % west
                marked_new(idx-m) = true;
            end
            if candidates(idx+m) && rop(img(idx+m), val) % east
                marked_new(idx+m) = true;
            end
            if candidates(idx-m-1) && rop(img(idx-m-1), val) % north-west
                marked_new(idx-m-1) = true;
            end
            if candidates(idx-m+1) && rop(img(idx-m+1), val) % south-west
                marked_new(idx-m+1) = true;
            end
            if candidates(idx+m-1) && rop(img(idx+m-1), val) % north-east
                marked_new(idx+m-1) = true;
            end
            if candidates(idx+m+1) && rop(img(idx+m+1), val) % south-east
                marked_new(idx+m+1) = true;
            end
        end
        if sum(marked_new(:))==0
            break;
        end
    end
    marked = bwmorph(marked, 'clean');
    %% add to the final image
    final = final + marked;
    objCtr = objCtr + 1;
    objects{objCtr} = marked(2:end-1, 2:end-1);
end

end

function [marked, needToCheck] = findNearbyMarked(img, marked, init_marked, ...
    marked_components, needToCheck)
dilmarked = imdilate(marked, ones(3));
idxList = find(init_marked & (dilmarked - marked));
for i = 1:numel(idxList)
    pos = idxList(i);
    for j = 1:marked_components.NumObjects
        if any(pos==marked_components.PixelIdxList{j})
            compMeanVal = mean(img(marked_components.PixelIdxList{j}));
            markedMeanVal = mean(img(marked));
            markedSD = std(img(marked));
            
%             markedSeed = marked & init_marked;
            markedSeedMeanVal = mean(find(marked & init_marked));
            
            if compMeanVal>=markedMeanVal+markedSD
%             if compMeanVal>=markedSeedMeanVal
%             if compMeanVal-markedSeedMeanVal<=markedSD
%                 disp(['markedMean=', num2str(markedMeanVal), ', std=', num2str(markedSD),... 
%                     ', intv=[', num2str(markedMeanVal-markedSD), ', ', num2str(markedMeanVal+markedSD), ...
%                     '], compMean=', num2str(compMeanVal)]);
                marked(marked_components.PixelIdxList{j}) = true;
                needToCheck(j) = false;
            end
        end
    end
end
end

function lut = makeCleanLut()
n4N = @(x) isequal(x, [0 0 0; ...
                       0 1 0; ...
                       0 1 0]);
n4S = @(x) isequal(x, [0 1 0; ...
                       0 1 0; ...
                       0 0 0]);
n4W = @(x) isequal(x, [0 0 0; ...
                       0 1 1; ...
                       0 0 0]);
n4E = @(x) isequal(x, [0 0 0; ...
                       1 1 0; ...
                       0 0 0]);
n8NE = @(x) isequal(x, [0 0 0; ...
                        0 1 0; ...
                        1 0 0]);
n8NW = @(x) isequal(x, [0 0 0; ...
                        0 1 0; ...
                        0 0 1]);
n8SW = @(x) isequal(x, [0 0 1; ...
                        0 1 0; ...
                        0 0 0]);
n8SE = @(x) isequal(x, [1 0 0; ...
                        0 1 0; ...
                        0 0 0]);

% spikes
sNE1 = @(x) isequal(x, [0 0 1; ...
                        0 1 1; ...
                        0 0 0]);
sNE2 = @(x) isequal(x, [0 1 1; ...
                        0 1 0; ...
                        0 0 0]);
sNW1 = @(x) isequal(x, [1 0 0; ...
                        1 1 0; ...
                        0 0 0]);
sNW2 = @(x) isequal(x, [1 1 0; ...
                        0 1 0; ...
                        0 0 0]);
sSW1 = @(x) isequal(x, [0 0 0; ...
                        1 1 0; ...
                        1 0 0]);
sSW2 = @(x) isequal(x, [0 0 0; ...
                        0 1 0; ...
                        1 1 0]);
sSE1 = @(x) isequal(x, [0 0 0; ...
                        0 1 1; ...
                        0 0 1]);
sSE2 = @(x) isequal(x, [0 0 0; ...
                        0 1 0; ...
                        0 1 1]);
                    
% pits
p4E = @(x) isequal(x, [1 1 0; ...
                       1 0 0; ...
                       1 1 0]);
p4N = @(x) isequal(x, [0 0 0; ...
                       1 0 1; ...
                       1 1 1]);
p4W = @(x) isequal(x, [0 1 1; ...
                       0 0 1; ...
                       0 1 1]);
p4S = @(x) isequal(x, [1 1 1; ...
                       1 0 1; ...
                       0 0 0]);
p8NE = @(x) isequal(x, [1 0 0; ...
                        1 0 0; ...
                        1 1 1]);
p8NW = @(x) isequal(x, [0 0 1; ...
                        0 0 1; ...
                        1 1 1]);
p8SW = @(x) isequal(x, [1 1 1; ...
                        0 0 1; ...
                        0 0 1]);
p8SE = @(x) isequal(x, [1 1 1; ...
                        1 0 0; ...
                        1 0 0]);
                   
lut = makelut(@(x) n4N(x) | n4W(x) | n4S(x) | n4E(x) | ...
          n8NE(x) | n8NW(x) | n8SW(x) | n8SE(x) | ...
          sNE1(x) | sNE2(x) | sNW1(x) | sNW2(x) | ...
          sSW1(x) | sSW2(x) | sSE1(x) | sSE2(x) | ...
          p4E(x) | p4N(x) | p4W(x) | p4S(x) | ...
          p8NE(x) | p8NW(x) | p8SW(x) | p8SE(x), 3);

end
