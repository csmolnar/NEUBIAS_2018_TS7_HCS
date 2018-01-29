function handles = DetectAllSpots(handles)

% Help for the Detect All Spots module:
% Category: Object Processing
%
% SHORT DESCRIPTION:
% Simplified spot detection function based on A-Trous wavelet transform and
% Hough transformation... 
% TODO: write this once the algorithm is stable and clean.
% *************************************************************************
%
% This module detect spots on the image based on the A-Trous wavelet
% transform.
%
% Settings:
%
% See also Identify primary Identify Secondary modules.

% CellProfiler is distributed under the GNU General Public License.
% See the accompanying file LICENSE for details.
%
% Developed by Peter Horvath 2010 and Abel Szkalisity 2016.

% $Revision: 5025 $

%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%
%drawnow

[~, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = What did you call the images you want to process?
%infotypeVAR01 = imagegroup
ImageName = char(handles.Settings.VariableValues{CurrentModuleNum,1});
%inputtypeVAR01 = popupmenu

%textVAR02 = What do you want to call the spots identified by this module?
%defaultVAR02 = Spots
%infotypeVAR02 = objectgroup indep
SpotObjectName = char(handles.Settings.VariableValues{CurrentModuleNum,2});

%textVAR03 = What is the radius of the smallest droplet that you want to detect?
%defaultVAR03 = 0
minSpotRadius = str2double(handles.Settings.VariableValues{CurrentModuleNum,3});
minSpotSize = minSpotRadius.^2 *pi;

%textVAR04 = What is the radius of the largest droplet that you want to detect?
%defaultVAR04 = 40
maxSpotRadius = str2double(handles.Settings.VariableValues{CurrentModuleNum,4});
maxSpotSize = maxSpotRadius.^2 *pi;

%textVAR05 = Noise removal factor (removes background + n times std).
%defaultVAR05 = 1.5
NoiseRemovalFactor = str2double(handles.Settings.VariableValues{CurrentModuleNum,5});

%textVAR06 = Static threshold for spot identificaiton.
%defaultVAR06 = 1e-010
FinalSpotThreshold = str2double(handles.Settings.VariableValues{CurrentModuleNum,6});

%textVAR07 = Define a circularity threshold for droplets!
%defaultVAR07 = 0.8
circularityThreshold = str2double(handles.Settings.VariableValues{CurrentModuleNum,7});

%textVAR08 = Specify a name for the outlines of the identified objects! (optional)
%defaultVAR08 = Do not save
%infotypeVAR08 = outlinegroup indep
SaveOutlines = char(handles.Settings.VariableValues{CurrentModuleNum,8});


%%%VariableRevisionNumber = 3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRELIMINARY CALCULATIONS & FILE HANDLING %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%drawnow

%%% Reads (opens) the image you want to analyze and assigns it to a
%%% variable.
OrigImage = CPretrieveimage(handles,ImageName,ModuleName,'MustBeGray','CheckScale');

[xs, ys] = size(OrigImage);
origMean = mean(OrigImage(:));
origStd = std(OrigImage(:));

imFilterSize = 5;
largestWaveletSize = 5; % It is empricially determined that the 5 level wavelet responds to quite big circles with radius ~15-16.
%another reason for choosing 5 as an upper limit is that greater values
%cause very large a-trous filter size which means that completely flat
%regions are also affected by far away droplets.

%waveletMaxSpotSize = min(maxSpotSize,largestWaveletSize);
%maxLevelCalc = (waveletMaxSpotSize-imFilterSize)/(imFilterSize-1);%this is the size of the gap
%maxLevelCalc = ceil(sqrt(maxLevelCalc)); %this is the input level size

%empricially determined wavelet response limits for size. Entry indices are
%the wavelet scales, starting from 2.
%sizeLimits = [-1,21,185,587,855];
sizeLimits = [-1,41,385,787,1055];

% objects with circularity below the threshold are removed.
%the circularity can go above 1 becuase the pixel space is discrete

w = a_trous(-OrigImage, largestWaveletSize,imFilterSize);

FinalLabelMatrixImage  = zeros(size(OrigImage));
identifiedRegionpropsPerLevel = cell(largestWaveletSize,3); % store area and perimeter next to it

%mask for shrinking
sr = strel('disk',1);

%Thresholded image for droplets area detection
%keepThreshold = 0.5;
%keepThreshold = 0.85;
%threshImage = thresholdPercentile(OrigImage,keepThreshold);
level = graythresh(OrigImage);
threshImage= im2bw(OrigImage,level);

for aTrousLevel = 2:largestWaveletSize
    
    %This will store the center points for the current scale level
    CurrentLabelMatrixImage = zeros(size(FinalLabelMatrixImage));
        
    product = w(:,:,  aTrousLevel);

    % coarse threshold
    spotthres = mean(product(:)) + NoiseRemovalFactor*std(product(:));   
    product(product < spotthres) = 0;   

    SizeOfSmoothingFilter = aTrousLevel-1;
    
    % shape based enhancement
    filter = fspecial('disk', SizeOfSmoothingFilter);
    product = imfilter(product, filter);
    %intensity based enhancement
    filter = fspecial('gaussian', SizeOfSmoothingFilter);
    product = imfilter(product, filter);
    
    %{
    identifiedRegionLabels = product;
    identifiedRegionLabels(product>0) = 1;
    CC = bwconncomp(identifiedRegionLabels);
    for i=1:CC.NumObjects
        identifiedRegionLabels(CC.PixelIdxList{i}) = i;
    end
    %}

    [~,IMAX,~,~] = extrema2(product);

    %REMOVE static threshold    
    IMAX = IMAX(product(IMAX) > FinalSpotThreshold);


    %% cleaning

    % this step could be improved by index trasform
    mask = zeros(size(OrigImage));
    mask(IMAX) = 1;
    [r c] = find(mask == 1);

    % delete spots on the boundary
    boundary = 10;
    todel = find(c<boundary); r(todel) = []; c(todel) = [];
    todel = find(r<boundary); r(todel) = []; c(todel) = [];
    todel = find(c>ys-boundary); r(todel) = []; c(todel) = [];
    todel = find(r>xs-boundary); r(todel) = []; c(todel) = [];


    %% remove spots too close to one another
    minDist = 2^aTrousLevel;
    minDistS = minDist^2;

    todel = [];

    filter = fspecial('disk', 3);

    avgint = imfilter(OrigImage, filter, 'replicate');

    for i=1:length(c)
        for j=i+1:length(c)
            dist = (r(i) - r(j))^2 + (c(i) - c(j))^2;
            if dist < minDistS
                if (avgint(r(i), c(i)) < avgint(r(j), c(j)))
                    todel = [todel i];
                else
                    todel = [todel j];                
                end;
            end;
        end;
    end;

    r(todel) = []; c(todel) = [];
    

    for i=1:length(r)
        CurrentLabelMatrixImage(r(i), c(i)) = i;
    end
    
    %Identify regions around the objects (???)
    
    %Identify automatically regions around the objects using CP's function.
    %propagatedImage = IdentifySecPropagateSubfunction(param1_primaryObjectsImageTheObjectsAreIdentifiedWithLabels,param2_TheOriginalImageOnWhichThePropagationWillBeDone,...
    %param3_ThresholdedOriginalImage,param4_regularizationFactorWhichTellsThatTheDividingLinesAreDrawnBasedOnIntensityOfTheSecondaryImageORBasedOnDistance);    
    
    [identifiedRegionLabelsCP,distances] = IdentifySecPropagateSubfunction(CurrentLabelMatrixImage,OrigImage,threshImage,0.05);
    %reduce thresholded boundaries by variance maxima detection
    [~,~,identifiedRegionLabels] = plotVarianceChange(identifiedRegionLabelsCP,distances,OrigImage,CurrentLabelMatrixImage,0);    
    
    %TODO: the filtering can be done on the initial variable don't have to
    %copy
    filteredIdentifiedRegionLabels = identifiedRegionLabels;
    %filter out by identified size and circularity
    r = regionprops(identifiedRegionLabels,{'Area','Perimeter'});
    area = cat(1,r.Area);
    perim = cat(1,r.Perimeter);
    indicesToDelete = find(area>sizeLimits(aTrousLevel) | area>maxSpotSize | area<minSpotSize);
    if aTrousLevel>2
        perim(perim == 0) = 1;
        indicesToDelete = union(indicesToDelete,find((area.*4.*pi ./ perim.^2)<circularityThreshold));
    end
    filteredIdentifiedRegionLabels(ismember(filteredIdentifiedRegionLabels,indicesToDelete)) = 0;
    fileteredArea = area(setdiff(1:length(area),indicesToDelete));
    filteredPerimeter = perim(setdiff(1:length(area),indicesToDelete));
    
    identifiedRegionpropsPerLevel{aTrousLevel,1} = filteredIdentifiedRegionLabels;
    identifiedRegionpropsPerLevel{aTrousLevel,2} = fileteredArea;
    identifiedRegionpropsPerLevel{aTrousLevel,3} = filteredPerimeter;
    
        
end

%The final label matrix image. First fill it out with the perfect circles.
FinalLabelMatrixImage = zeros(size(OrigImage));

%The constantly growing label counter
actualLabel = 1;

%Detect perfectly circular objects by Hough transform. (Hough transform works reliably only on radii greater than 10.)
radiiRanges = 10:10:(maxSpotRadius-1);
radiiRanges(end+1) = maxSpotRadius;
centers = cell(1,length(radiiRanges));
radii = cell(1,length(radiiRanges));
for i=1:length(radiiRanges)-1
    [centers{i},radii{i}] = imfindcircles(OrigImage,radiiRanges(i:i+1));
    for j=1:size(centers{i},1)
        mask = drawCircle(zeros(size(OrigImage)),centers{i}(j,2),centers{i}(j,1),radii{i}(j),1);
        intensities = OrigImage(logical(mask));
        if mean(intensities) > origMean + origStd * NoiseRemovalFactor
            FinalLabelMatrixImage = drawCircle(FinalLabelMatrixImage,centers{i}(j,2),centers{i}(j,1),radii{i}(j),actualLabel);
            actualLabel = actualLabel + 1;
        end
    end
end

%after identifying spots on each level remove overlaps from top to bottom.

for i=largestWaveletSize:-1:2
    labels = unique(identifiedRegionpropsPerLevel{i,1});
    if labels(1) == 0
        labels(1) = [];
    end
    for j=1:length(labels)
        mask = identifiedRegionpropsPerLevel{i,1} == labels(j);
        newCirc = identifiedRegionpropsPerLevel{i,2}(j).*4.*pi./(identifiedRegionpropsPerLevel{i,3}(j).^2);
        %shrinked = imerode(mask,sr);
        %newCirc = circularityByMask(mask,shrinked);
        %if the area is clean yet
        if ~any(FinalLabelMatrixImage(:) & mask(:))
            FinalLabelMatrixImage(mask) = actualLabel;
            actualLabel = actualLabel + 1;
        %if there is an overlap       
        elseif i>2
            overlappingObjectIndices = unique(FinalLabelMatrixImage.*mask);
            overlappingObjectIndices = overlappingObjectIndices(2:end);
            circularities = zeros(1,length(overlappingObjectIndices));
            oldMask = zeros(size(FinalLabelMatrixImage));
            for k=1:length(overlappingObjectIndices)
                actMask = FinalLabelMatrixImage == overlappingObjectIndices(k);
                shrinked = imerode(actMask,sr);
                circularities(k) = circularityByMask(actMask,shrinked);
                oldMask = oldMask | actMask;
            end            
            %if the new object is more circular than the old ones then
            if newCirc>max(circularities) 
                %remove old objects                
                FinalLabelMatrixImage(oldMask) = 0;
                FinalLabelMatrixImage(mask) = actualLabel;
                actualLabel = actualLabel + 1;
            end
        end
    end
end

%Make the outlines.
StructuringElement = strel('square',3);
bitBiggerObjects = imdilate(FinalLabelMatrixImage,StructuringElement);
outlines = bitBiggerObjects - FinalLabelMatrixImage;
LogicalOutlines = outlines > 0;    

% Story END and CP variable savings

if ~isfield(handles.Measurements,SpotObjectName)
    handles.Measurements.(SpotObjectName) = {};
end

%%% Saves the final, segmented label matrix image of secondary objects to
%%% the handles structure so it can be used by subsequent modules.
fieldname = ['Segmented',SpotObjectName];
handles.Pipeline.(fieldname) = FinalLabelMatrixImage;

fieldname = ['SmallRemovedSegmented',SpotObjectName];
handles.Pipeline.(fieldname) = FinalLabelMatrixImage;

handles = CPsaveObjectCount(handles, SpotObjectName, FinalLabelMatrixImage);
handles = CPsaveObjectLocations(handles, SpotObjectName, FinalLabelMatrixImage);

%%% Saves images to the handles structure so they can be saved to the hard
%%% drive, if the user requested.
try
    if ~strcmpi(SaveOutlines,'Do not save')        
        handles.Pipeline.(SaveOutlines) = LogicalOutlines;
    end
catch e
    error(['The object outlines were not calculated by the ', ModuleName, ' module, so these images were not saved to the handles structure. The Save Images module will therefore not function on these images. This is just for your information - image processing is still in progress, but the Save Images module will fail if you attempted to save these images.'])
end


%A_trous transform
function w = a_trous(in, level,filterSize)

in = double(in);

[sizex sizey] = size(in);

c = zeros(sizex, sizey, level);

w = zeros(sizex, sizey, level);

c(:, :, 1) = in;

for i=2:level

    c(:, :, i) = calculateC(in, i-2,filterSize);

    w(:,:,i-1) = c(:, :, i) - c(:, :, i-1);

end;

w(:,:,i) = c(:,:,i);


function out = calculateC(in, level,filterSize)

[sizex sizey] = size(in);

rlevel = 2^level;

if sizex < 2 * rlevel || sizey < 2 * rlevel
    
    disp('Too huge level value');
    
    return;
    
end;

%% boundary mirroring
inb = zeros(sizex + 4*rlevel, sizey + 4*rlevel);
inb(2*rlevel+1:2*rlevel+sizex, 2*rlevel+1:2*rlevel+sizey) = in;

inb(2*rlevel+1:2*rlevel+sizex, 1:2*rlevel) = in(1:sizex, 2*rlevel:-1:1);
inb(2*rlevel+1:2*rlevel+sizex, 2*rlevel+sizey+1:sizey+4*rlevel) = in(1:sizex, sizey:-1:sizey-2*rlevel+1);

inb(1:2*rlevel, 1:sizey + 4*rlevel) = inb(4*rlevel:-1:2*rlevel+1, 1:sizey + 4*rlevel);
inb(2*rlevel+sizex+1:sizex+4*rlevel, 1:sizey + 4*rlevel) = inb(2*rlevel+sizex:-1:sizex+1, 1:sizey + 4*rlevel);
filter = createfilter_2d(1, filterSize, level^2);
%Why the convolution is done on the in image? Or at least then why the inb
%computed at all?
out = conv2(double(inb),filter,'same');

% cut borders

cutsize = (size(filter, 1)-1) / 2;

out = out(cutsize+1:cutsize+sizex, cutsize+1:cutsize+sizey);


function filt = createfilter_2d(sigma, size, gap)

if nargin == 2
    gap = 0;
end;

filt = zeros((size-1)*(gap+1)+1, (size-1)*(gap+1)+1);

for i=1:size
    for j=1:size
        middle = size/2;
        dist = sqrt((middle-i+0.5)^2+(middle-j+0.5)^2);
        filt((i-1)*(gap+1)+1, (j-1)*(gap+1)+1) = normpdf(dist, 0, sigma);
    end;
end;

filt = filt ./ sum(filt(:));

%thresold the image according to a percentile to remove from the low
%intensity region
function threshImage = thresholdPercentile(image,percent)
    [~,idx] = sort(image(:));
    threshold = image(idx(round(numel(image)*percent)));
    threshImage = image > threshold;
    

function [xmax,imax,xmin,imin] = extrema(x)
%EXTREMA   Gets the global extrema points from a time series.
%   [XMAX,IMAX,XMIN,IMIN] = EXTREMA(X) returns the global minima and maxima
%   points of the vector X ignoring NaN's, where
%    XMAX - maxima points in descending order
%    IMAX - indexes of the XMAX
%    XMIN - minima points in descending order
%    IMIN - indexes of the XMIN
%
%   DEFINITION (from http://en.wikipedia.org/wiki/Maxima_and_minima):
%   In mathematics, maxima and minima, also known as extrema, are points in
%   the domain of a function at which the function takes a largest value
%   (maximum) or smallest value (minimum), either within a given
%   neighbourhood (local extrema) or on the function domain in its entirety
%   (global extrema).
%
%   Example:
%      x = 2*pi*linspace(-1,1);
%      y = cos(x) - 0.5 + 0.5*rand(size(x)); y(40:45) = 1.85; y(50:53)=NaN;
%      [ymax,imax,ymin,imin] = extrema(y);
%      plot(x,y,x(imax),ymax,'g.',x(imin),ymin,'r.')
%
%   See also EXTREMA2, MAX, MIN

%   Written by
%   Lic. on Physics Carlos Adrián Vargas Aguilera
%   Physical Oceanography MS candidate
%   UNIVERSIDAD DE GUADALAJARA
%   Mexico, 2004
%
%   nubeobscura@hotmail.com

% From       : http://www.mathworks.com/matlabcentral/fileexchange
% File ID    : 12275
% Submited at: 2006-09-14
% 2006-11-11 : English translation from spanish.
% 2006-11-17 : Accept NaN's.
% 2007-04-09 : Change name to MAXIMA, and definition added.


xmax = [];
imax = [];
xmin = [];
imin = [];

% Vector input?
Nt = numel(x);
if Nt ~= length(x)
    error('Entry must be a vector.')
end

% NaN's:
inan = find(isnan(x));
indx = 1:Nt;
if ~isempty(inan)
    indx(inan) = [];
    x(inan) = [];
    Nt = length(x);
end

% Difference between subsequent elements:
dx = diff(x);

% Is an horizontal line?
if ~any(dx)
    return
end

% Flat peaks? Put the middle element:
a = find(dx~=0);              % Indexes where x changes
lm = find(diff(a)~=1) + 1;    % Indexes where a do not changes
d = a(lm) - a(lm-1);          % Number of elements in the flat peak
a(lm) = a(lm) - floor(d/2);   % Save middle elements
a(end+1) = Nt;

% Peaks?
xa  = x(a);             % Serie without flat peaks
b = (diff(xa) > 0);     % 1  =>  positive slopes (minima begin)
% 0  =>  negative slopes (maxima begin)
xb  = diff(b);          % -1 =>  maxima indexes (but one)
% +1 =>  minima indexes (but one)
imax = find(xb == -1) + 1; % maxima indexes
imin = find(xb == +1) + 1; % minima indexes
imax = a(imax);
imin = a(imin);

nmaxi = length(imax);
nmini = length(imin);

% Maximum or minumim on a flat peak at the ends?
if (nmaxi==0) && (nmini==0)
    if x(1) > x(Nt)
        xmax = x(1);
        imax = indx(1);
        xmin = x(Nt);
        imin = indx(Nt);
    elseif x(1) < x(Nt)
        xmax = x(Nt);
        imax = indx(Nt);
        xmin = x(1);
        imin = indx(1);
    end
    return
end

% Maximum or minumim at the ends?
if (nmaxi==0)
    imax(1:2) = [1 Nt];
elseif (nmini==0)
    imin(1:2) = [1 Nt];
else
    if imax(1) < imin(1)
        imin(2:nmini+1) = imin;
        imin(1) = 1;
    else
        imax(2:nmaxi+1) = imax;
        imax(1) = 1;
    end
    if imax(end) > imin(end)
        imin(end+1) = Nt;
    else
        imax(end+1) = Nt;
    end
end
xmax = x(imax);
xmin = x(imin);

% NaN's:
if ~isempty(inan)
    imax = indx(imax);
    imin = indx(imin);
end

% Same size as x:
imax = reshape(imax,size(xmax));
imin = reshape(imin,size(xmin));

% Descending order:
[~,inmax] = sort(-xmax); clear temp
xmax = xmax(inmax);
imax = imax(inmax);
[xmin,inmin] = sort(xmin);
imin = imin(inmin);


% Carlos Adrián Vargas Aguilera. nubeobscura@hotmail.com

function [xymax,smax,xymin,smin] = extrema2(xy,varargin)
%EXTREMA2   Gets the extrema points from a surface.
%   [XMAX,IMAX,XMIN,IMIN] = EXTREMA2(X) returns the maxima and minima
%   elements of the matriz X ignoring NaN's, where
%    XMAX - maxima points in descending order (the bigger first and so on)
%    IMAX - linear indexes of the XMAX
%    XMIN - minima points in descending order
%    IMIN - linear indexes of the XMIN.
%   The program uses EXTREMA.
%
%   The extrema points are searched only through the column, the row and
%   the diagonals crossing each matrix element, so it is not a perfect
%   mathematical program and for this reason it has an optional argument.
%   The user should be aware of these limitations.
%
%   [XMAX,IMAX,XMIN,IMIN] = EXTREMA2(X,1) does the same but without
%   searching through the diagonals (less strict and perhaps the user gets
%   more output points).
%
%   DEFINITION (from http://en.wikipedia.org/wiki/Maxima_and_minima):
%   In mathematics, maxima and minima, also known as extrema, are points in
%   the domain of a function at which the function takes a largest value
%   (maximum) or smallest value (minimum), either within a given
%   neighbourhood (local extrema) or on the function domain in its entirety
%   (global extrema).
%
%   Note: To change the linear index to (i,j) use IND2SUB.
%
%   Example:
%      [x,y] = meshgrid(-2:.2:2,3:-.2:-2);
%      z = x.*exp(-x.^2-y.^2); z(10,7)= NaN; z(16:19,13:17) = NaN;
%      surf(x,y,z), shading interp
%      [zmax,imax,zmin,imin] = extrema2(z);
%      hold on
%       plot3(x(imax),y(imax),zmax,'bo',x(imin),y(imin),zmin,'ro')
%       for i = 1:length(zmax)
%        text(x(imax(i)),y(imax(i)),zmax(i),['  ' num2str(zmax(i))])
%       end
%       for i = 1:length(zmin)
%        text(x(imin(i)),y(imin(i)),zmin(i),['  ' num2str(zmin(i))])
%       end
%      hold off
%
%   See also EXTREMA, MAX, MIN

%   Written by
%   Lic. on Physics Carlos Adrián Vargas Aguilera
%   Physical Oceanography MS candidate
%   UNIVERSIDAD DE GUADALAJARA
%   Mexico, 2005
%
%   nubeobscura@hotmail.com

% From       : http://www.mathworks.com/matlabcentral/fileexchange
% File ID    : 12275
% Submited at: 2006-09-14
% 2006-11-11 : English translation from spanish.
% 2006-11-17 : Accept NaN's.
% 2006-11-22 : Fixed bug in INDX (by JaeKyu Suhr)
% 2007-04-09 : Change name to MAXIMA2, and definition added.

M = size(xy);
if length(M) ~= 2
    error('Entry must be a matrix.')
end
N = M(2);
M = M(1);

% Search peaks through columns:
[smaxcol,smincol] = extremos(xy);

if isempty(smaxcol)
    xymax = []; 
    smax = []; 
    xymin = []; 
    smin = []; 
    return;
end;

% Search peaks through rows, on columns with extrema points:
im = unique([smaxcol(:,1);smincol(:,1)]); % Rows with column extrema
[smaxfil,sminfil] = extremos(xy(im,:).');

% Convertion from 2 to 1 index:
smaxcol = sub2ind([M,N],smaxcol(:,1),smaxcol(:,2));
smincol = sub2ind([M,N],smincol(:,1),smincol(:,2));
smaxfil = sub2ind([M,N],im(smaxfil(:,2)),smaxfil(:,1));
sminfil = sub2ind([M,N],im(sminfil(:,2)),sminfil(:,1));

% Peaks in rows and in columns:
smax = intersect(smaxcol,smaxfil);
smin = intersect(smincol,sminfil);

% Search peaks through diagonals?
if nargin==1
    % Check peaks on down-up diagonal:
    [iext,jext] = ind2sub([M,N],unique([smax;smin]));
    [sextmax,sextmin] = extremos_diag(iext,jext,xy,1);

    % Check peaks on up-down diagonal:
    smax = intersect(smax,[M; (N*M-M); sextmax]);
    smin = intersect(smin,[M; (N*M-M); sextmin]);

    % Peaks on up-down diagonals:
    [iext,jext] = ind2sub([M,N],unique([smax;smin]));
    [sextmax,sextmin] = extremos_diag(iext,jext,xy,-1);

    % Peaks on columns, rows and diagonals:
    smax = intersect(smax,[1; N*M; sextmax]);
    smin = intersect(smin,[1; N*M; sextmin]);
end

% Extrema points:
xymax = xy(smax);
xymin = xy(smin);

% Descending order:
[~,inmax] = sort(-xymax); clear temp
xymax = xymax(inmax);
smax = smax(inmax);
[xymin,inmin] = sort(xymin);
smin = smin(inmin);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [smax,smin] = extremos(matriz)
% Peaks through columns or rows.

smax = [];
smin = [];

for n = 1:length(matriz(1,:))
    [~,imaxfil,~,iminfil] = extrema(matriz(:,n)); clear temp
    if ~isempty(imaxfil)     % Maxima indexes
        imaxcol = repmat(n,length(imaxfil),1);
        smax = [smax; imaxfil imaxcol];
    end
    if ~isempty(iminfil)     % Minima indexes
        imincol = repmat(n,length(iminfil),1);
        smin = [smin; iminfil imincol];
    end
end


function [sextmax,sextmin] = extremos_diag(iext,jext,xy,A)
% Peaks through diagonals (down-up A=-1)

[M,N] = size(xy);
if A==-1
    iext = M-iext+1;
end
[iini,jini] = cruce(iext,jext,1,1);
[iini,jini] = ind2sub([M,N],unique(sub2ind([M,N],iini,jini)));
[ifin,jfin] = cruce(iini,jini,M,N);
sextmax = [];
sextmin = [];
for n = 1:length(iini)
    ises = iini(n):ifin(n);
    jses = jini(n):jfin(n);
    if A==-1
        ises = M-ises+1;
    end
    s = sub2ind([M,N],ises,jses);
    [~,imax,~,imin] = extrema(xy(s)); clear temp
    sextmax = [sextmax; s(imax)'];
    sextmin = [sextmin; s(imin)'];
end


function [i,j] = cruce(i0,j0,I,J)
% Indexes where the diagonal of the element io,jo crosses the left/superior
% (I=1,J=1) or right/inferior (I=M,J=N) side of an MxN matrix.

arriba = 2*(I*J==1)-1;

si = (arriba*(j0-J) > arriba*(i0-I));
i = (I - (J+i0-j0)).*si + J+i0-j0;
j = (I+j0-i0-(J)).*si + J;


% Carlos Adrián Vargas Aguilera. nubeobscura@hotmail.com

function [vari,img,newRegionLabels] = plotVarianceChange(identifiedRegionLabels,distances,intensityImage,centerImage,toPlot)
%script to plot out variance change inside the identified particles

    %for all identified subregion
    sorted = sort(centerImage(:));
    positives = sorted(sorted>0);
    vari = cell(1,length(positives));
    k = 1;
    img = zeros(size(intensityImage));
    newRegionLabels = zeros(size(intensityImage));
    for i=positives'                 
        origIdx = find(identifiedRegionLabels == i);    
        dist = distances(origIdx);
        inten = intensityImage(origIdx);
        [~,sortedIdx] = sort(dist);
        vari{k} = zeros(1,length(sortedIdx));                   
        
        actMean = inten(sortedIdx(1));        
        soFarDiv = 0;
        vari{k}(1) = 0;
        for j=2:length(vari{k})
            [soFarDiv,actMean] = onlineDiffFromMeanSquared(soFarDiv,actMean,inten(sortedIdx(j)),j);
            vari{k}(j) = soFarDiv / j;
            img(origIdx(sortedIdx(j))) = vari{k}(j);
        end        
                
        %find the first local maxima 
        %maxIndices = localMaxima(vari{k},50);
        %if isempty(maxIndices)
        %    maxIndices = length(vari{k});
        %end
        [~,maxIndices] = max(vari{k});
        %The maximum variance region detection seems to overshoot a bit,
        %therefore we emprically identified a reduction factor, which is a
        %linear function fitted to the following training data:
        % x = [14,29,34,35,97,459];
        % y = [14,18,16,22,43,154];
        % p = polyfit(x,y,1);
        p = [0.316238724039588 9.292088723592572];
        maxIndicesReduced = min(round(polyval(p,maxIndices)),maxIndices);
        newRegionLabels(origIdx(sortedIdx(1:maxIndicesReduced(1)))) = i;
        
        %plot if requested
        if toPlot
            f = figure();
            mask = identifiedRegionLabels == i;
            imshow(mask);
            g = figure();
            plot(vari{k});
            s = input('Do you want to save this figure? [y/n]','s');
            if s=='y'
                saveas(gcf,['variPlot_' num2str(i,'%03d') '.png']);
            elseif s=='s'
                break;
            end
            close(g);
            close(f);
        end
        k = k+1;
    end

function newMean = onlineMean(prevMean,xn,n)    
    newMean = prevMean + (xn-prevMean)/n;
    
function [newM,xnA] = onlineDiffFromMeanSquared(prevM,prevMean,xn,n)
    xnA = onlineMean(prevMean,xn,n);
    newM = prevM + (xn-prevMean)*(xn-xnA);    

 function [circularity,shrinked] = circularityByMask(mask,shrinked,area)
%Circularity is measured by area*4*pi / perimeter.^2. Theoretically this
%value is 1 for the circle and smaller for all other graphical object but
%in practice it can go above 1.
%INPUT:
%   mask        A binary image with the object
%   area        The area of the object can be given to fasten the
%               calculation but it is optional.
%   shirnked    A shrinked version of the mask, for calculating perimeter.
%   
%OUTPUT:
%   circularity The circularity value of the binary object in mask.

if nargin<3
    area = sum(mask(:));
end

contours = mask-shrinked;
contourChain = find(contours);
perimeter = length(contourChain);
circularity = area*4*pi / (perimeter.^2);

function img = drawCircle(img,x,y,r,numberToFill)
%fills up with one if the last parameter is not specified
%in case of index out of bound it just jumps.

if nargin<5
    numberToFill = 1;
end

alpha = 0:0.01:2*pi;
for i=1:length(alpha)
    for j = 0:r/100:r        
        if (round(x+j*cos(alpha(i))) <= size(img,1) && round(x+j*cos(alpha(i)))>0  && round(y+j*sin(alpha(i)))<=size(img,2) && round(y+j*sin(alpha(i)))>0)
            img(round(x+j*cos(alpha(i))),round(y+j*sin(alpha(i)))) = numberToFill;
        end
    end
end
    