%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Reads (opens) the image you want to analyze and assigns it to a
%%% variable.
OrigImage = CPretrieveimage(handles,ImageName,ModuleName,'MustBeGray','CheckScale');

[xs, ys] = size(OrigImage);

%% compute the 'a trous' wavelet transform and compute the planes product,
%% find local maxima
aTrousLevel = 2;
w = a_trous(-OrigImage, aTrousLevel+1);

%product(:,:) = w(:,:, aTrousLevel);%prod(w(:,:, :), 3);
product = prod(w(:,:,  aTrousLevel:aTrousLevel), 3);

% remove noise
spotthres = mean(product(:)) + NoiseRemovalFactor*std(product(:));

noiseValues = find(product < spotthres);

product(noiseValues) = 0;

filter = fspecial('disk', SizeOfSmoothingFilter);

product = imfilter(product, filter);

filter = fspecial('gaussian', SizeOfSmoothingFilter);

product = imfilter(product, filter);

[XMAX,IMAX,XMIN,IMIN] = extrema2(product);

%disp(['Median spot intensity: ' num2str(median(product(IMAX)))]);

llist = find(product(IMAX) > FinalSpotThreshold);

IMAX = IMAX(llist);


%% cleaning

% this step could be improved by index trasform
mask = zeros(size(OrigImage));
mask(IMAX) = 1;
[r c] = find(mask == 1);

% delete spots on the bounday
boundary = 10;
todel = find(c<boundary); r(todel) = []; c(todel) = [];
todel = find(r<boundary); r(todel) = []; c(todel) = [];
todel = find(c>ys-boundary); r(todel) = []; c(todel) = [];
todel = find(r>xs-boundary); r(todel) = []; c(todel) = [];

% 
% % shrink to 1 pixel
% ShrinkMatrix = logical(size(OrigImage));
% 
% for i=1:length(r)
%     ShrinkMatrix(r(i), c(i)) = 1;
% end;
% 
% ShrinkMatrix2 = bwmorph(ShrinkMatrix,'shrink', Inf);
% 
% [r c] = find(ShrinkMatrix2 == 1);

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

FinalLabelMatrixImage  = zeros(size(OrigImage));

for i=1:length(r)
    FinalLabelMatrixImage(r(i), c(i)) = i;
end

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
        LogicalOutlines = FinalLabelMatrixImage > 0;
        handles.Pipeline.(SaveOutlines) = LogicalOutlines;
    end
catch
    error(['The object outlines were not calculated by the ', ModuleName, ' module, so these images were not saved to the handles structure. The Save Images module will therefore not function on these images. This is just for your information - image processing is still in progress, but the Save Images module will fail if you attempted to save these images.'])
end