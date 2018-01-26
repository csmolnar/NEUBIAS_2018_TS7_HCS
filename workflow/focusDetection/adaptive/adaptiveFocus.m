function focusedCompositeImage = adaptiveFocus(planeImagesStack)
% composite image preparations step

global options;

if isfield(options.segmentation, 'expectedObjectSize')
    expectedObjectSize = options.segmentation.expectedObjectSize;
else
    expectedObjectSize = 30;
end

if isempty(planeImagesStack)
    error('HCS:focus','Undefined input planes.');
end

gaussianFilter = fspecial('gaussian', [7*expectedObjectSize 7*expectedObjectSize], expectedObjectSize);

[h,w,planeNum] = size(planeImagesStack);

convolvedMap = zeros(h,w,planeNum);

for i=1:planeNum
    in = planeImagesStack(:,:,i);
    [dx, dy] = gradient(double(in));
    gradMap = dx.^2 + dy.^2;
    convolvedMap(:, :, i) = imfilter(gradMap, gaussianFilter);
end

[~,depthMap] = max(convolvedMap, [], 3);

focusedCompositeImage = cast(zeros(h,w),class(planeImagesStack));

for i=1:h
    for j=1:w
        focusedCompositeImage(i, j) = planeImagesStack(i, j, depthMap(i, j));
    end
end

end