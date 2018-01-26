function bestPlane = selectFocusPlane(planeImagesStack)
% choose the best focused image based on total image gradient

[~,~,planeNum] = size(planeImagesStack);

tg = zeros(planeNum,1);

for i=1:planeNum
    in = planeImagesStack(:,:,i);
    [dx, dy] = gradient(double(in));
    tg(i) = sum(dx(:).^2) + sum(dy(:).^2);
end
    
[~,bestPlaneIdx] = max(tg);

bestPlane = squeeze(planeImagesStack(:,:,bestPlaneIdx));