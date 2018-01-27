function [initMLPhi] = createInitMLPhiSquares(Nx, Ny, parameters, layerNum)

d = parameters(1).d;
d = int32(round(d/2));

if layerNum < 2
    layerNum = 2;
end

% bgSquare = zeros(2*d, 2*d) - 1;
bgSquare = zeros(2*d, 2*d) ;
fgSquare = zeros(2*d, 2*d) + 1;

if layerNum==2
    
    initMLPhi = zeros(Nx, Ny, layerNum);
    
    tile = [bgSquare fgSquare; fgSquare bgSquare];
    [hx,wx] = size(tile);
    
    xSize = round( Nx / hx );
    ySize = round( Ny / wx );
    
    layer1 = repmat(tile, xSize+1, ySize+1);
    
    initMLPhi(:,:,1) = layer1(1:Nx, 1:Ny);
    
    tile = [fgSquare bgSquare; bgSquare fgSquare];
    
    [hx,wx] = size(tile);
    
    xSize = round( Nx / hx );
    ySize = round( Ny / wx );
    
    layer2 = repmat(tile, xSize+1, ySize+1);
    
    initMLPhi(:,:,2) = layer2(1:Nx, 1:Ny);
    
else
    
%     tile = zeros(2*d, 2*d*layerNum) - 1;
    tile = zeros(2*d, 2*d*layerNum) ;
    tile(:, 1:2*d) = 1;
    
    xSize = round( Nx/(2*d*layerNum) );
    
    bigTileRow1 = repmat( tile, 1, xSize+2);
    
%     bigTileRow2 = zeros(size(bigTileRow1)) - 1;
    bigTileRow2 = zeros(size(bigTileRow1)) ;
    
    bigTileRow2(:, 3*d+1:end) = bigTileRow1( :, 1:end-3*d );
    
    ySize = round( Ny/(4*d) );
    
    patternAll = repmat( [bigTileRow1; bigTileRow2], ySize+1, 1);
    
    initMLPhi = zeros(Nx, Ny, layerNum);
    
    for i=1:layerNum
%         Nx+(i-1)*2*d
        initMLPhi(:,:,i) = patternAll(1:Nx, 1+(i-1)*2*d:Ny+(i-1)*2*d);
    end
       
end

end