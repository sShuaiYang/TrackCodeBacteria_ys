function mask= inMaskShow(bioTree,pixelIdxListIn,i)
xSize=bioTree{1}.imageSize(1);
ySize=bioTree{1}.imageSize(2);
mask=false(xSize,ySize);
mask(pixelIdxListIn{i})=true;
figure;
imshow(mask);
end