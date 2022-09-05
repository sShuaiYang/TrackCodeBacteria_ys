function mask= outMaskShow(bioTree,pixelIdxListOut,i)
xSize=bioTree{1}.imageSize(1);
ySize=bioTree{1}.imageSize(2);
pixelList=pixelIdxListOut{i};
mask=false(xSize,ySize,size(pixelList,2));
for iTrace=1:size(pixelList,2)
    tempMask=false(xSize,ySize);
    tempMask(pixelList{iTrace})=true;
    mask(:,:,iTrace)= tempMask;
end
end