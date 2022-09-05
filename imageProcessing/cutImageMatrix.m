function bwImageStack=cutImageMatrix(CC,nObjective, maskIamge)
xSize=CC.ImageSize(1);
ySize=CC.ImageSize(1);
nPixelIdxList=CC.PixelIdxList(1,nObjective);
nPixelIdx=nPixelIdxList{1};
startFrame=fix(nPixelIdx(1)/(xSize*ySize))+1;
if nPixelIdx(length(nPixelIdx))/(xSize*ySize)~=fix(nPixelIdx(length(nPixelIdx))/(xSize*ySize));
endFrame=fix(nPixelIdx(length(nPixelIdx))/(xSize*ySize))+1;
else
endFrame=fix(nPixelIdx(length(nPixelIdx))/(xSize*ySize));   
end
maskStack1=maskIamge(:,:,(startFrame:endFrame));
maskIamge(CC.PixelIdxList{nObjective}) = 0;
maskStack2=maskIamge(:,:,(startFrame:endFrame));
bwImageStack=(maskStack1)&(~maskStack2);
end