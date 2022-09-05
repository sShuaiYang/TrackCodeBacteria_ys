in=[];
out=[];
for i=1:numel(afterDivide)
[regionNum,xyMin,regionImage]=findRegionNum(afterDivide{i},pictureSize);
in{i}.BWImage=regionImage{1};
in{i}.xyMin=xyMin;
end
[regionNum,xyMin,regionImage]=findRegionNum(pixelIdxListOut{1}{4},pictureSize);
for i=1:numel(regionImage)
out{i}.xyMin=xyMin;
out{i}.BWImage=regionImage{i};
end

eachInfo=afterDivide;
maskImage=pixelIdxListOut{1,1}{1};
tic;[afterDivide,canDivideorNot]=basicDivideSolution(eachInfo,maskImage,pictureSize);toc
[regionNum,xyMin,regionImage]=findRegionNum(afterDivide{1},pictureSize);
[regionNum,xyMin,regionImage]=findRegionNum(afterDivide{2},pictureSize);
[regionNum,xyMin,regionImage]=findRegionNum(maskImage,pictureSize);