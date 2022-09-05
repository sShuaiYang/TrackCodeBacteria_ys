function resultMap=traceStatisticMap(bioTree,bacteriaFrameInfo,step)
picSize=bioTree{1}.imageSize;
resultMap=zeros(picSize);
for iframe=1:step:numel(bioTree)
    bacteriaInfo=bacteriaFrameInfo{iframe}.bacteriaInfo;
    for ibac=1:size(bacteriaInfo,1)
        iBacInfo=bacteriaInfo(ibac,:);
        if iBacInfo(3)==0
            pixelIdxList=bioTree{iBacInfo(1)}.root{iBacInfo(2)}.traceInfo.pixelIdxList{iBacInfo(4)};
        else
            pixelIdxList=bioTree{iBacInfo(1)}.node{iBacInfo(2)}.Out{iBacInfo(3)}.traceInfo.pixelIdxList{iBacInfo(4)};
        end
        [xyMin,BWImage]=idx2Xy(pixelIdxList,picSize);
        pixelIdxList=xy2Idx(xyMin,BWImage,picSize);
        resultMap(pixelIdxList)=resultMap(pixelIdxList)+1;
    end
end
resultMap=resultMap./max(max(resultMap));
resultMap=im2uint8(resultMap);
end
function [xyMin,BWImageGain]=idx2Xy(pixelIdxList,pictureSize)
% this function is used to convert linearIdx to a small BWImage
xSize=pictureSize(1);
yresult=ceil(pixelIdxList/xSize);
xresult=round(pixelIdxList-(yresult-1)*xSize);
xMin=min(xresult);
xMax=max(xresult);
yMin=min(yresult);
yMax=max(yresult);
yresult2=yresult-yMin+1;
xresult2=xresult-xMin+1;
BWImage=false(xMax-xMin+1,yMax-yMin+1);
pixelIdxList2=xresult2+(yresult2-1)*(xMax-xMin+1);
BWImage(pixelIdxList2)=true;
BWImageGain=false(xMax-xMin+3,yMax-yMin+3);
xyMin=[xMin,yMin];
BWImageGain(2:end-1,2:end-1)=BWImage;
BWImageGain=imfill(BWImageGain,'holes');
end
function pixelIdxList=xy2Idx(xyMin,BWImageGain,pictureSize)
%this function is used to convert a small BWImage to its original pixelIdx
BWImage=BWImageGain(2:end-1,2:end-1);
% BWImage=bwmorph(BWImage,'remove');
pixelIdxListOri=find(BWImage==1);
smallxSize=size(BWImage,1);
yresult=ceil(pixelIdxListOri/smallxSize);
xresult=pixelIdxListOri-(yresult-1)*smallxSize;
xresult2=xresult+xyMin(1)-1;
yresult2=yresult+xyMin(2)-1;
xSize=pictureSize(1);
pixelIdxList=xresult2+(yresult2-1)*xSize;
end