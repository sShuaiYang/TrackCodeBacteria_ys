function imageStack=transProbabilityFigure(bioTree,bacteriaInfo,linkRate,coreBranch)
% 新network某个特定branch的概率分布
% bacteriaInfo是从clusterTree里得到的
% linkRate 就是某个branch的概率信息
imageStack1=uint8(zeros(bioTree{1}.imageSize));
imageStack2=uint8(zeros(bioTree{1}.imageSize));
imageStack3=uint8(zeros(bioTree{1}.imageSize));
linkRate=log(linkRate);
if ~all(linkRate==-inf)
    linkRate=linkRate+12;
%     linkRate=linkRate-min(linkRate(linkRate~=-inf));
    linkRate(linkRate==-inf)=0;
    linkRate=100/12*linkRate;
    colorAll=colormap(jet(130));
    colorAll(1:30,:)=[];
%     linkRate=100/max(linkRate)*linkRate;
%     colorAll=colormap(jet(ceil(max(linkRate))));
else
    linkRate(linkRate==-inf)=0;
end
for iBac=1:size(bacteriaInfo,1)
    iLine=bacteriaInfo(iBac,:);
    if iLine(3)>0
        pixelIdxList=bioTree{bacteriaInfo(iBac,1)}.node{bacteriaInfo(iBac,2)}.Out{bacteriaInfo(iBac,3)}.traceInfo.pixelIdxList{bacteriaInfo(iBac,4)};
    end
    if iLine(3)==0
        pixelIdxList=bioTree{bacteriaInfo(iBac,1)}.root{bacteriaInfo(iBac,2)}.traceInfo.pixelIdxList{bacteriaInfo(iBac,4)};
    end
    if iLine(3)==-1
        pixelIdxList=bioTree{bacteriaInfo(iBac,1)}.leavies{bacteriaInfo(iBac,2)}.leafPixelDetail;
    end
    [xyMin,BWImage]=idx2Xy(pixelIdxList,bioTree{1}.imageSize);
    pixelIdxList=xy2Idx(xyMin,BWImage,bioTree{1}.imageSize);
    if linkRate(iBac)==0
        BWImage=bwmorph(BWImage,'remove');
        pixelIdxList=xy2Idx(xyMin,BWImage,bioTree{1}.imageSize);
        imageStack1(pixelIdxList)=255;
        imageStack2(pixelIdxList)=255;
        imageStack3(pixelIdxList)=255;
    end
    if ceil(linkRate(iBac))~=0
        bacColor=colorAll(ceil(linkRate(iBac)),:);
        imageStack1(pixelIdxList)=bacColor(1)*255;
        imageStack2(pixelIdxList)=bacColor(2)*255;
        imageStack3(pixelIdxList)=bacColor(3)*255;
    else
        if iLine(5)==coreBranch
            imageStack1(pixelIdxList)=255;
            imageStack3(pixelIdxList)=255;
        else
            imageStack3(pixelIdxList)=255;
        end
    end
end
imageStack=cat(3,imageStack1,imageStack2,imageStack3);
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
pixelIdxList2=round(xresult2+(yresult2-1)*(xMax-xMin+1));
BWImage(pixelIdxList2)=true;
BWImageGain=false(xMax-xMin+3,yMax-yMin+3);
BWImageGain(2:end-1,2:end-1)=BWImage;
BWImageGain=imfill(BWImageGain,'holes');
% BWImageGain=bwmorph(BWImageGain,'dilate');
xyMin=[xMin,yMin];
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