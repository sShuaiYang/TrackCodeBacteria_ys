function finalImage=oneGeneDiffusionImage(bioTree,clusterTree,resultAll)
branchNum=13;
finalImage=zeros([bioTree{1}.imageSize,3],'uint8');
image1=finalImage(:,:,1);
image2=finalImage(:,:,2);
image3=finalImage(:,:,3);
iframe=4200;
nodeInfo=clusterTree(iframe).nodeInfo;
finalMatrix=resultAll{iframe/200}.finalMatrix;
branchMatrix=clusterTree(iframe).branchMatrix;
finalMatrix=finalMatrix/160000;
finalMatrix=finalMatrix+branchMatrix;
iRow=finalMatrix(:,branchNum);
logRow=log10(iRow);
logRow(logRow==-Inf)=min(ceil(logRow(logRow~=-Inf)))-1;
logRow=logRow-min(logRow);
colorGet=colormap(jet(10*max(ceil(logRow))));
for i=1:size(nodeInfo,1)
    if nodeInfo(i,3)==-1        
        pixelIdxList=bioTree{nodeInfo(i,1)}.leavies{nodeInfo(i,2)}.leafPixelDetail;
    end
    if nodeInfo(i,3)==0
        pixelIdxList=bioTree{nodeInfo(i,1)}.root{nodeInfo(i,2)}.traceInfo.pixelIdxList{nodeInfo(i,4)};
    end
    if nodeInfo(i,3)>=1
        pixelIdxList=bioTree{nodeInfo(i,1)}.node{nodeInfo(i,2)}.Out{nodeInfo(i,3)}.traceInfo.pixelIdxList{nodeInfo(i,4)};
    end
    [xyMin,BWImageGain]=idx2Xy(pixelIdxList,bioTree{1}.imageSize);
    pixelIdxList=xy2Idx(xyMin,BWImageGain,bioTree{1}.imageSize);
    if logRow(i)==0
        color=[1,1,1];
    else
        color=colorGet(ceil(10*logRow(i)),:);
    end
    image1(pixelIdxList)=color(1)*255;
    image2(pixelIdxList)=color(2)*255;
    image3(pixelIdxList)=color(3)*255;
end
finalImage=cat(3,image1,image2,image3);
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
xyMin=[xMin,yMin];
end
function pixelIdxList=xy2Idx(xyMin,BWImageGain,pictureSize)
%this function is used to convert a small BWImage to its original pixelIdx
BWImage=BWImageGain(2:end-1,2:end-1);
pixelIdxListOri=find(BWImage==1);
smallxSize=size(BWImage,1);
yresult=ceil(pixelIdxListOri/smallxSize);
xresult=pixelIdxListOri-(yresult-1)*smallxSize;
xresult2=xresult+xyMin(1)-1;
yresult2=yresult+xyMin(2)-1;
xSize=pictureSize(1);
pixelIdxList=xresult2+(yresult2-1)*xSize;
end
