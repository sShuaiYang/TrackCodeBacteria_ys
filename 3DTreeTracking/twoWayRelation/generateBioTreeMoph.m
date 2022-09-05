function generateBioTreeMoph(bioTree,bacteriaFrameInfo)
orderBranch=7;
for i=1:numel(bacteriaFrameInfo)
    bacteriaFrameInfo{i}.bacteriaInfo(~ismember(bacteriaFrameInfo{i}.bacteriaInfo(:,5),orderBranch),:)=[];
end
pixelIdxList=[];
for i=1
    for j=1:size(bacteriaFrameInfo{i}.bacteriaInfo,1)
        bacInfo=bacteriaFrameInfo{i}.bacteriaInfo(j,:);
        if bacInfo(3)==0
            pixelIdxList=[pixelIdxList;bioTree{bacInfo(1)}.root{bacInfo(2)}.traceInfo.pixelIdxList{bacInfo(4)}];
        end
        if bacInfo(3)~=0
            pixelIdxList=[pixelIdxList;bioTree{bacInfo(1)}.node{bacInfo(2)}.Out{bacInfo(3)}.traceInfo.pixelIdxList{bacInfo(4)}];
        end
    end
    pixelIdxList=unique(pixelIdxList);
end
regionInfo=getSumRegion(pixelIdxList,bioTree{1}.imageSize);
color=[1,0,0];
dirFile=strcat('C:\Users\jzy\Desktop\test Morph\ori',num2str(orderBranch));
mkdir(dirFile)
for iframe=1:20:numel(bacteriaFrameInfo)
% for iframe=1
    pixelIdxList=[];
    finalImage=zeros(regionInfo.xMax-regionInfo.xMin+3,regionInfo.yMax-regionInfo.yMin+3,3);
    for j=1:size(bacteriaFrameInfo{iframe}.bacteriaInfo,1)
        bacInfo=bacteriaFrameInfo{iframe}.bacteriaInfo(j,:);
        if bacInfo(3)==0
            pixelIdxList=[pixelIdxList;bioTree{bacInfo(1)}.root{bacInfo(2)}.traceInfo.pixelIdxList{bacInfo(4)}];
        end
        if bacInfo(3)~=0
            pixelIdxList=[pixelIdxList;bioTree{bacInfo(1)}.node{bacInfo(2)}.Out{bacInfo(3)}.traceInfo.pixelIdxList{bacInfo(4)}];
        end
    end
    pixelIdxList=unique(pixelIdxList);
%     imageSize=[2160,2560];
%     p=false(imageSize);
%     p(pixelIdxList)=1;
%     imshow(p)
    [~,BWImageGain]=idx2Xy(pixelIdxList,bioTree{1}.imageSize,regionInfo);
    image1=finalImage(:,:,1);
    image2=finalImage(:,:,2);
    image3=finalImage(:,:,3);
    image1(BWImageGain)=255*color(1);
    image2(BWImageGain)=255*color(2);
    image3(BWImageGain)=255*color(3);
    finalImage=cat(3,image1,image2,image3);
    finalImage=im2uint8(finalImage);
    imwrite(finalImage,strcat(dirFile,'\',num2str(iframe),'.tif'))
end
end
function regionInfo=getSumRegion(pixelIdxList,pictureSize)
% get the min region for input and output images
xSize=pictureSize(1);
yresult=ceil(pixelIdxList/xSize);
xresult=pixelIdxList-(yresult-1)*xSize;
% regionInfo.xMin=min(xresult);
% regionInfo.xMax=max(xresult);
% regionInfo.yMin=min(yresult);
% regionInfo.yMax=max(yresult);
regionInfo.xMin=1000;
regionInfo.xMax=1900;
regionInfo.yMin=250;
regionInfo.yMax=1100;
end
function [xyMin,BWImageGain]=idx2Xy(pixelIdxList,pictureSize,regionInfo)
% this function is used to convert linearIdx to a small BWImage
xSize=pictureSize(1);
yresult=ceil(pixelIdxList/xSize);
xresult=round(pixelIdxList-(yresult-1)*xSize);
xMin=regionInfo.xMin;
xMax=regionInfo.xMax;
yMin=regionInfo.yMin;
yMax=regionInfo.yMax;
yresult2=yresult-yMin+1;
properY=yresult2>=1 & yresult2<=yMax-yMin+1;
xresult2=xresult-xMin+1;
properX=xresult2>=1 & xresult2<=xMax-xMin+1;
properAll=~properX | ~properY;
xresult2(properAll)=[];
yresult2(properAll)=[];
BWImage=false(xMax-xMin+1,yMax-yMin+1);
pixelIdxList2=xresult2+(yresult2-1)*(xMax-xMin+1);
BWImage(pixelIdxList2)=true;
BWImageGain=false(xMax-xMin+3,yMax-yMin+3);
xyMin=[xMin,yMin];
BWImageGain(2:end-1,2:end-1)=BWImage;
BWImageGain=imfill(BWImageGain,'holes');
BWImageGain=bwmorph(BWImageGain,'dilate');
end