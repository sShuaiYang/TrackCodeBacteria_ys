function makeDegreeDemo(bioTree,dirFile,dirFileTiff2mat,clusterTree,type)
% dirFinalTree=strcat(dirFile,'\bioTreeResult\allTree\finalTree');
% measureResult=dir(dirFinalTree);
% measureResultNum=numel(measureResult)-2;
% for i=1:measureResultNum
%     load(strcat(dirFinalTree,'\',strcat('_',num2str(i),'_')));
%     if i==1
%         bioTree=bioTreeStack;
%     else
%         bioTree(end+1:end+size(bioTreeStack,2))=bioTreeStack;
%     end
% end
% if strcmp(type,'wspf')
%     load(strcat(dirFile,'\forbiddenImage.mat'))
%     bioTree=bioTreeCorrection(bioTree,forbiddenImage);
% else
%     bioTree=bioTreeCorrection(bioTree);
% end
endMatrix=clusterTree{end}.distMatrix;
degree=[];
for i=1:size(endMatrix,1)
    iLine=endMatrix(i,:);
    degree(i)=numel(iLine(iLine==1));
end
% maxDegree=max(degree);
maxDegree=300;
colorBuild=zeros(maxDegree,3);
colorBuild(:,1)=0.5;
colorBuild(1:150,:)=colormap(jet(150));
dirImages=strcat(dirFileTiff2mat,'\tiff2matlab');
dirSave=strcat(dirFile,'\degreeDemo');
mkdir(dirSave);
cd(dirImages);
nameList=dir(dirImages);
vars = whos('-file', nameList(3).name);
stackSize=vars.size(3);
startFrame=401;
step=200;
endFrame=numel(clusterTree);
stackNum=fix(startFrame/(stackSize));
if stackNum~=startFrame/(stackSize)
    fileName=nameList(stackNum+3).name;
else
    fileName=nameList(stackNum+2).name;
end
imageStack=loadImageStack(fileName);
for iframe=startFrame:step:endFrame;
    disp(iframe);
    stackNumNext=fix(iframe/(stackSize));
    if stackNumNext==iframe/(stackSize)
        stackNumNext=stackNumNext-1;
    end
    if stackNumNext~=stackNum
        clear imageStack;
        fileName=nameList(stackNumNext+3).name;
        imageStack=loadImageStack(fileName);
        stackNum=stackNumNext;
    end
    demoImage=imoverlay(imageStack(:,:,iframe-stackNum*stackSize),colorBuild,clusterTree{iframe},bioTree,iframe);
    saveFile2=strcat(dirSave,'\',num2str(iframe),'.tif');
    imwrite(demoImage,saveFile2);
end
end
function demoImage=imoverlay(imageStack,colorBuild,clusterDetail,bioTree,iframe)
distMatrix=full(clusterDetail.distMatrix);
bacteriaList=clusterDetail.bacteriaList;
orderNum=(1:size(bacteriaList,1))';
orderNum=orderNum((bacteriaList(:,1)+bacteriaList(:,4)-1)==iframe,:);
properList=bacteriaList((bacteriaList(:,1)+bacteriaList(:,4)-1)==iframe,:);
demoImage=cat(3,imageStack,imageStack,imageStack);
for i=1:size(properList,1)
    iList=properList(i,:);
    if iList(3)==0
        pixelIdxList=bioTree{iList(1)}.root{iList(2)}.traceInfo.pixelIdxList{iList(4)};
    end
    if iList(3)~=0
        pixelIdxList=bioTree{iList(1)}.node{iList(2)}.Out{iList(3)}.traceInfo.pixelIdxList{iList(4)};
    end
    [xyMin,BWImageGain]=idx2Xy(pixelIdxList,bioTree{1}.imageSize);
    pixelIdxList=xy2Idx(xyMin,BWImageGain,bioTree{1}.imageSize);
    distMatrixLine=distMatrix(orderNum(i),:);
    iDegree=numel(distMatrixLine(distMatrixLine==1));
    if iDegree>=1
        imageStack1=demoImage(:,:,1);
        imageStack2=demoImage(:,:,2);
        imageStack3=demoImage(:,:,3);
        imageStack1(pixelIdxList)=255*colorBuild(iDegree,1);
        imageStack2(pixelIdxList)=255*colorBuild(iDegree,2);
        imageStack3(pixelIdxList)=255*colorBuild(iDegree,3);
        demoImage=cat(3,imageStack1,imageStack2,imageStack3);
    end
end
end
function [xyMin,BWImageGain]=idx2Xy(pixelIdxList,pictureSize)
% this function is used to convert linearIdx to a small BWImage
xSize=pictureSize(1);
yresult=ceil(pixelIdxList/xSize);
xresult=pixelIdxList-(yresult-1)*xSize;
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
pixelIdxListOri=find(BWImage==1);
if size(pixelIdxListOri,2)~=1
    pixelIdxListOri=pixelIdxListOri';
end
smallxSize=size(BWImage,1);
yresult=ceil(pixelIdxListOri/smallxSize);
xresult=pixelIdxListOri-(yresult-1)*smallxSize;
xresult2=xresult+xyMin(1)-1;
yresult2=yresult+xyMin(2)-1;
xSize=pictureSize(1);
pixelIdxList=xresult2+(yresult2-1)*xSize;
end
function imageStack=loadImageStack(fileName)
tempImage=load(fileName);
imageStack=tempImage.imageStack;
end