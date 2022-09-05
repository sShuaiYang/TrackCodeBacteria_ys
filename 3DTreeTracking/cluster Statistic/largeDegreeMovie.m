function largeDegreeMovie(bioTree,clusterTree)
% 亮点程序：直接在格子图上做出了直线，可以直接用imwrite
finalClusterTree=clusterTree{end};
[bioTree,branchList,~,~]=myBiograph_new2(bioTree);
[bioTree,~,allList]=divisionFinder(bioTree,branchList);
distMatrix=full(finalClusterTree.distMatrix);
coreBranch=[3,10,2];
bacteriaList=finalClusterTree.bacteriaList;
branchList=branchList(coreBranch,:);
% coreBranch=[];
% bacteriaList=[];
% for i=1:size(distMatrix,2)
%     distInfo=distMatrix(:,i);
%     eachDegree=numel(distInfo(distInfo==1));
%     if eachDegree>800
%         coreBranch=[coreBranch;finalClusterTree.bacteriaList(i,5)];
%         bacteriaList=[bacteriaList;finalClusterTree.bacteriaList(i,:)];
%     end
% end
% coreBranch=sort(coreBranch);
% diffDist=diff(coreBranch);
% diffDist=[1;diffDist];
% coreBranch(diffDist==0)=[];
% branchList=branchList(coreBranch,:);
makeBranchDemo1(bioTree,branchList,allList,201,200,size(bioTree,2),bacteriaList,coreBranch,clusterTree);
end
function makeBranchDemo1(bioTree,branchList,allList,startFrame,stepFrame,endFrame,bacteriaList,coreBranch,clusterTree)
dirFile=uigetdir();
glueColor=[0.2,0.2,0.2];
branchColor=zeros(size(branchList,1),3);
nodeOrRoot=branchList(:,3);
coreBranch1=numel(nodeOrRoot(nodeOrRoot==1));
branchColor(1:coreBranch1,:)=colormap(hsv(coreBranch1));
branchColor(coreBranch1+1:end,1)=1;
% branchColor=zeros(size(branchList,1),3);
% branchColor(2,:)=[1,0,0];
% branchColor(5,:)=[0,1,0];
% branchColor(8,:)=[0,0,1];
% branchColor(11,:)=[1,1,0];
% branchColor(16,:)=[1,0,1];
% branchColor(18,:)=[0,1,1];
% branchColor(22,:)=[0.5,1,0.7];
% dirImages=uigetdir();
% dirSave=uigetdir();
dirImages=strcat(dirFile,'\tiff2matlab');
dirSave=strcat(dirFile,'\largeClusterDemo2');
mkdir(dirSave)
cd(dirImages);
nameList=dir(dirImages);
vars = whos('-file', nameList(3).name);
stackSize=vars.size(3);
stackNum=fix(startFrame/(stackSize));
if stackNum~=startFrame/(stackSize)
    fileName=nameList(stackNum+3).name;
else
    fileName=nameList(stackNum+2).name;
end
imageStack=loadImageStack(fileName);
% imageStack=imageStack(xCrop,yCrop,:);
iImage=1;
branchColor1(1,:)=[1,1,0];
branchColor1(2,:)=[1,0,1];
branchColor1(3,:)=[0,1,1];
for iframe=startFrame:stepFrame:endFrame;
    overLayer=2;
    disp(iframe);
    stackNumNext=fix(iframe/(stackSize));
    if stackNumNext==iframe/(stackSize)
        stackNumNext=stackNumNext-1;
    end
    if stackNumNext~=stackNum
        clear imageStack;
        fileName=nameList(stackNumNext+3).name;
        imageStack=loadImageStack(fileName);
%         imageStack=imageStack(xCrop,yCrop,:);
        stackNum=stackNumNext;
    end
%     demoImage=imoverlay(imageStack(:,:,iframe-stackNum*stackSize),glueMask(:,:,iImage),glueColor);
    demoImage=imoverlay(imageStack(:,:,iframe-stackNum*stackSize),[],glueColor);
    iImage=iImage+1;
    linkMaskPre=false(bioTree{1}.imageSize);
    for iBranch=1:size(branchList,1)
        [branchMask,isMask,linkMask]=getCoreBranchMask(bioTree,branchList,iBranch,iframe,bacteriaList(bacteriaList(:,5)==coreBranch(iBranch),:),clusterTree{iframe});
        if isMask==true
%             linkMaskPre=linkMask | linkMaskPre;
            if overLayer==1
                demoImage=imoverlay(imageStack(:,:,iframe-stackNum*stackSize),branchMask,branchColor(iBranch,:));
                overLayer=2;
            end
            if overLayer==2
                colorEach=branchColor(iBranch,:);
                demoImage=imoverlay(demoImage,linkMask,branchColor1(iBranch,:));
                demoImage1=demoImage(:,:,1);
                demoImage1(branchMask)=255*colorEach(1);
                demoImage2=demoImage(:,:,2);
                demoImage2(branchMask)=255*colorEach(2);
                demoImage3=demoImage(:,:,3);
                demoImage3(branchMask)=255*colorEach(3);
                demoImage(:,:,2)=demoImage2;
                demoImage(:,:,1)=demoImage1;
                demoImage(:,:,3)=demoImage3;
%                 demoImage=imoverlay(demoImage,branchMask,branchColor(iBranch,:));
%                 demoImage1=demoImage(:,:,1);
%                 demoImage1(linkMaskPre)=255;
%                 demoImage2=demoImage(:,:,2);
%                 demoImage2(linkMaskPre)=255;
%                 demoImage(:,:,2)=demoImage2;
%                 demoImage(:,:,1)=demoImage1;
            end
        end
    end
    [rootMask,isRootMask,leafMask,isLeafMask,divisionMask,isDivisionMask]=getEventMask(bioTree,allList,iframe,stepFrame);
%     if isRootMask==true
%         demoImage=imoverlay(demoImage,rootMask,rootColor);
%     end
%     if isLeafMask==true
%         demoImage=imoverlay(demoImage,leafMask,leafColor);
%     end
%     if isDivisionMask==true
%        demoImage=imoverlay(demoImage,divisionMask,divisionColor);
%     end
%     imshow(demoImage);
%     plotAllEventinList(iframe,allList,rootColor,leafColor,divisionColor,hyperNodeColor);
    saveFile2=strcat(dirSave,'\',num2str(iframe),'.tif');
%     saveas(h,saveFile2,'tif');
    imwrite(demoImage,saveFile2);
% imwrite(imageStack(:,:,iframe-stackNum*stackSize),saveFile2);
end
end
function plotAllEventinList(iframe,allList,rootColor,leafColor,divisionColor,hyperNodeColor)
allRoot=allList.allRoot;
allLeaf=allList.allLeaf;
allNode=allList.allNode;
for iroot=1:size(allRoot,1)
    if allRoot(iroot,1)>iframe
        break;
    end
    x1=allRoot(iroot,4);
    y1=allRoot(iroot,5);
    plot(x1,y1,'MarkerFaceColor',rootColor,'MarkerSize',2,'Marker','o','MarkerEdgeColor','none','LineStyle','none');hold on;
end
for iLeaf=1:size(allLeaf,1)
    if allLeaf(iLeaf,1)>iframe
        break;
    end
    x1=allLeaf(iLeaf,4);
    y1=allLeaf(iLeaf,5);
    plot(x1,y1,'MarkerFaceColor',leafColor,'MarkerSize',2,'Marker','o','MarkerEdgeColor','none','LineStyle','none');hold on;
end
for iNode=1:size(allNode,1)
    if allNode(iNode,1)>iframe
        break;
    end 
    if allNode(iNode,4)==true
        x1=allNode(iNode,7);
        y1=allNode(iNode,8);
        plot(x1,y1,'MarkerFaceColor',divisionColor,'MarkerSize',2,'Marker','o','MarkerEdgeColor','none','LineStyle','none');hold on;
    end
    if allNode(iNode,5)==true
        x1=allNode(iNode,7);
        y1=allNode(iNode,8);
        plot(x1,y1,'MarkerFaceColor',hyperNodeColor,'MarkerSize',2,'Marker','o','MarkerEdgeColor','none','LineStyle','none');hold on;
    end
end
end
function imageStack=loadImageStack(fileName)
tempImage=load(fileName);
imageStack=tempImage.imageStack;
end
function [branchMask,isMask,linkMask]=getCoreBranchMask(bioTree,branchList,branchIndex,frame,bacteriaList,distInfo)
xSize=bioTree{1}.imageSize(1);
ySize=bioTree{1}.imageSize(2);
imageSize=[xSize,ySize];
branchMask=false(xSize,ySize);
linkMask=false(xSize,ySize);
isMask=false;
branchInfo=branchList(branchIndex,1:3);
if branchInfo(3)==1
    allRoot=bioTree{branchInfo(1)}.node{branchInfo(2)}.allRoot;
    allNode=bioTree{branchInfo(1)}.node{branchInfo(2)}.allNode;
    for iNode=1:size(allNode,1)
        nodeInfo=allNode(iNode,:);
        for iOut=1:size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out,2)
            bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.needFill=0;
        end
    end
    for i=1:size(bacteriaList,1)
        bacteriaInfo=bacteriaList(i,:);
        while bioTree{bacteriaInfo(1)}.node{bacteriaInfo(2)}.In{1}.isNode==1
            nodeInfo=bacteriaInfo(1:3);
            bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.needFill=1;
            bacteriaInfo=bioTree{bacteriaInfo(1)}.node{bacteriaInfo(2)}.In{1}.nodeInfo;
        end
        nodeInfo=bacteriaInfo(1:3);
        bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.needFill=1;
        rootInfo=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{1}.rootInfo;
        iRoot=1;
        while size(allRoot,1)~=1
            if isequal(rootInfo,allRoot(iRoot,:))
               iRoot=iRoot+1;
            else
                allRoot(iRoot,:)=[];
                iRoot=iRoot-1;
            end
        end
        if ~isequal(rootInfo,allRoot)
            allRoot=rootInfo;
        end
    end        
end
if branchInfo(3)==0
    allRoot=bioTree{branchInfo(1)}.root{branchInfo(2)}.allRoot;
    allNode=[];
end
for iroot=1:size(allRoot,1)
    rootInfo=allRoot(iroot,1:2);
    if rootInfo(1)>frame
        continue;
    end
    if rootInfo(1)<=frame
        if bioTree{rootInfo(1)}.root{rootInfo(2)}.is2Node==1
            nextnodeInfo=bioTree{rootInfo(1)}.root{rootInfo(2)}.nodeInfo;
            if nextnodeInfo(1)<=frame
                continue;
            end
            if nextnodeInfo(1)>frame
                pixelIdxList=bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.pixelIdxList{frame-rootInfo(1)+1};
                pixelIdxList=pixelTransFilled(pixelIdxList,imageSize);
                branchMask(pixelIdxList)=1;
                isMask=true;
                eachBacteriaInfo=[rootInfo(1),rootInfo(2),0,frame-rootInfo(1)+1];
                linkMask=getLinkLine(xSize,eachBacteriaInfo,distInfo,linkMask,bioTree);
            end
        end
        if bioTree{rootInfo(1)}.root{rootInfo(2)}.is2Node==0
            leafInfo=bioTree{rootInfo(1)}.root{rootInfo(2)}.leafInfo;
            if leafInfo(1)<frame
                continue;
            end
            if leafInfo(1)>=frame
                pixelIdxList=bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.pixelIdxList{frame-rootInfo(1)+1};
                pixelIdxList=pixelTransFilled(pixelIdxList,imageSize);
                branchMask(pixelIdxList)=1;
                isMask=true;
                eachBacteriaInfo=[rootInfo(1),rootInfo(2),0,frame-rootInfo(1)+1];
                linkMask=getLinkLine(xSize,eachBacteriaInfo,distInfo,linkMask,bioTree);
            end
        end
    end
end
for iNode=1:size(allNode,1)
    nodeInfo=allNode(iNode,1:2);
    if nodeInfo(1)>frame
        continue;
    end
    if nodeInfo(1)<=frame
        for iOut=1:size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out,2)
            if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.is2Node==true
                nextnodeInfo=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.nodeInfo;
                if nextnodeInfo(1)<=frame
                    continue;
                end
                if nextnodeInfo(1)>frame
                    pixelIdxList=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.traceInfo.pixelIdxList{frame-nodeInfo(1)+1};
                    if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.needFill==1
                        pixelIdxList=pixelTransFilled(pixelIdxList,imageSize);
                        eachBacteriaInfo=[nodeInfo(1),nodeInfo(2),iOut,frame-nodeInfo(1)+1];
                        linkMask=getLinkLine(xSize,eachBacteriaInfo,distInfo,linkMask,bioTree);
                    end
                    branchMask(pixelIdxList)=1;
                    isMask=true;
                end
            end
            if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.is2Node==false
                leafInfo=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.leafInfo;
                if leafInfo(1)<frame
                    continue;
                end
                if leafInfo(1)>=frame
                    pixelIdxList=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.traceInfo.pixelIdxList{frame-nodeInfo(1)+1};
                    if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.needFill==1
                        pixelIdxList=pixelTransFilled(pixelIdxList,imageSize);
                        eachBacteriaInfo=[nodeInfo(1),nodeInfo(2),iOut,frame-nodeInfo(1)+1];
                        linkMask=getLinkLine(xSize,eachBacteriaInfo,distInfo,linkMask,bioTree);
                    end
                    branchMask(pixelIdxList)=1;
                    isMask=true;
                end
            end
        end
    end
end
% if isMask==true
%     branchMask=imfill(branchMask,'holes');
% end
end
function pixelIdxList=pixelTransFilled(pixelIdxList,imageSize)
[xyMin,BWImageGain]=idx2Xy(pixelIdxList,imageSize);
BWImageGain=imfill(BWImageGain,'holes');
pixelIdxList=xy2Idx(xyMin,BWImageGain,imageSize);
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
end
function pixelIdxListinLine=pixelInLine(centroid1,centroid2,xSize)
if abs((centroid2(1)-centroid1(1)))>abs((centroid2(2)-centroid1(2)))
    if fix(centroid1(1))<=fix(centroid2(1))
        xinLine=fix(centroid1(1)):1:fix(centroid2(1));
    else
        xinLine=fix(centroid1(1)):-1:fix(centroid2(1));
    end
    if fix(centroid2(1)-centroid1(1))~=0
        k=((centroid2(2)-centroid1(2))/(centroid2(1)-centroid1(1)));
        b=centroid1(2)-centroid1(1)*k;
        yinLine=fix(xinLine.*k+b);
    else
        yinLine=fix(centroid2(2));
    end
    pixelIdxListinLine=xinLine.*xSize+yinLine;
    return;
end

if abs((centroid2(1)-centroid1(1)))<=abs((centroid2(2)-centroid1(2)))
    if fix(centroid1(2))<=fix(centroid2(2))
        yinLine=fix(centroid1(2)):1:fix(centroid2(2));
    else
        yinLine=fix(centroid1(2)):-1:fix(centroid2(2));
    end
    if fix(centroid2(1)-centroid1(1))~=0
        k=((centroid2(2)-centroid1(2))/(centroid2(1)-centroid1(1)));
        b=centroid1(2)-centroid1(1)*k;
        xinLine=fix((yinLine-b)./k);
    else
        xinLine=fix(centroid1(1));
    end
    pixelIdxListinLine=xinLine.*xSize+yinLine;
    return;
end
end
function linkMask=getLinkLine(xSize,bacteriaInfo,distInfo,linkMask,bioTree)
if bacteriaInfo(3)==0
    aimCentroid=bioTree{bacteriaInfo(1)}.root{bacteriaInfo(2)}.traceInfo.measurment{bacteriaInfo(4)}(1).Centroid;
else
    aimCentroid=bioTree{bacteriaInfo(1)}.node{bacteriaInfo(2)}.Out{bacteriaInfo(3)}.traceInfo.measurment{bacteriaInfo(4)}(1).Centroid;
end
for iBac=1:size(distInfo.bacteriaList,1)
    if isequal(bacteriaInfo,distInfo.bacteriaList(iBac,1:4))
        linkInfo=distInfo.distMatrix(iBac,:);
    end
end
for iBac=1:numel(linkInfo)
    if linkInfo(iBac)==1
        linkBacteria=distInfo.bacteriaList(iBac,:);
        if linkBacteria(3)==0
            linkCentroid=bioTree{linkBacteria(1)}.root{linkBacteria(2)}.traceInfo.measurment{linkBacteria(4)}(1).Centroid;
        else
            linkCentroid=bioTree{linkBacteria(1)}.node{linkBacteria(2)}.Out{linkBacteria(3)}.traceInfo.measurment{linkBacteria(4)}(1).Centroid;
        end
        pixelIdxListinLine=pixelInLine(aimCentroid,linkCentroid,xSize);
        linkMask(pixelIdxListinLine)=1;
    end
end
end