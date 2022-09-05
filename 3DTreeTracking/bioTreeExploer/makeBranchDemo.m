function makeBranchDemo(bioTree,branchList,allList,startFrame,stepFrame,endFrame,dirFile,bacteriaFrameInfo,glueMask)
% function makeBranchDemo(bioTree,branchList,glueMask,allList,startFrame,stepFrame,endFrame)
myBranch=[3,4,6];
rootColor=[1,0,0];
divisionColor=[1,1,0];
leafColor=[0,1,0];
hyperNodeColor=[0,1,1];
glueColor=[0.5,0,0.5];
branchColor=zeros(size(branchList,1),3);
nodeOrRoot=branchList(:,3);
coreBranch=numel(nodeOrRoot(nodeOrRoot==1));
coreBranch=8;
branchColor(1:coreBranch,:)=colormap(hsv(coreBranch));
branchColor(coreBranch+1:end,1)=1;
dirImages=strcat(dirFile,'\tiff2matlab');
dirSave=strcat(dirFile,'\demo3');
mkdir(dirSave);
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
iImage=1;
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
        stackNum=stackNumNext;
    end
    demoImage=imoverlay(imageStack(:,:,iframe-stackNum*stackSize),glueMask(:,:,iImage),glueColor);
%     demoImage=imoverlay(imageStack(1:end,:,iframe-stackNum*stackSize),[],glueColor);
    iImage=iImage+1;
    if nargin==8
        for iBranch=1:size(branchList,1)
            [branchMask,isMask]=getBranchMask(bioTree,branchList,iBranch,iframe);
            branchMask=bwmorph(branchMask,'remove');
            if isMask==true
                if overLayer==1
                    demoImage=imoverlay(imageStack(:,:,iframe-stackNum*stackSize),branchMask,branchColor(iBranch,:));
                    overLayer=2;
                end
                if overLayer==2
                    demoImage=imoverlay(demoImage,branchMask,branchColor(iBranch,:));
                end
            end
        end
    else if nargin==9
            properFrame=bacteriaFrameInfo{iframe};
            for i=1:size(properFrame.bacteriaInfo)
                if ismember(properFrame.bacteriaInfo(i,5),myBranch)
                    if properFrame.bacteriaInfo(i,3)==0
                        rootInfo=properFrame.bacteriaInfo(i,:);
                        pixelIdxList=bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.pixelIdxList{rootInfo(4)};
                        iBranch=rootInfo(5);
                    else
                        nodeInfo=properFrame.bacteriaInfo(i,:);
                        pixelIdxList=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.pixelIdxList{nodeInfo(4)};
                        iBranch=nodeInfo(5);
                    end
                    [xyMin,BWImage]=idx2Xy(pixelIdxList,bioTree{1}.imageSize);
                    pixelIdxList=xy2Idx(xyMin,BWImage,bioTree{1}.imageSize); %add by jzy fill holes
                    for j=1:3
                        image=demoImage(:,:,j);
                        image(pixelIdxList)=255*branchColor(iBranch,j);
                        demoImage(:,:,j)=image;
                    end
                end
            end
        end
    end
    [rootMask,isRootMask,leafMask,isLeafMask,divisionMask,isDivisionMask]=getEventMask(bioTree,allList,iframe,stepFrame);
    if isRootMask==true
        demoImage=imoverlay(demoImage,rootMask,rootColor);
    end
    if isLeafMask==true
        demoImage=imoverlay(demoImage,leafMask,leafColor);
    end
    if isDivisionMask==true
       demoImage=imoverlay(demoImage,divisionMask,divisionColor);
    end
%     imshow(demoImage);
%     plotAllEventinList(iframe,allList,rootColor,leafColor,divisionColor,hyperNodeColor);
    saveFile2=strcat(dirSave,'\',num2str(iframe),'.tif');
%     saveas(h,saveFile2,'tif');
    demoImage=demoImage(:,end:-1:1,:);
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