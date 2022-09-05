function [newTrace,canDivideorNot,bioTree]=fullNodeType1Tracking_deep(bioTree,nodeInfo)%further track the N input and one Output node, return N new complete bacterial trace
imageSize=bioTree{1}.imageSize;
pixelIdxListIn=getInputMask(bioTree,nodeInfo);
for iIn=1:size(pixelIdxListIn,2)
    newTrace{iIn}.traceInfo.pixelIdxList=[];
    newTrace{iIn}.is2Node=[];
end
newRootInfo=0;
for iTrace=1:size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1}.traceInfo.pixelIdxList,2)
    [pixelIdxListNew,canDivideorNot,newRoot]=basicDivideSolution_deep(pixelIdxListIn,bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1}.traceInfo.pixelIdxList{iTrace},imageSize);
    if newRoot==1
        if numel(newRootInfo)==1 && newRootInfo==0
            newRootInfo=iTrace+nodeInfo(1)-1;
        else
            newRootInfo=[newRootInfo;iTrace+nodeInfo(1)-1];
        end
    end
    if canDivideorNot==1
        leafNum=size(bioTree{iTrace+nodeInfo(1)-2}.leavies,2);
        [newTrace,pixelIdxListIn,isBreak,error]=getNewTrace(newTrace,pixelIdxListNew,iTrace+nodeInfo(1)-2,leafNum);
        if error==1
            canDivideorNot=0;
            newTrace=[];
            return
        end
        if isBreak==true
            for iIn=1:size(newTrace,2)
                if  isempty(newTrace{iIn}.is2Node)
                    if iTrace<size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1}.traceInfo.pixelIdxList,2)
                        newTrace{iIn}.traceInfo.pixelIdxList=[newTrace{iIn}.traceInfo.pixelIdxList,bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1}.traceInfo.pixelIdxList(iTrace+1:end)];
                    end
                end
            end
            break;
        end
    else
        return
    end
end
testis2Node=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1}.is2Node;
if testis2Node==true
    testnodeInfo=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1}.nodeInfo;
    inNum=size(bioTree{testnodeInfo(1)}.node{testnodeInfo(2)}.In,2);
end
realInNode=0;
realInLeaf=0;
for iIn=1:size(newTrace,2)
%     if newTrace{iIn}.is2Node==0
%         if isempty(newTrace{iIn}.traceInfo.pixelIdxList)
%             pixelIdxListIn=getInputMask(bioTree,nodeInfo);
%             newTrace{iIn}.traceInfo.pixelIdxList{1}=pixelIdxListIn{1,iIn};
%         end
%     end
    if isempty(newTrace{iIn}.is2Node)
        newTrace{iIn}.is2Node=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1}.is2Node;
        if newTrace{iIn}.is2Node==true
            newTrace{iIn}.nodeInfo=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1}.nodeInfo;
            realInNode=realInNode+1;
            if realInNode~=1
                newTrace{iIn}.nodeInfo(1,3)=inNum+realInNode-1;
%                 newNodeInfo=newTrace{iIn}.nodeInfo;
            end
        end
        if newTrace{iIn}.is2Node==false
            newTrace{iIn}.leafInfo=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1}.leafInfo;
            if realInLeaf~=0
                newTrace{iIn}.leafInfo(1,2)=size(bioTree{newTrace{iIn}.leafInfo(1)}.leavies,2)+realInLeaf;
            end
            realInLeaf=realInLeaf+1;
        end
    end
end
for iRoot=size(newRootInfo,1)
    if size(newRootInfo,1)==1 && newRootInfo==0
        return
    end
    oriRootNum=size(bioTree{newRootInfo(iRoot)}.root,2);
    bioTree{1,newRootInfo(iRoot)}.root{1,oriRootNum+1}.is2Node=newTrace{end-size(newRootInfo,1)+iRoot}.is2Node;
    if bioTree{1,newRootInfo(iRoot)}.root{1,oriRootNum+1}.is2Node==1
        bioTree{1,newRootInfo(iRoot)}.root{1,oriRootNum+1}.leafInfo=[];
        nextNodeInfo=newTrace{end-size(newRootInfo,1)+iRoot}.nodeInfo;
        bioTree{1,newRootInfo(iRoot)}.root{1,oriRootNum+1}.nodeInfo=nextNodeInfo;
        bioTree{1,newRootInfo(iRoot)}.root{1,oriRootNum+1}.rootPixelDetail=newTrace{end-size(newRootInfo,1)+iRoot}.traceInfo.pixelIdxList{1};
        bioTree{1,newRootInfo(iRoot)}.root{1,oriRootNum+1}.traceInfo.pixelIdxList=newTrace{end-size(newRootInfo,1)+iRoot}.traceInfo.pixelIdxList;
        bioTree{nextNodeInfo(1)}.node{nextNodeInfo(2)}.In{nextNodeInfo(3)}.isNode=0;
        bioTree{nextNodeInfo(1)}.node{nextNodeInfo(2)}.In{nextNodeInfo(3)}.nodeInfo=[];
        bioTree{nextNodeInfo(1)}.node{nextNodeInfo(2)}.In{nextNodeInfo(3)}.rootInfo=[newRootInfo(iRoot),oriRootNum+1];
    end
    if bioTree{1,newRootInfo(iRoot)}.root{1,oriRootNum+1}.is2Node==0
        bioTree{1,newRootInfo(iRoot)}.root{1,oriRootNum+1}.rootInfo=[];
        nextLeafInfo=newTrace{end-size(newRootInfo,1)+iRoot}.leafInfo;
        bioTree{1,newRootInfo(iRoot)}.root{1,oriRootNum+1}.leafInfo=nextLeafInfo;
        bioTree{1,newRootInfo(iRoot)}.root{1,oriRootNum+1}.rootPixelDetail=newTrace{end-size(newRootInfo,1)+iRoot}.traceInfo.pixelIdxList{1};
        bioTree{1,newRootInfo(iRoot)}.root{1,oriRootNum+1}.traceInfo.pixelIdxList=newTrace{end-size(newRootInfo,1)+iRoot}.traceInfo.pixelIdxList;
        bioTree{nextLeafInfo(1)}.leavies{nextLeafInfo(2)}.is2Node=0;
        bioTree{nextLeafInfo(1)}.leavies{nextLeafInfo(2)}.nodeInfo=[];
        bioTree{nextLeafInfo(1)}.leavies{nextLeafInfo(2)}.rootInfo=[newRootInfo(iRoot),oriRootNum+1];
        bioTree{nextLeafInfo(1)}.leavies{nextLeafInfo(2)}.leaviesPixelDetail=newTrace{end-size(newRootInfo,1)+iRoot}.traceInfo.pixelIdxList{end};
    end
end
newTrace(end-size(newRootInfo,1)+1:end)=[];
end
function pixelIdxListIn=getInputMask(bioTree,nodeInfo)
for iIn=1:size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In,2)
    if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.isNode==true
        nodeInfo_pre=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.nodeInfo;
        pixelIdxListIn{iIn}=bioTree{nodeInfo_pre(1)}.node{nodeInfo_pre(2)}.Out{nodeInfo_pre(3)}.traceInfo.pixelIdxList{end};
    end
    if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.isNode==false
        rootInfo=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.rootInfo;
        if rootInfo(1)==nodeInfo(1)
            pixelIdxListIn{iIn}=bioTree{rootInfo(1)}.root{rootInfo(2)}.rootPixelDetail;
        end
        if  rootInfo(1)<nodeInfo(1)
            pixelIdxListIn{iIn}=bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.pixelIdxList{end};
        end
    end
end
end
function [newTrace,pixelIdxListIn,isBreak,error]=getNewTrace(newTrace,pixelIdxListNew,iframe,leafNum)
emptyCount=0;
nleaf=0;
for iIn=size(newTrace,2)+1:size(pixelIdxListNew,2)
    newTrace{iIn}.traceInfo.pixelIdxList=[];
    newTrace{iIn}.is2Node=[];
end
for iIn=1:size(pixelIdxListNew,2)
    if ~isempty(pixelIdxListNew{iIn});
        if ~isempty(newTrace{iIn}.is2Node)
            error=1;
            isBreak=[];
            pixelIdxListIn=[];
            return
        end
        newTrace{iIn}.traceInfo.pixelIdxList=[newTrace{iIn}.traceInfo.pixelIdxList,pixelIdxListNew(iIn)];
    else
        if isempty(newTrace{iIn}.is2Node)
            nleaf=nleaf+1;
            newTrace{iIn}.is2Node=false;
            newTrace{iIn}.leafInfo=[iframe,leafNum+nleaf];
            emptyCount= emptyCount+1;
        end
    end
end
if size(newTrace,2)-emptyCount==1
    isBreak=true;
else
    isBreak=false;
end
error=0;
pixelIdxListIn=pixelIdxListNew;
end

%% basicDivideSolution_deep
function [afterDivide,canDivideorNot,newRoot]=basicDivideSolution_deep(eachInfo,maskImage,pictureSize)
inputNum=numel(eachInfo);
realIn=0;
newRoot=0;
for i=1:inputNum
[~,xyMin,regionImage]=findRegionNum(eachInfo{i},pictureSize);
in{i}.xyMin=xyMin;
in{i}.BWImage=regionImage{1};
if ~isempty(in{i}.BWImage)
    realIn=realIn+1;
end
end
for i=1:size(in,2)
    image=imfill(in{i}.BWImage,'holes');
    if numel(image(image==1))>350
        afterDivide=[];
        canDivideorNot=0;
        return
    end
end
[outputNum,xyMin,regionImage]=findRegionNum(maskImage,pictureSize);
for i=1:outputNum
    out{i}.xyMin=xyMin;
    out{i}.BWImage=regionImage{i};
end
if realIn<outputNum
    if realIn+1==outputNum
        [inArea,~]=getSumArea(in,inputNum);
        [outArea,minArea]=getSumArea(out,outputNum);
        if (outArea-inArea)/minArea>0.6
            [afterDivide,canDivideorNot,newRoot]=n2mNodeSolution(in,out,pictureSize);
            return
        end
    end
    canDivideorNot=false;
    afterDivide=[];
    return
end
if realIn==outputNum
    [afterDivide,canDivideorNot,~]=n2mNodeSolution(in,out,pictureSize);
    return
end
if realIn>outputNum
    [inArea,minArea]=getSumArea(in,inputNum);
    [outArea,~]=getSumArea(out,outputNum);
    if (inArea-outArea)/(realIn-outputNum)/minArea>0.6
        [afterDivide,canDivideorNot]=n2mNodeSolution(in,out,pictureSize);
    else
        for i=1:size(out,2)
            out{1}.BWImage=out{1}.BWImage|out{i}.BWImage;
        end
        [afterDivide,canDivideorNot]=n21NodeDivide(in,out,pictureSize);
        memberNum=0;
        if canDivideorNot==1
            for i=1:size(afterDivide,2)
                if ~isempty(afterDivide{i})
                    memberNum=memberNum+1;
                end
            end
            if memberNum==0
                canDivideorNot=0;
            end
        end
    end
end
end
%% this  four function is used to divide N to 1 Node,by rotating the model and erode the out image
function [afterDivide,canDivideorNot]=n21NodeDivide(in,out,pictureSize)
sita=-30:2:30;
outInfo=out{1};
outImage=outInfo.BWImage;
outImage=imfill(outImage,'holes');
xyMin=outInfo.xyMin;
outEndImage=bwmorph(outImage,'thin',inf);
outEndImage=bwmorph(outEndImage,'endpoints');
[xInfo,yInfo]=find(outEndImage==1);
outEndInfo(:,1)=xInfo;
outEndInfo(:,2)=yInfo;
if size(outEndInfo,1)>=numel(in)
    [afterDivide,canDivideorNot]=ASolution(in,outInfo,pictureSize,outEndInfo);
else
    [afterDivide,canDivideorNot]=BSolution(in,sita,outInfo,outImage,pictureSize,xyMin);
end
end
%% this function is used to divide N to M Node directly,choose the nearest bacteria
function [afterDivide,isN2mNode,newRoot]=n2mNodeSolution(in,out,pictureSize)
% used to deal with N-to-M node and divide them,here N>=M
threShold=30;
inNum=numel(in);
for i=1:inNum
    if ~isempty(in{i}.BWImage)
    Info=regionprops(in{i}.BWImage,'Centroid','MajorAxisLength','MinorAxisLength');
    xyMin=in{i}.xyMin;
    inInfo(i,:)=[Info.Centroid(2)+xyMin(1)-1,Info.Centroid(1)+xyMin(2)-1,Info.MajorAxisLength,Info.MinorAxisLength];
    else
        inInfo(i,:)=zeros(1,4);
    end
end
outNum=numel(out);
for i=1:numel(out)
    Info=regionprops(out{i}.BWImage,'Centroid','MajorAxisLength','MinorAxisLength');
    xyMin=out{i}.xyMin;
    outInfo(i,:)=[Info.Centroid(2)+xyMin(1)-1,Info.Centroid(1)+xyMin(2)-1,Info.MajorAxisLength,Info.MinorAxisLength];
end
distanceMat=zeros(max(inNum,outNum));
distanceMat(1:inNum,1:outNum)=pdist2(inInfo,outInfo);
possibleSol=perms(1:max(inNum,outNum));
caculateDist=zeros(size(possibleSol,1),1);
for i=1:size(possibleSol,1)
    caculateDist(i)=trace(distanceMat(:,possibleSol(i,:)));
end
[minDist,Num]=min(caculateDist);
if minDist>threShold*outNum
    isN2mNode=false;
    afterDivide=[];
    newRoot=0;
else
    isN2mNode=true;
    correctVec=possibleSol(Num,:);
    if inNum>=outNum
        newRoot=0;
    else
        newRoot=1;
    end
    for i=1:max(inNum,outNum)
        if correctVec(i)>outNum
            afterDivide{i}=[];
        else
            afterDivide{i}=xy2Idx(out{correctVec(i)}.xyMin,out{correctVec(i)}.BWImage,pictureSize);
        end
    end
end
end
%% here are some basic functions
function [regionNum,xyMin,regionImage]=findRegionNum(pixelIdxList,pictureSize)
% find how many regions there are in a input/output
if isempty(pixelIdxList)
    regionNum=0;
    xyMin=zeros(1,2);
    regionImage{1}=[];
else
    [xyMin,BWImage]=idx2Xy(pixelIdxList,pictureSize);
    CC=bwconncomp(BWImage);
    regionNum=CC.NumObjects;
    for i=1:regionNum
        pixelIdxList2=CC.PixelIdxList{i};
        BW=false(CC.ImageSize);
        BW(pixelIdxList2)=1;
        regionImage{i}=BW;
    end
end
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
function result=AminusRow(A,row)
%used to make a differ between matrixA and a row vector
rowNum=size(A,1);
B=repmat(row,rowNum,1);
result=A-B;
end
function [sumArea,minArea]=getSumArea(info,num)
areaInfo=zeros(num,1);
for i=1:num
    filledImage=imfill(info{i}.BWImage,'holes');
    areaInfo(i)=numel(filledImage(filledImage==1));
end
sumArea=sum(areaInfo);
minArea=min(areaInfo(areaInfo~=0));
end
function newA=gainSameSize(A,B,pictureSize)
A.BWImage=imfill(A.BWImage,'holes');
% change picture in A to the same size that in B
pixelIdxList=xy2Idx(A.xyMin,A.BWImage,pictureSize);
new=false(pictureSize);
new(pixelIdxList)=true;
bSize=size(B.BWImage);
newA=new(B.xyMin(1)-1:B.xyMin(1)+bSize(1)-2,B.xyMin(2)-1:B.xyMin(2)+bSize(2)-2);
end
function I=getCentroid(pictureSize,pixelIdxList)
xSize=pictureSize(1);
yresult=ceil(pixelIdxList/xSize);
xresult=round(pixelIdxList-(yresult-1)*xSize);
yCentroid=round(mean(yresult));
xCentroid=round(mean(xresult));
centroidLinearIndex=xCentroid+(yCentroid-1)*xSize;
% centroidLinearIndex=round(mean(pixelIdxList));
I=false(pictureSize);
I(centroidLinearIndex)=1;
end