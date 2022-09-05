function [correctionRateEach,divAll]=bioTreeCorrectRate(bioTree)
minDivisionTime=800;  %% original 800
[bioTree,branchList,~,~]=myBiograph_new2(bioTree);
coreBranch=branchList(:,3)==1;
branchList=branchList(coreBranch,:);
correctionRateEach=[];
for ibranch=1:size(branchList,1)
    iList=branchList(ibranch,:);
    allNode=bioTree{iList(1)}.node{iList(2)}.allNode;
    nodeNum=size(allNode,1);
    for iNode=1:size(allNode,1)
        allNode(iNode,4)=morphCorrection(bioTree,allNode(iNode,:));
        bioTree{allNode(iNode,1)}.node{allNode(iNode,2)}.divisionNode=allNode(iNode,4);
    end
    allNode=allNode(allNode(:,4)==1,:);
    for iNode=1:size(allNode,1)
        nodeInfo=allNode(iNode,1:2);
        if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{1}.isNode==0
            continue
        end
        preNode=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{1}.nodeInfo;
        if bioTree{preNode(1)}.node{preNode(2)}.divisionNode==0
            continue
        end
        if nodeInfo(1)-preNode(1)<=minDivisionTime
            allNode(iNode,4)=0;
            bioTree{nodeInfo(1)}.node{nodeInfo(2)}.divisionNode=0;
        end
    end
    isDiv=allNode(:,4);
    divisionNum=numel(isDiv(isDiv==1));
    correctionRateEach=[correctionRateEach;nodeNum,divisionNum,divisionNum/nodeNum];
end
nodeNum=0;
divNode=0;
for iframe=1:size(bioTree,2)
    for iNode=1:size(bioTree{iframe}.node,2)
        nodeNum=nodeNum+1;
        if bioTree{iframe}.node{iNode}.divisionNode==1
            divNode=divNode+1;
        end
    end
end
divAll=[nodeNum,divNode,divNode/nodeNum];
end
function rightNode=morphCorrection(bioTree,nodeInfo)
    pictureSize=bioTree{1}.imageSize;
    if ~(size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In,2)==1 && size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out,2)==2)
        rightNode=0;
        return
    end
    rightNode=1;
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
BWImageGain=imfill(BWImageGain,'holes');
end