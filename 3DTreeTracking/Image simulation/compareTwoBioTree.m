function score=compareTwoBioTree(bioTree,bioTree1)
[score(1,1),score(1,2)]=nodeComparation(bioTree,bioTree1);
[score(2,1),score(2,2)]=IrrelevantCellComparation(bioTree,bioTree1);
score(3,1)=sum(score(:,1));
score(3,2)=sum(score(:,2));
[score(4,1),score(4,2)]=AttachingCaseComparation(bioTree,bioTree1);
[score(5,1),score(5,2)]=DetachingCaseComparation(bioTree,bioTree1);
score(:,3)=score(:,2)./score(:,1);
end

%% AttachingCase Comparation
function [allCase,rightCase]=AttachingCaseComparation(bioTree,bioTree1)
attachingInfo=findAttachingCase(bioTree);
attachingInfo1=findAttachingCase(bioTree1);
rightCase=0;
allCase=0;
for iRoot=1:size(attachingInfo1)
    rootInfo1=attachingInfo1(iRoot,:);
    possibleNum=find(attachingInfo(:,1)==rootInfo1(1));
    if ~isempty(possibleNum)
        matchNum=findMatchRoot(bioTree,attachingInfo(possibleNum,:),bioTree1,rootInfo1);
        if ~isempty(matchNum)
            rightCase=rightCase+1;
        end
    end
    allCase=allCase+1;
end
end
function attachingInfo=findAttachingCase(bioTree)
num=0;
for iframe=2:size(bioTree,2)
    if ~isempty(bioTree{iframe}.root)
        for iRoot=1:size(bioTree{iframe}.root,2)
            num=num+1;
            attachingInfo(num,:)=[iframe,iRoot];
        end
    end
end
end

%% DetachingCase Comparation
function [allCase,rightCase]=DetachingCaseComparation(bioTree,bioTree1)
detachingInfo=findDetachingCase(bioTree);
detachingInfo1=findDetachingCase(bioTree1);
rightCase=0;
allCase=0;
for ileaf=1:size(detachingInfo1)
    leafInfo1=detachingInfo1(ileaf,:);
    possibleNum=find(detachingInfo(:,1)==leafInfo1(1));
    if ~isempty(possibleNum)
        matchNum=findMatchleaf(bioTree,detachingInfo(possibleNum,:),bioTree1,leafInfo1);
        if ~isempty(matchNum)
            rightCase=rightCase+1;
        end
    end
    allCase=allCase+1;
end
end    
function detachingInfo=findDetachingCase(bioTree)
num=0;
for iframe=1:size(bioTree,2)-1
    if ~isempty(bioTree{iframe}.leavies)
        for iLeaf=1:size(bioTree{iframe}.leavies,2)
            num=num+1;
            detachingInfo(num,:)=[iframe,iLeaf];
        end
    end
end
end
function matchNum=findMatchleaf(bioTree,possibleLeaf,bioTree1,leafInfo1)
imageSize=bioTree{1}.imageSize;
for i=1:size(possibleLeaf,1)
    leafInfo=possibleLeaf(i,:);
    leafPixel=bioTree{leafInfo(1)}.leavies{leafInfo(2)}.leaviesPixelDetail;
    image=false(imageSize);
    image(leafPixel)=true;
    cc=regionprops(image,'Centroid');
    possibleInfo(i,:)=cc.Centroid;
end
leafPixel=bioTree1{leafInfo1(1)}.leavies{leafInfo1(2)}.leaviesPixelDetail;
image=false(imageSize);
image(leafPixel)=true;
cc=regionprops(image,'Centroid');
targetInfo(1,:)=cc.Centroid;
radius=pdist2(possibleInfo,targetInfo);
if min(radius)<10
    matchNum=find(radius==min(radius));
else
    matchNum=[];
end
end

%% IrrelevantCell Comparation
function [allFrame,rightFrame]=IrrelevantCellComparation(bioTree,bioTree1)
infoList=findIrrelevantCell(bioTree);
infoList1=findIrrelevantCell(bioTree1);
rightFrame=0;
allFrame=0;
for iRoot=1:size(infoList1)
    rootInfo1=infoList1(iRoot,:);
    possibleNum=find(infoList(:,1)==rootInfo1(1));
    if ~isempty(possibleNum)
        matchNum=findMatchRoot(bioTree,infoList(possibleNum,:),bioTree1,rootInfo1);
        if ~isempty(matchNum)
            rootInfo=infoList(possibleNum(matchNum),:);
            if bioTree{rootInfo(1)}.root{rootInfo(2)}.leafInfo(1)==bioTree1{rootInfo1(1)}.root{rootInfo1(2)}.leafInfo(1)
                rightFrame=rightFrame+size(bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.pixelIdxList,2);
            end
        end
    end
    allFrame=allFrame+size(bioTree1{rootInfo1(1)}.root{rootInfo1(2)}.traceInfo.pixelIdxList,2);
end
end
function matchNum=findMatchRoot(bioTree,possibleRoot,bioTree1,rootInfo1)
imageSize=bioTree{1}.imageSize;
for i=1:size(possibleRoot,1)
    rootInfo=possibleRoot(i,:);
    rootPixel=bioTree{rootInfo(1)}.root{rootInfo(2)}.rootPixelDetail;
    image=false(imageSize);
    image(rootPixel)=true;
    cc=regionprops(image,'Centroid');
    possibleInfo(i,:)=cc.Centroid;
end
rootPixel=bioTree1{rootInfo1(1)}.root{rootInfo1(2)}.rootPixelDetail;
image=false(imageSize);
image(rootPixel)=true;
cc=regionprops(image,'Centroid');
targetInfo(1,:)=cc.Centroid;
radius=pdist2(possibleInfo,targetInfo);
if min(radius)<10
    matchNum=find(radius==min(radius));
else
    matchNum=[];
end
end
function infoList=findIrrelevantCell(bioTree)
infoList=[];
num=0;
for iFrame=1:size(bioTree,2)
   if ~isempty(bioTree{iFrame}.root)
       for iRoot=1:size(bioTree{iFrame}.root,2)
           num=num+1;
           if bioTree{iFrame}.root{iRoot}.is2Node==0
               infoList=[infoList;iFrame,iRoot];
           end
       end
   end
end
end
%% branch Comparation
function [allFrame,rightFrame]=nodeComparation(bioTree,bioTree1)
[bioTree,branchList,~,~]=myBiograph_new2(bioTree);
branchList(branchList(:,3)==0,:)=[];
[bioTree1,branchList1,~,~]=myBiograph_new2(bioTree1);
branchList1(branchList1(:,3)==0,:)=[];
rightFrame=0;
allFrame=0;
for ibranch=1:size(branchList1,1)
    rightNode=[];
    branchNode1=branchList1(ibranch,:);
    possibleNum=find(branchList(:,1)==branchNode1(1));
    nodeList1=changeList(bioTree1,bioTree1{branchNode1(1)}.node{branchNode1(2)}.allNode);
    if ~isempty(possibleNum)
%          nodeList=changeList(bioTree{branchNode(1)}.node{branchNode(2)}.allNode);
%          frameNum=branchsFrame(bioTree1,nodeList);
        matchNum=findMatchNode(bioTree,branchList(possibleNum,:),bioTree1,branchNode1);
        if ~isempty(matchNum)
            branchNode=branchList(possibleNum(matchNum),:);
            nodeList=bioTree{branchNode(1)}.node{branchNode(2)}.allNode;
            for i=1:size(nodeList1);
                eachNode1=nodeList1(i,:);
                possibleNodeNum=find(nodeList(:,1)==eachNode1(1));
                if ~isempty(possibleNodeNum)
                    possibleNode=changeList(bioTree,nodeList(possibleNodeNum,:));
                    matchNum=findEachMatchNode(bioTree,possibleNode,bioTree1,eachNode1);
                    if ~isempty(matchNum)
                        eachNode=possibleNode(matchNum,:);
                        nextFrame=findNextFrame(bioTree,eachNode);
                        nextFrame1=findNextFrame(bioTree1,eachNode1);
                        if nextFrame==nextFrame1;
                            rightNode=[rightNode;eachNode];
                        end
                    end
                end
            end
        end
    end
%     disp(ibranch)
    frameNum1=branchsFrame(bioTree1,nodeList1);
    allFrame=allFrame+frameNum1;
    if ~isempty(rightNode)
        frameNum=branchsFrame(bioTree,rightNode);
        rightFrame=rightFrame+frameNum;
    end
end
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
function frameNum=branchsFrame(bioTree,nodeList)
frameNum=0;
if bioTree{nodeList(1,1)}.node{nodeList(1,2)}.In{1}.isNode==0
rootInfo=bioTree{nodeList(1,1)}.node{nodeList(1,2)}.In{1}.rootInfo;
frameNum=frameNum+size(bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.pixelIdxList,2);
end
for iNode=1:size(nodeList,1)
    nodeInfo=nodeList(iNode,:);
    nodeSize=size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.pixelIdxList,2);
    frameNum=frameNum+nodeSize;
end
end
function nodeList=changeList(bioTree,allNode)
nodeList=[];
for i=1:size(allNode,1)
    nodeInfo=allNode(i,:);
    for t=1:size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out,2);
        nodeOutInfo=[nodeInfo(1),nodeInfo(2),t];
        nodeList=[nodeList;nodeOutInfo];
    end
end
end
function matchNum=findMatchNode(bioTree,possibleNode,bioTree1,branchNode)
imageSize=bioTree{1}.imageSize;
for i=1:size(possibleNode,1)
    pixelIdxListIn=getInputMask(bioTree,possibleNode(i,:));
    pixelIdxListIn=pixelIdxListIn{end};
    image=false(imageSize);
    image(pixelIdxListIn)=true;
    cc=regionprops(image,'Centroid');
    possibleInfo(i,:)=cc.Centroid;
end
pixelIdxListIn=getInputMask(bioTree1,branchNode);
pixelIdxListIn=pixelIdxListIn{end};
image=false(imageSize);
image(pixelIdxListIn)=true;
cc=regionprops(image,'Centroid');
targetInfo(1,:)=cc.Centroid;
radius=pdist2(possibleInfo,targetInfo);
if min(radius)<10
    matchNum=find(radius==min(radius));
else
    matchNum=[];
end
end
function matchNum=findEachMatchNode(bioTree,possibleNode,bioTree1,eachNode1)
imageSize=bioTree{1}.imageSize;
for i=1:size(possibleNode,1)
    nodeInfo=possibleNode(i,:);
    pixelIdxListOut=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.pixelIdxList{1};
    image=false(imageSize);
    image(pixelIdxListOut)=true;
    cc=regionprops(image,'Centroid');
    possibleInfo(i,:)=cc.Centroid;
end
pixelIdxListOut=bioTree1{eachNode1(1)}.node{eachNode1(2)}.Out{eachNode1(3)}.traceInfo.pixelIdxList{1};
image=false(imageSize);
image(pixelIdxListOut)=true;
cc=regionprops(image,'Centroid');
targetInfo(1,:)=cc.Centroid;
radius=pdist2(possibleInfo,targetInfo);
if min(radius)<10
    matchNum=find(radius==min(radius));
else
    matchNum=[];
end
end
function nextFrame=findNextFrame(bioTree,nodeInfo)
if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.is2Node==1
    nextNode=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.nodeInfo;
    nextFrame=nextNode(1);
else
    nextLeaf=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.leafInfo;
    nextFrame=nextLeaf(1);
end
end