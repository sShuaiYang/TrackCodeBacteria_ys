function [linkMatrix,centroidInfo,leafList,leafNum,allList,getLeaf,allLeaf,linkTwoPoint]=clusterLikeOneBranchPhytree(bioTree,iBra,ibranch,branchList)
% give mark to each bioTree node and find the right node
minDivisionTime=600;  %% original 800
% [bioTree,branchList,~,~]=myBiograph_new2(bioTree);
coreBranch=branchList(:,3)==1;
branchList=branchList(coreBranch,:);
for ibac=1:size(branchList,1)
    iList=branchList(ibac,:);
    allNode=bioTree{iList(1)}.node{iList(2)}.allNode;
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
    bioTree{iList(1)}.node{iList(2)}.allNode=allNode; %将是不是divisionNode的信息放在allNode的第四列
end

% traceMatrix means the current set of divide or leaf
% (:,3)=0 means not need to trace, =1 means the opposite
% (:,4)=-1 means root, 0 leaf ,1 node, 2 final leaf

% 三种节点  node和root都称为node, 直连的node为leaf, 非直连的为otherLeaf
if ibranch(3)~=0
    allRoot=bioTree{ibranch(1)}.node{ibranch(2)}.allRoot(1,:);
    allLeaf=bioTree{ibranch(1)}.node{ibranch(2)}.allLeaf;
    allNode=bioTree{ibranch(1)}.node{ibranch(2)}.allNode;
    isDiv=allNode(:,4);
    divisionNum=numel(isDiv(isDiv==1));
    if divisionNum==size(allNode,1)
        getLeaf=1;
    else
        getLeaf=0;
    end
else
    allRoot=bioTree{ibranch(1)}.root{ibranch(2)}.allRoot(1,:);
    allLeaf=bioTree{ibranch(1)}.root{ibranch(2)}.allLeaf;
    getLeaf=1;
end
allLeaf(:,3)=0;
allLeaf(:,4)=0;
allRoot(:,3)=1;
allRoot(:,4)=-1;
traceMatrix{1}=allRoot;
i=1;
nodeList=[];
while  ~all(traceMatrix{i}(:,3)==0)
    i=i+1;
    traceMatrix{i}=[];
    lastList=traceMatrix{i-1};
    for iNum=1:size(lastList,1)
        if i==2
            nextIsNode=bioTree{lastList(iNum,1)}.root{lastList(iNum,2)}.is2Node;
            if nextIsNode==0 
                nextLeaf=bioTree{lastList(iNum,1)}.root{lastList(iNum,2)}.leafInfo;
                traceMatrix{i}=[traceMatrix{i};nextLeaf(1),nextLeaf(2),0,0];
            end
            if nextIsNode==1
                nextNode=bioTree{lastList(iNum,1)}.root{lastList(iNum,2)}.nodeInfo;
                traceMatrix{i}=[traceMatrix{i};nextNode(1),nextNode(2),1,1];
            end
        end
        if i>=3
            if lastList(iNum,3)==0
                traceMatrix{i}=[traceMatrix{i};lastList(iNum,:)];
            end
            if lastList(iNum,3)==1
                outNum=size(bioTree{lastList(iNum,1)}.node{lastList(iNum,2)}.Out,2);
                if ~(size(bioTree{lastList(iNum,1)}.node{lastList(iNum,2)}.Out,2)==2 && size(bioTree{lastList(iNum,1)}.node{lastList(iNum,2)}.In,2)==1)
                    traceMatrix{i}=[traceMatrix{i};lastList(iNum,1),lastList(iNum,2),0,lastList(iNum,4)];
                    continue
                end
                isDivNode=[];
                isFinalLeaf=[];
                for iOut=1:2
                    nextIsNode=bioTree{lastList(iNum,1)}.node{lastList(iNum,2)}.Out{iOut}.is2Node;
                    if nextIsNode==0
                        nextLeaf=bioTree{lastList(iNum,1)}.node{lastList(iNum,2)}.Out{iOut}.leafInfo;
                        isDivNode=[isDivNode;1];
                        isFinalLeaf=[isFinalLeaf;nextLeaf(1)==numel(bioTree)];
                    else
                        isFinalLeaf=[isFinalLeaf;0];
                        nextNode=bioTree{lastList(iNum,1)}.node{lastList(iNum,2)}.Out{iOut}.nodeInfo;
                        if bioTree{nextNode(1)}.node{nextNode(2)}.divisionNode==1 
                            isDivNode=[isDivNode;1];
                        else
                            isDivNode=[isDivNode;0];
                        end
                    end
                end
                if ~all(isDivNode==1) || (any(isFinalLeaf==1) && getLeaf==0)
                    traceMatrix{i}=[traceMatrix{i};lastList(iNum,1),lastList(iNum,2),0,lastList(iNum,4)];
                    continue
                end
                noChangeNum=0;
                for iOut=1:outNum
                    nextIsNode=bioTree{lastList(iNum,1)}.node{lastList(iNum,2)}.Out{iOut}.is2Node;
                    if nextIsNode==0
                        leafInfo=bioTree{lastList(iNum,1)}.node{lastList(iNum,2)}.Out{iOut}.leafInfo;
                        if ~isempty(traceMatrix{i}) && any(leafInfo(1)==traceMatrix{i}(:,1) & leafInfo(2)==traceMatrix{i}(:,2) & 0==traceMatrix{i}(:,4))
                            noChangeNum=noChangeNum+1;
                            if noChangeNum==1
                                traceMatrix{i}=[traceMatrix{i};lastList(iNum,1),lastList(iNum,2),0,lastList(iNum,4)];
                            end
                        else
                            traceMatrix{i}=[traceMatrix{i};leafInfo(1),leafInfo(2),0,0];
                        end
                    end
                    if nextIsNode==1
                        nodeInfo=bioTree{lastList(iNum,1)}.node{lastList(iNum,2)}.Out{iOut}.nodeInfo;
                        if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.branchIndex~=iBra
                            noChangeNum=noChangeNum+1;
                            if noChangeNum==1
                                traceMatrix{i}=[traceMatrix{i};lastList(iNum,1),lastList(iNum,2),0,lastList(iNum,4)];
                            end
                            continue
                        end
                        if ~isempty(traceMatrix{i}) && any(nodeInfo(1)==traceMatrix{i}(:,1) & nodeInfo(2)==traceMatrix{i}(:,2) & 1==traceMatrix{i}(:,4))
                            noChangeNum=noChangeNum+1;
                            if noChangeNum==1
                                traceMatrix{i}=[traceMatrix{i};lastList(iNum,1),lastList(iNum,2),0,lastList(iNum,4)];
                            end
                        else
                            if isempty(nodeList) || ~isempty(nodeList) && ~(any(nodeList(:,1)==nodeInfo(1) & nodeList(:,2)==nodeInfo(2)))
                                traceMatrix{i}=[traceMatrix{i};nodeInfo(1),nodeInfo(2),1,1];
                            end
                        end
                    end
                end
            end
        end
    end
    for iNum=1:size(traceMatrix{i},1)
        if traceMatrix{i}(iNum,3)==1 && traceMatrix{i}(iNum,4)==1
            nodeList=[nodeList;traceMatrix{i}(iNum,1),traceMatrix{i}(iNum,2),1,1];
        end
    end
end
nodeList=[allRoot;nodeList];
[~,sortOrder]=sort(nodeList(:,1));
nodeList=nodeList(sortOrder,:);

% if size(nodeList,1)==leafNum-1
nodeList=nodeList;
leafList=traceMatrix{end}(traceMatrix{end}(:,3)==0,:);
leafNum=size(leafList,1);

deleteNode=[];
for iNode=1:size(nodeList,1)
    nodeInfo=nodeList(iNode,:);
    if leafList(leafList(:,1)==nodeInfo(1) & leafList(:,2)==nodeInfo(2),4)==1
        deleteNode=[deleteNode;iNode];
    end
end
nodeList(deleteNode,:)=[];
nodeNum=size(nodeList,1);
% 判断一下，如果是100%正确的tree就画出leaf，否则一个leaf都不连，全部都是最后的leaf
if getLeaf==0;
    allLeaf=allLeaf(allLeaf(:,1)==numel(bioTree),:);
else 
    allLeaf=leafList;
end

leafOrder=1:leafNum;
nodeOrder=1:nodeNum;
linkMatrix=zeros(leafNum+nodeNum);
for iNode=1:nodeNum
    if nodeList(iNode,4)==1
        nodeInfo=nodeList(iNode,:);
        outNum=size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out,2);
        for iOut=1:outNum
            nextIsNode=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.is2Node;
            if nextIsNode==0
                leafInfo=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.leafInfo;
                if ~(getLeaf==0 && leafInfo(1)==numel(bioTree))
                    orderNum=leafOrder(leafInfo(1)==leafList(:,1) & leafInfo(2)==leafList(:,2));
                    linkMatrix(iNode,nodeNum+orderNum)=1;
                end
            end
            if nextIsNode==1
                nextNode=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.nodeInfo;
                orderNum=nodeOrder(nextNode(1)==nodeList(:,1) & nextNode(2)==nodeList(:,2));
                linkMatrix(min(orderNum,iNode),max(orderNum,iNode))=1;
                if isempty(orderNum)
                    orderNum=leafOrder(nextNode(1)==leafList(:,1) & nextNode(2)==leafList(:,2));
                    if leafList(orderNum,4)==1                        
                        linkMatrix(iNode,nodeNum+orderNum)=1;
                    end
                end
            end
        end
    end
    if nodeList(iNode,4)==-1
        rootInfo=nodeList(iNode,:);
        nextIsNode=bioTree{rootInfo(1)}.root{rootInfo(2)}.is2Node;
        if nextIsNode==0
            leafInfo=bioTree{rootInfo(1)}.root{rootInfo(2)}.leafInfo;
            orderNum=leafOrder(leafInfo(1)==leafList(:,1) & leafInfo(2)==leafList(:,2));
            linkMatrix(iNode,nodeNum+orderNum)=1;
        end
        if nextIsNode==1
            nextNode=bioTree{rootInfo(1)}.root{rootInfo(2)}.nodeInfo;
            orderNum=nodeOrder(nextNode(1)==nodeList(:,1) & nextNode(2)==nodeList(:,2));
            linkMatrix(min(orderNum,iNode),max(orderNum,iNode))=1;
        end
    end
end
if ~isempty(nodeList)
    centroidInfo=[nodeList(:,1);leafList(:,1)];
else
    centroidInfo=leafList(:,1);
end
if getLeaf==0
    lastNode=find(leafList(:,4)==1);
    linkTwoPoint=[lastNode(1)+nodeNum,lastNode(end)+nodeNum];
else
    linkTwoPoint=[];
end
    allList=[nodeList;leafList];
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