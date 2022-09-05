function [bioTree,branchList,allList]=divisionFinder(bioTree,branchList)
divisionThreshold=0.2;
singleCellSizeUpLimited=600;
frameThreshold=10;  %% change by jzy 5.19 pre(450)
hyperNodeNum=zeros(size(branchList,1),1);
divisionNum=zeros(size(branchList,1),1);
allListNode=[];
allListRoot=[];
allListLeaf=[];
for iList=1:size(branchList,1)
    branchInfo=branchList(iList,1:3);
    if branchInfo(3)==1
        allNode=bioTree{branchInfo(1)}.node{branchInfo(2)}.allNode;
        allRoot=bioTree{branchInfo(1)}.node{branchInfo(2)}.allRoot;
        allLeaf=bioTree{branchInfo(1)}.node{branchInfo(2)}.allLeaf;
%         [divisionCondition,hyperNodeList]=isDivision(bioTree,allNode,divisionThreshold,frameThreshold,singleCellSizeUpLimited);
        [divisionCondition,hyperNodeList]=isDivsionNode(bioTree,allNode);
        allNode(:,4)=divisionCondition;
        allNode(:,5)=hyperNodeList;
        bioTree{branchInfo(1)}.node{branchInfo(2)}.allNode=allNode;
    else
        if branchInfo(3)==0
            allNode=bioTree{branchInfo(1)}.root{branchInfo(2)}.allNode;
            allRoot=bioTree{branchInfo(1)}.root{branchInfo(2)}.allRoot;
            allLeaf=bioTree{branchInfo(1)}.root{branchInfo(2)}.allLeaf;
            divisionCondition=[];
            hyperNodeList=[];
        end
    end
hyperNodeNum(iList,1)=sum(hyperNodeList);
divisionNum(iList,1)=sum(divisionCondition);
branchIndexRoot=ones(size(allRoot,1),1).*iList;
branchIndexLeaf=ones(size(allLeaf,1),1).*iList;
branchIndexNode=ones(size(allNode,1),1).*iList;
allListNode=[allListNode;[allNode,branchIndexNode]];
allListRoot=[allListRoot;[allRoot,branchIndexRoot]];
allListLeaf=[allListLeaf;[allLeaf,branchIndexLeaf]];
end
branchList(:,5)=divisionNum;
branchList(:,6)=hyperNodeNum;
% [allListRoot,allListLeaf,allListNode]=getCentroid(bioTree,allListRoot,allListLeaf,allListNode);
allList.allRoot=allListRoot;
allList.allLeaf=allListLeaf;
allList.allNode=allListNode;
end
function [divisionCondition1,hyperNodeList]=isDivision(bioTree,allNode,divisionThreshold,frameThreshold,singleCellSizeUpLimited)
hyperNodeList=false(size(allNode,1),1);
divisionCondition1=false(size(allNode,1),1);
divisionCondition2=true(size(allNode,1),1);
for iList=1:size(allNode,1)
    nodeInfo=allNode(iList,:);
    hyperNode=isHyperNode(bioTree,nodeInfo);
    if hyperNode==true
        hyperNodeList(iList,1)=true;
        bioTree{nodeInfo(1)}.node{nodeInfo(2)}.isHyperNode=true;
        continue;
    end
    if divisionCondition2(iList,1)==false
        continue;
    end
    if ~(size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In,2)==1 && size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out,2)<=2)
        continue;
    end
    if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{1}.isNode==true
        rootCount=0;
        parentNodeInfo=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{1}.nodeInfo;
        parentHyperNode=isHyperNode(bioTree,parentNodeInfo);
        if parentHyperNode==true
            continue;
        end
        for iIn=1:size(bioTree{parentNodeInfo(1)}.node{parentNodeInfo(2)}.In,2)
            if bioTree{parentNodeInfo(1)}.node{parentNodeInfo(2)}.In{iIn}.isNode==false
                rootCount=rootCount+1;
            end
        end
        if rootCount>1
            continue;
        end
        if size(bioTree{parentNodeInfo(1)}.node{parentNodeInfo(2)}.In,2)==2 && (nodeInfo(1)-parentNodeInfo(1))<=frameThreshold
             continue;
        end
    end
    isTrue=filledAreaCheck(bioTree,nodeInfo,divisionThreshold,singleCellSizeUpLimited);
    if isTrue==true;
        divisionCondition1(iList,1)=true;
        divisionCondition2=findSubNodeinAllNode(bioTree,nodeInfo,allNode,divisionCondition2,frameThreshold);
        continue;
    end
end
end
function hyperNode=isHyperNode(bioTree,nodeInfo)
branchIndexList=[];
for iIn=1:size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In,2)
    branchIndex=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.branchIndex;
    branchIndexList=[branchIndexList,branchIndex];
end
if size(unique(branchIndexList),2)==1
   hyperNode=false;
else
   hyperNode=true;
end
end
function isTrue=filledAreaCheck(bioTree,nodeInfo,divisionThreshold,singleCellSizeUpLimited)
if size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out,2)==2
pixelList1=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1}.traceInfo.pixelIdxList{1};
pixelList2=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{2}.traceInfo.pixelIdxList{1};
if isequal(pixelList1,pixelList2)
    if size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1}.traceInfo.measurment{1},1)~=2
        isTrue=false;
        return;
    end
    aera1=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1}.traceInfo.measurment{1}(1).FilledArea;
    aera2=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1}.traceInfo.measurment{1}(2).FilledArea;
    if abs(aera1-aera2)/(aera1+aera2)<=divisionThreshold && (aera1+aera2)<=singleCellSizeUpLimited
        isTrue=true;
        return;
    else
       isTrue=false;
        return;
    end
end
if ~(size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1}.traceInfo.measurment{1},1)==1 && size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{2}.traceInfo.measurment{1},1)==1)
    isTrue=false;
    return;
end
aera1=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1}.traceInfo.measurment{1}(1).FilledArea;
aera2=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{2}.traceInfo.measurment{1}(1).FilledArea;
if abs(aera1-aera2)/(aera1+aera2)<=divisionThreshold && (aera1+aera2)<=singleCellSizeUpLimited
    isTrue=true;
    return;
else
    isTrue=false;
    return;
end
end
if size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out,2)==1
    if size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1}.traceInfo.measurment{1},1)~=2
        isTrue=false;
        return;
    end
    aera1=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1}.traceInfo.measurment{1}(1).FilledArea;
    aera2=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1}.traceInfo.measurment{1}(2).FilledArea;
    if abs(aera1-aera2)/(aera1+aera2)<=divisionThreshold && (aera1+aera2)<=singleCellSizeUpLimited
        isTrue=true;
        return;
    else
        isTrue=false;
        return;
    end
end
end
function  divisionCondition2=findSubNodeinAllNode(bioTree,nodeInfo,allNode,divisionCondition2,frameThreshold)
nodeList=nodeInfo;
allSubList=[];
while true
subNodeList=subNodeListSearching(bioTree,nodeList,frameThreshold);
if isempty(subNodeList)
    break;
end
subNodeList=listRuduction(subNodeList);
[~,indexList]=findNodeInList(allNode,subNodeList,divisionCondition2);
nodeList=allNode(indexList,1:3);
allSubList=[allSubList;nodeList];
end
allSubList=listRuduction(allSubList);
[divisionCondition2,~]=findNodeInList(allNode,allSubList,divisionCondition2);
end
function subNodeList=subNodeListSearching(bioTree,nodeList,frameThreshold)
subNodeList=[];
for iList=1:size(nodeList,1)
    nodeInfo=nodeList(iList,:);
    for iOut=1:size( bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out,2)
        if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.is2Node==true
            subNodeInfo=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.nodeInfo;
            if subNodeInfo(1)-nodeInfo(1)<=frameThreshold
                subNodeList=[subNodeList;subNodeInfo];
            end
        end
    end
end
end
function [divisionCondition2,tempList3]=findNodeInList(allNode,allSubList,divisionCondition2)
if isempty(allSubList)
    tempList3=[];
    return;
end
for iList=1:size(allSubList,1)
tempList1=allSubList(iList,1);
tempList2=allSubList(iList,2);
tempList3=(allNode(:,1)==tempList1 & allNode(:,2)==tempList2);
divisionCondition2=divisionCondition2 & (~tempList3);
end
end
function [allListRoot,allListLeaf,allListNode]=getCentroid(bioTree,allListRoot,allListLeaf,allListNode)
rootCentroid=zeros(size(allListRoot,1),2);
leafCentroid=zeros(size(allListLeaf,1),2);
nodeCentroid=zeros(size(allListNode,1),2);
for iroot=1:size(allListRoot,1)
    rootInfo=allListRoot(iroot,:);
    rootCentroid(iroot,1)=bioTree{rootInfo(1)}.root{rootInfo(2)}.rootMeasurment(1).Centroid(1);
    rootCentroid(iroot,2)=bioTree{rootInfo(1)}.root{rootInfo(2)}.rootMeasurment(1).Centroid(2);
end
for iLeaf=1:size(allListLeaf,1)
    leafInfo=allListLeaf(iLeaf,:);
    leafCentroid(iLeaf,1)=bioTree{leafInfo(1)}.leavies{leafInfo(2)}.leafMeasurment(1).Centroid(1);
    leafCentroid(iLeaf,2)=bioTree{leafInfo(1)}.leavies{leafInfo(2)}.leafMeasurment(1).Centroid(2);
end
for iNode=1:size(allListNode,1)
    nodeInfo=allListNode(iNode,:);
    nodeCentroid(iNode,1)=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1}.traceInfo.measurment{1}(1).Centroid(1);
    nodeCentroid(iNode,2)=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1}.traceInfo.measurment{1}(1).Centroid(2);
end
allListRoot(:,4:5)=rootCentroid;
allListLeaf(:,4:5)=leafCentroid;
allListNode(:,7:8)=nodeCentroid;
[~,sortRoot]=sort(allListRoot(:,1));
[~,sortLeaf]=sort(allListLeaf(:,1));
[~,sortNode]=sort(allListNode(:,1));
allListRoot=allListRoot(sortRoot,:);
allListLeaf=allListLeaf(sortLeaf,:);
allListNode=allListNode(sortNode,:);
end
function nodeListafter=listRuduction(nodeListbefore)
nodeListafter=[];
while ~isempty(nodeListbefore)
    nodeListafter=[nodeListafter;nodeListbefore(1,:)];
    tempListFrame=nodeListbefore(:,1)-nodeListbefore(1,1);
    tempListNode=nodeListbefore(:,2)-nodeListbefore(1,2);
    nodeListbefore=nodeListbefore(tempListFrame~=0|tempListNode~=0,:);
end
end
function [divisionCondition,hyperNodeList]=isDivsionNode(bioTree,allNode)
hyperNodeList=zeros(size(allNode,1),1);
divisionCondition=zeros(size(allNode,1),1);
for iList=1:size(allNode,1)
    nodeInfo=allNode(iList,:);
    if size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In,2)==1 && size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out,2)==2
%         frameOut=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out;
%         bacteriaSize=[frameOut{1}.traceInfo.measurment{1}.FilledArea,frameOut{2}.traceInfo.measurment{1}.FilledArea];
%         if ~(max(bacteriaSize)-min(bacteriaSize)>150 || max(bacteriaSize)>min(bacteriaSize)*2 )
            divisionCondition(iList)=1;
%             continue
%         end
    end
end
end