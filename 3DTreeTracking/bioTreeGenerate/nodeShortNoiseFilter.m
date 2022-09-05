function [bioTree,nodeNum]=nodeShortNoiseFilter(bioTree,frameShift,frameThreshold) %this function use to remove all Out in i Node to connect in j node, might be need further debug
nodeNumberBefore=nodeCounter(bioTree,frameShift);
nodeNum=nodeNumberBefore;
% disp(strcat('orignal node number=',num2str(nodeNumber)));
while true
    fprintf('.');
for iframe=1+frameShift:size(bioTree,2)
    if ~isempty(bioTree{iframe}.node)
        countNode=1;
        while true
            [nodetest,nextNodeInfo]=isTrueNode(bioTree,iframe,countNode,frameThreshold);
            if nodetest==false
                bioTree=removeFakeNode(bioTree,iframe,countNode,nextNodeInfo);
                bioTree=nodeCorrect(bioTree,nextNodeInfo);
            end
            if countNode>=size(bioTree{iframe}.node,2)
                break;
            end
            countNode=countNode+1;
        end
    end
end
nodeNumberAfter=nodeCounter(bioTree,frameShift);
if nodeNumberBefore==nodeNumberAfter
    nodeNum=[nodeNum,nodeNumberAfter];
    break;
else
    nodeNumberBefore=nodeNumberAfter;
end
% disp(strcat('afterRefine node number=',num2str(nodeNumber)));
end
end
function [nodetest,nextNodeInfo]=isTrueNode(bioTree,iframe,iNode,frameThreshold)
nextNodeInfo=[];
parentNodeInfo=[];
for iOut=1:size(bioTree{iframe}.node{iNode}.Out,2)
    if bioTree{iframe}.node{iNode}.Out{iOut}.is2Node==false
        nodetest=true;
        return;
    end
    if bioTree{iframe}.node{iNode}.Out{iOut}.is2Node==true
        sunNodeInfo=bioTree{iframe}.node{iNode}.Out{iOut}.nodeInfo;
        if (sunNodeInfo(1)-iframe)>frameThreshold
            nodetest=true;
            return;
        else
            nextNodeInfo=[nextNodeInfo;sunNodeInfo];
        end
    end
end
nextNodeInfo=listRuduction(nextNodeInfo);
if size(nextNodeInfo,1)>1
    nodetest=true;
    return;
end
if size(nextNodeInfo,1)==1
    for iIn=1:size(bioTree{nextNodeInfo(1)}.node{nextNodeInfo(2)}.In,2)
        if bioTree{nextNodeInfo(1)}.node{nextNodeInfo(2)}.In{iIn}.isNode==false
            nodetest=true;
            return;
        end
        if bioTree{nextNodeInfo(1)}.node{nextNodeInfo(2)}.In{iIn}.isNode==true
            sunNodeInfo=bioTree{nextNodeInfo(1)}.node{nextNodeInfo(2)}.In{iIn}.nodeInfo;
            parentNodeInfo=[parentNodeInfo;sunNodeInfo];
        end
    end
    parentNodeInfo=listRuduction(parentNodeInfo);
    if size(parentNodeInfo,1)==1
        nodetest=false;
        return;
    else
        nodetest=true;
        return;
    end
end
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
function bioTree=removeFakeNode(bioTree,iframe,iNode,nextNodeInfo)
sizeStack1=size(bioTree{iframe}.node{iNode}.Out{1}.traceInfo.pixelIdxList,2);
for itrace=1:sizeStack1
    sumPixelList=[];
    sumMeasurment=[];
    for iNodeOut=1:size(bioTree{iframe}.node{iNode}.Out,2)
        sumPixelList=[sumPixelList;bioTree{iframe}.node{iNode}.Out{iNodeOut}.traceInfo.pixelIdxList{itrace}];
        sumMeasurment=[sumMeasurment;bioTree{iframe}.node{iNode}.Out{iNodeOut}.traceInfo.measurment{itrace}];
    end
    outPixelIdxList{itrace}=sumPixelList;
    outMeasurment{itrace}=sumMeasurment;
end
bioTree{iframe}.node{iNode}.Out=[];
for iOut=1:size(bioTree{nextNodeInfo(1)}.node{nextNodeInfo(2)}.Out,2)
    sizeStack2=size(bioTree{nextNodeInfo(1)}.node{nextNodeInfo(2)}.Out{iOut}.traceInfo.pixelIdxList,2);
    bioTree{iframe}.node{iNode}.Out{iOut}=bioTree{nextNodeInfo(1)}.node{nextNodeInfo(2)}.Out{iOut};
    bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.pixelIdxList(1:sizeStack1)=outPixelIdxList;
    bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.pixelIdxList(sizeStack1+1:sizeStack1+sizeStack2)=bioTree{nextNodeInfo(1)}.node{nextNodeInfo(2)}.Out{iOut}.traceInfo.pixelIdxList;
    bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.measurment(1:sizeStack1)=outMeasurment;
    bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.measurment(sizeStack1+1:sizeStack1+sizeStack2)=bioTree{nextNodeInfo(1)}.node{nextNodeInfo(2)}.Out{iOut}.traceInfo.measurment;
    if bioTree{iframe}.node{iNode}.Out{iOut}.is2Node==false
        leafInfo=bioTree{iframe}.node{iNode}.Out{iOut}.leafInfo;
        bioTree{leafInfo(1)}.leavies{leafInfo(2)}.nodeInfo=[iframe,iNode,iOut];
    end
    if bioTree{iframe}.node{iNode}.Out{iOut}.is2Node==true
        subNodeInfo=bioTree{iframe}.node{iNode}.Out{iOut}.nodeInfo;
        bioTree{subNodeInfo(1)}.node{subNodeInfo(2)}.In{subNodeInfo(3)}.nodeInfo=[iframe,iNode,iOut];
    end
end
end
function bioTree=nodeCorrect(bioTree,nextNodeInfo)
bioTree{nextNodeInfo(1)}.node(nextNodeInfo(2))=[];

if ~isempty(bioTree{nextNodeInfo(1)}.node)
    for iNode=1:size(bioTree{nextNodeInfo(1)}.node,2)
        for iNodeIn=1:size(bioTree{nextNodeInfo(1)}.node{iNode}.In,2)
            if bioTree{nextNodeInfo(1)}.node{iNode}.In{iNodeIn}.isNode==false
                rootInfo=bioTree{nextNodeInfo(1)}.node{iNode}.In{iNodeIn}.rootInfo;
                bioTree{rootInfo(1)}.root{rootInfo(2)}.nodeInfo=[nextNodeInfo(1),iNode,iNodeIn];
            end
            if bioTree{nextNodeInfo(1)}.node{iNode}.In{iNodeIn}.isNode==true
                nodeInfo=bioTree{nextNodeInfo(1)}.node{iNode}.In{iNodeIn}.nodeInfo;
                bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.nodeInfo=[nextNodeInfo(1),iNode,iNodeIn];
            end
        end
        for iNodeOut=1:size(bioTree{nextNodeInfo(1)}.node{iNode}.Out,2)
            if bioTree{nextNodeInfo(1)}.node{iNode}.Out{iNodeOut}.is2Node==false
                leafInfo=bioTree{nextNodeInfo(1)}.node{iNode}.Out{iNodeOut}.leafInfo;
                bioTree{leafInfo(1)}.leavies{leafInfo(2)}.nodeInfo=[nextNodeInfo(1),iNode,iNodeOut];
            end
            if bioTree{nextNodeInfo(1)}.node{iNode}.Out{iNodeOut}.is2Node==true
                nodeInfo=bioTree{nextNodeInfo(1)}.node{iNode}.Out{iNodeOut}.nodeInfo;
                bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{nodeInfo(3)}.nodeInfo=[nextNodeInfo(1),iNode,iNodeOut];
            end
        end
    end
end
end
function nodeNumber=nodeCounter(bioTree,frameShift)
nodeNumber=0;
for iframe=1+frameShift:size(bioTree,2)
    if ~isempty(bioTree{iframe}.node)
        nodeNumber=nodeNumber+size(bioTree{iframe}.node,2);
    end
end
end