function [bioTree,nodeNumber]=nodeRefineBioTree(bioTree,frameShift,frameThreshold) %this function use to remove all Out in i Node to connect in j node, might be need further debug
nodeNumber0=nodeCounter(bioTree,frameShift);
for iframe=1+frameShift:size(bioTree,2)
    countNode=0;
    nodeClear=[];
    if ~isempty(bioTree{iframe}.node)
        for iNode=1:size(bioTree{iframe}.node,2)
            [nodetest1,nextnodeInfo1]=isTrueNodeOut(bioTree{iframe}.node{iNode},iframe,frameThreshold);
            [nodetest2,nextnodeInfo2]=isTrueNodeIn(bioTree{iframe}.node{iNode},iframe,frameThreshold);
            if nodetest1==1 && nodetest2==0
                bioTree=removeFakeNode1(bioTree{iframe}.node{iNode},nextnodeInfo1,bioTree);
                nodeClear=[nodeClear,iNode];
                countNode=countNode+1;
            end
            if nodetest2==1 && nodetest1==0
                bioTree=removeFakeNode2(bioTree{iframe}.node{iNode},nextnodeInfo2,bioTree);
                nodeClear=[nodeClear,iNode];
                countNode=countNode+1;
            end
        end
        if countNode>0
            bioTree{iframe}.node(nodeClear)=[];
            bioTree=nodeCorrect(bioTree,iframe);
        end
    end
end
nodeNumber1=nodeCounter(bioTree,frameShift);
nodeNumber=[nodeNumber0,nodeNumber1];
end
function [nodetest,nextnodeInfo]=isTrueNodeOut(node,iframe,frameThreshold)
nextNodeList=[];
for iOut=1:size(node.Out,2)
    if node.Out{iOut}.is2Node==false
        nodetest=0;
        nextnodeInfo=[];
        return;
    else
        if (node.Out{iOut}.nodeInfo(1)-iframe)>frameThreshold
            nodetest=0;
            nextnodeInfo=[];
            return;
        else
            nextNodeList=[nextNodeList;node.Out{iOut}.nodeInfo];
        end
    end
end
if size(nextNodeList,1)==1
    nodetest=1;
    nextnodeInfo=nextNodeList;
    return;
else
    nodetest=0;
    nextnodeInfo=[];
    return;
end
end
function [nodetest,nextnodeInfo]=isTrueNodeIn(node,iframe,frameThreshold)
nextNodeList=[];
for iIn=1:size(node.In,2)
    if node.In{iIn}.isNode==false
        nodetest=0;
        nextnodeInfo=[];
        return;
    else
        if (iframe-node.In{iIn}.nodeInfo(1))>frameThreshold
            nodetest=0;
            nextnodeInfo=[];
            return;
        else
            nextNodeList=[nextNodeList;node.In{iIn}.nodeInfo];
        end
    end
end
if size(nextNodeList,1)==1
    nodetest=1;
    nextnodeInfo=nextNodeList;
    return;
else
    nodetest=0;
    nextnodeInfo=[];
    return;
end
end
function bioTree=removeFakeNode1(node,nextnodeInfo,bioTree)
for itrace=1:size(node.Out{1}.traceInfo.pixelIdxList,2)
    sumPixelList=[];
    sumMeasurment=[];
    for iNodeOut=1:size(node.Out,2)
        sumPixelList=[sumPixelList;node.Out{iNodeOut}.traceInfo.pixelIdxList{itrace}];
        sumMeasurment=[sumMeasurment;node.Out{iNodeOut}.traceInfo.measurment{itrace}];
    end
    outPixelIdxList{itrace}=sumPixelList;
    outMeasurment{itrace}=sumMeasurment;
end
for iNodeIn=1:size(node.In,2)
    if iNodeIn==1
        if node.In{iNodeIn}.isNode==false
            rootInfo=node.In{iNodeIn}.rootInfo;
            bioTree{nextnodeInfo(1)}.node{nextnodeInfo(2)}.In{nextnodeInfo(3)}.isNode=false;
            bioTree{nextnodeInfo(1)}.node{nextnodeInfo(2)}.In{nextnodeInfo(3)}.rootInfo=rootInfo;
            bioTree{nextnodeInfo(1)}.node{nextnodeInfo(2)}.In{nextnodeInfo(3)}.nodeInfo=[];
            bioTree{rootInfo(1)}.root{rootInfo(2)}.nodeInfo=[nextnodeInfo(1),nextnodeInfo(2),nextnodeInfo(3)];
            bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.pixelIdxList(end+1:end+size(outPixelIdxList,2))=outPixelIdxList;
            bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.measurment(end+1:end+size(outPixelIdxList,2))=outMeasurment;
            continue;
        end
        if node.In{iNodeIn}.isNode==true
            nodeInfo=node.In{iNodeIn}.nodeInfo;            
            bioTree{nextnodeInfo(1)}.node{nextnodeInfo(2)}.In{nextnodeInfo(3)}.nodeInfo=nodeInfo;
            bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.nodeInfo=[nextnodeInfo(1),nextnodeInfo(2),nextnodeInfo(3)];
            bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.pixelIdxList(end+1:end+size(outPixelIdxList,2))=outPixelIdxList;
            bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.measurment(end+1:end+size(outPixelIdxList,2))=outMeasurment;
            continue;
        end
    end
    if iNodeIn>1
        if node.In{iNodeIn}.isNode==false
            rootInfo=node.In{iNodeIn}.rootInfo;
            bioTree{nextnodeInfo(1)}.node{nextnodeInfo(2)}.In{end+1}.isNode=false;
            bioTree{nextnodeInfo(1)}.node{nextnodeInfo(2)}.In{end}.rootInfo=rootInfo;
            bioTree{rootInfo(1)}.root{rootInfo(2)}.nodeInfo=[nextnodeInfo(1),nextnodeInfo(2),size(bioTree{nextnodeInfo(1)}.node{nextnodeInfo(2)}.In,2)];
            bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.pixelIdxList(end+1:end+size(outPixelIdxList,2))=outPixelIdxList;
            bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.measurment(end+1:end+size(outPixelIdxList,2))=outMeasurment;
            continue;
        end
        if node.In{iNodeIn}.isNode==true
            nodeInfo=node.In{iNodeIn}.nodeInfo;
            bioTree{nextnodeInfo(1)}.node{nextnodeInfo(2)}.In{end+1}.isNode=true;
            bioTree{nextnodeInfo(1)}.node{nextnodeInfo(2)}.In{end}.nodeInfo=nodeInfo;
            bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.nodeInfo=[nextnodeInfo(1),nextnodeInfo(2),size(bioTree{nextnodeInfo(1)}.node{nextnodeInfo(2)}.In,2)];
            bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.pixelIdxList(end+1:end+size(outPixelIdxList,2))=outPixelIdxList;
            bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.measurment(end+1:end+size(outPixelIdxList,2))=outMeasurment;
            continue;
        end
    end
end
end
function bioTree=removeFakeNode2(node,nextnodeInfo,bioTree)

for iNodeOut=1:size(node.Out,2)
    outPixelIdxList=node.Out{iNodeOut}.traceInfo.pixelIdxList;
    outMeasurment=node.Out{iNodeOut}.traceInfo.measurment;
    if iNodeOut==1
        if node.Out{iNodeOut}.is2Node==false
            leafInfo=node.Out{iNodeOut}.leafInfo;
            bioTree{nextnodeInfo(1)}.node{nextnodeInfo(2)}.Out{nextnodeInfo(3)}.is2Node=false;
            bioTree{nextnodeInfo(1)}.node{nextnodeInfo(2)}.Out{nextnodeInfo(3)}.leafInfo=leafInfo;
            bioTree{nextnodeInfo(1)}.node{nextnodeInfo(2)}.Out{nextnodeInfo(3)}.nodeInfo=[];
            bioTree{leafInfo(1)}.leavies{leafInfo(2)}.nodeInfo=[nextnodeInfo(1),nextnodeInfo(2),nextnodeInfo(3)];
            bioTree{nextnodeInfo(1)}.node{nextnodeInfo(2)}.Out{nextnodeInfo(3)}.traceInfo.pixelIdxList(end+1:end+size(outPixelIdxList,2))=outPixelIdxList;
            bioTree{nextnodeInfo(1)}.node{nextnodeInfo(2)}.Out{nextnodeInfo(3)}.traceInfo.measurment(end+1:end+size(outPixelIdxList,2))=outMeasurment;
            continue;
        end
        if node.Out{iNodeOut}.is2Node==true
            nodeInfo=node.Out{iNodeOut}.nodeInfo;
            bioTree{nextnodeInfo(1)}.node{nextnodeInfo(2)}.Out{nextnodeInfo(3)}.nodeInfo=nodeInfo;
            bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{nodeInfo(3)}.nodeInfo=[nextnodeInfo(1),nextnodeInfo(2),nextnodeInfo(3)];
            bioTree{nextnodeInfo(1)}.node{nextnodeInfo(2)}.Out{nextnodeInfo(3)}.traceInfo.pixelIdxList(end+1:end+size(outPixelIdxList,2))=outPixelIdxList;
            bioTree{nextnodeInfo(1)}.node{nextnodeInfo(2)}.Out{nextnodeInfo(3)}.traceInfo.measurment(end+1:end+size(outPixelIdxList,2))=outMeasurment;
            continue;
        end
    end
    if iNodeOut>1
        if node.Out{iNodeOut}.is2Node==false
            outIndex=size(bioTree{nextnodeInfo(1)}.node{nextnodeInfo(2)}.Out,2)+1;
            leafInfo=node.Out{iNodeOut}.leafInfo;
            bioTree{nextnodeInfo(1)}.node{nextnodeInfo(2)}.Out{outIndex}.is2Node=false;
            bioTree{nextnodeInfo(1)}.node{nextnodeInfo(2)}.Out{outIndex}.leafInfo=leafInfo;
            bioTree{leafInfo(1)}.leavies{leafInfo(2)}.nodeInfo=[nextnodeInfo(1),nextnodeInfo(2),outIndex];
            bioTree{nextnodeInfo(1)}.node{nextnodeInfo(2)}.Out{outIndex}.traceInfo.pixelIdxList=bioTree{nextnodeInfo(1)}.node{nextnodeInfo(2)}.Out{nextnodeInfo(3)}.traceInfo.pixelIdxList;
            bioTree{nextnodeInfo(1)}.node{nextnodeInfo(2)}.Out{outIndex}.traceInfo.measurment=bioTree{nextnodeInfo(1)}.node{nextnodeInfo(2)}.Out{nextnodeInfo(3)}.traceInfo.measurment;
            bioTree{nextnodeInfo(1)}.node{nextnodeInfo(2)}.Out{outIndex}.traceInfo.pixelIdxList(end+1:end+size(outPixelIdxList,2))=outPixelIdxList;
            bioTree{nextnodeInfo(1)}.node{nextnodeInfo(2)}.Out{outIndex}.traceInfo.measurment(end+1:end+size(outPixelIdxList,2))=outMeasurment;
            continue;
        end
        if node.Out{iNodeOut}.is2Node==true
            outIndex=size(bioTree{nextnodeInfo(1)}.node{nextnodeInfo(2)}.Out,2)+1;
            nodeInfo=node.Out{iNodeOut}.nodeInfo;
            bioTree{nextnodeInfo(1)}.node{nextnodeInfo(2)}.Out{outIndex}.is2Node=true;
            bioTree{nextnodeInfo(1)}.node{nextnodeInfo(2)}.Out{outIndex}.nodeInfo=nodeInfo;
            bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{nodeInfo(3)}.nodeInfo=[nextnodeInfo(1),nextnodeInfo(2),outIndex];
            bioTree{nextnodeInfo(1)}.node{nextnodeInfo(2)}.Out{outIndex}.traceInfo.pixelIdxList=bioTree{nextnodeInfo(1)}.node{nextnodeInfo(2)}.Out{nextnodeInfo(3)}.traceInfo.pixelIdxList;
            bioTree{nextnodeInfo(1)}.node{nextnodeInfo(2)}.Out{outIndex}.traceInfo.measurment=bioTree{nextnodeInfo(1)}.node{nextnodeInfo(2)}.Out{nextnodeInfo(3)}.traceInfo.measurment;
            bioTree{nextnodeInfo(1)}.node{nextnodeInfo(2)}.Out{outIndex}.traceInfo.pixelIdxList(end+1:end+size(outPixelIdxList,2))=outPixelIdxList;
            bioTree{nextnodeInfo(1)}.node{nextnodeInfo(2)}.Out{outIndex}.traceInfo.measurment(end+1:end+size(outPixelIdxList,2))=outMeasurment;
            continue;
        end
    end
end
end

function bioTree=nodeCorrect(bioTree,iframe)
if ~isempty(bioTree{iframe}.node)
    for iNode=1:size(bioTree{iframe}.node,2)
        for iNodeIn=1:size(bioTree{iframe}.node{iNode}.In,2)
            if bioTree{iframe}.node{iNode}.In{iNodeIn}.isNode==false
                rootInfo=bioTree{iframe}.node{iNode}.In{iNodeIn}.rootInfo;
                bioTree{rootInfo(1)}.root{rootInfo(2)}.nodeInfo=[iframe,iNode,iNodeIn];
            end
            if bioTree{iframe}.node{iNode}.In{iNodeIn}.isNode==true
                nodeInfo=bioTree{iframe}.node{iNode}.In{iNodeIn}.nodeInfo;
                bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.nodeInfo=[iframe,iNode,iNodeIn];
            end
        end
        for iNodeOut=1:size(bioTree{iframe}.node{iNode}.Out,2)
            if bioTree{iframe}.node{iNode}.Out{iNodeOut}.is2Node==false
                leafInfo=bioTree{iframe}.node{iNode}.Out{iNodeOut}.leafInfo;
                bioTree{leafInfo(1)}.leavies{leafInfo(2)}.nodeInfo=[iframe,iNode,iNodeOut];
            end
            if bioTree{iframe}.node{iNode}.Out{iNodeOut}.is2Node==true
                nodeInfo=bioTree{iframe}.node{iNode}.Out{iNodeOut}.nodeInfo;
                bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{nodeInfo(3)}.nodeInfo=[iframe,iNode,iNodeOut];
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