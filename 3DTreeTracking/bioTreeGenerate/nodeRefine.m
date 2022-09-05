function [bioTree,nodeNumber]=nodeRefine(bioTree,frameShift) %this function use to remove 1 In to 1 Out
nodeNumber0=nodeCounter(bioTree,frameShift);
% disp(strcat('orignal node number=',num2str(nodeNumber)));
for iframe=1+frameShift:size(bioTree,2)
%     disp(iframe);
%     if iframe==860
%         a=1;
%     end
    countNode=0;
    nodeClear=[];
    if ~isempty(bioTree{iframe}.node)
        for iNode=1:size(bioTree{iframe}.node,2)
            [nodetest,nextnodeInfo,nextleafInfo]=isTrueNode(bioTree{iframe}.node{iNode});
            if nodetest==false
                bioTree=removeFakeNode(bioTree{iframe}.node{iNode},nextnodeInfo,nextleafInfo,bioTree);
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
% disp(strcat('afterRefine node number=',num2str(nodeNumber)));
nodeNumber=[nodeNumber0,nodeNumber1];
end
function [nodetest,nextnodeInfo,nextleafInfo]=isTrueNode(node)
if size(node.Out,2)==1 && size(node.In,2)==1 
    nodetest=false;
    if node.Out{1}.is2Node==false
        nextleafInfo=node.Out{1}.leafInfo;
        nextnodeInfo=[];
        return;
    end
    if node.Out{1}.is2Node==true
        nextleafInfo=[];
        nextnodeInfo=node.Out{1}.nodeInfo;
        return;
    end
else
   nodetest=true;
   nextleafInfo=[];
   nextnodeInfo=[];             
end
end
function bioTree=removeFakeNode(node,nextnodeInfo,nextleafInfo,bioTree)
outPixelIdxList=node.Out{1}.traceInfo.pixelIdxList;
outMeasurment=node.Out{1}.traceInfo.measurment;
if node.In{1}.isNode==false
    rootInfo=node.In{1}.rootInfo;
    if node.Out{1}.is2Node==false
        bioTree{nextleafInfo(1)}.leavies{nextleafInfo(2)}.is2Node=false;
        bioTree{nextleafInfo(1)}.leavies{nextleafInfo(2)}.rootInfo=rootInfo;
        bioTree{nextleafInfo(1)}.leavies{nextleafInfo(2)}.nodeInfo=[];
        bioTree{rootInfo(1)}.root{rootInfo(2)}.is2Node=false;
        bioTree{rootInfo(1)}.root{rootInfo(2)}.nodeInfo=[];
        bioTree{rootInfo(1)}.root{rootInfo(2)}.leafInfo=nextleafInfo;
        bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.pixelIdxList(end+1:end+size(outPixelIdxList,2))=outPixelIdxList;
        bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.measurment(end+1:end+size(outPixelIdxList,2))=outMeasurment;
    end
    if node.Out{1}.is2Node==true
        bioTree{nextnodeInfo(1)}.node{nextnodeInfo(2)}.In{nextnodeInfo(3)}.isNode=false;
        bioTree{nextnodeInfo(1)}.node{nextnodeInfo(2)}.In{nextnodeInfo(3)}.rootInfo=rootInfo;
        bioTree{nextnodeInfo(1)}.node{nextnodeInfo(2)}.In{nextnodeInfo(3)}.nodeInfo=[];
        bioTree{rootInfo(1)}.root{rootInfo(2)}.nodeInfo=nextnodeInfo;
        bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.pixelIdxList(end+1:end+size(outPixelIdxList,2))=outPixelIdxList;
        bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.measurment(end+1:end+size(outPixelIdxList,2))=outMeasurment;
    end
end
if node.In{1}.isNode==true
    nodeInfo=node.In{1}.nodeInfo;
    if node.Out{1}.is2Node==false
        bioTree{nextleafInfo(1)}.leavies{nextleafInfo(2)}.nodeInfo=nodeInfo;
        bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.is2Node=false;
        bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.leafInfo=nextleafInfo;
        bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.nodeInfo=[];
        bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.pixelIdxList(end+1:end+size(outPixelIdxList,2))=outPixelIdxList;
        bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.measurment(end+1:end+size(outPixelIdxList,2))=outMeasurment;
    end
    if node.Out{1}.is2Node==true
        bioTree{nextnodeInfo(1)}.node{nextnodeInfo(2)}.In{nextnodeInfo(3)}.isNode=true;
        bioTree{nextnodeInfo(1)}.node{nextnodeInfo(2)}.In{nextnodeInfo(3)}.rootInfo=[];
        bioTree{nextnodeInfo(1)}.node{nextnodeInfo(2)}.In{nextnodeInfo(3)}.nodeInfo=nodeInfo;
        bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.nodeInfo=nextnodeInfo;
        bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.pixelIdxList(end+1:end+size(outPixelIdxList,2))=outPixelIdxList;
        bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.measurment(end+1:end+size(outPixelIdxList,2))=outMeasurment;
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