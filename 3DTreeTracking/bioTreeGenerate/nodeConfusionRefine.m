function bioTree=nodeConfusionRefine(bioTree,frameShift) %this function use to reconnect the nodeIn equal nodeOUT, need to debug.
nodeNumber=nodeCounter(bioTree,frameShift);
disp(strcat('orignal node number=',num2str(nodeNumber)));
for iframe=1+frameShift:size(bioTree,2)
    if ~isempty(bioTree{iframe}.node)
        clearNode=[];
        for iNode=1:size(bioTree{iframe}.node,2)
            [linkIn,linkOut,isRemoveable]=findLinker(bioTree,iframe,iNode);
            if isRemoveable==true
                bioTree=confuseLinker(bioTree,iframe,iNode,linkIn,linkOut);
                clearNode=[clearNode,iNode];
            else
                continue;
            end    
        end
        if ~isempty(clearNode)
            bioTree{iframe}.node(clearNode)=[];
            bioTree=nodeCorrect(bioTree,iframe);
        end
    end
end
nodeNumber=nodeCounter(bioTree,frameShift);
disp(strcat('orignal node number=',num2str(nodeNumber)));
end
function [linkIn,linkOut,isRemoveable]=findLinker(bioTree,iframe,iNode)
linkIn=[];
linkOut=[];
nodeInXY=[];
nodeOutXY=[];
if size(bioTree{iframe}.node{iNode}.In,2)~=size(bioTree{iframe}.node{iNode}.Out,2)
    isRemoveable=false;
    return;
end
for iIn=1:size(bioTree{iframe}.node{iNode}.In,2)
    if bioTree{iframe}.node{iNode}.In{iIn}.isNode==false
        rootInfo=bioTree{iframe}.node{iNode}.In{iIn}.rootInfo;
        for iTrace=1:size(bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.measurment,2)
            if size(bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.measurment{end-iTrace+1},1)==1
                nodeInXY=[nodeInXY;bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.measurment{end-iTrace+1}.Centroid];
                break;
            end
        end
        
    else
        nodeInfo=bioTree{iframe}.node{iNode}.In{iIn}.nodeInfo;
        for iTrace=1:size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.measurment,2)
            if size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.measurment{end-iTrace+1},1)==1
                nodeInXY=[nodeInXY;bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.measurment{end-iTrace+1}.Centroid];
                break;
            end
        end
        
    end
end
for iOut=1:size(bioTree{iframe}.node{iNode}.Out,2)
    for iTrace=1:size(bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.measurment,2)
        if size(bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.measurment{iTrace},1)==1
            nodeOutXY=[nodeOutXY;bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.measurment{iTrace}.Centroid];
            break;
        end
    end
end

if isempty(nodeInXY) || isempty(nodeOutXY)
    isRemoveable=false;
    return;
else
    if size(nodeInXY,1)==size(bioTree{iframe}.node{iNode}.In,2) && size(nodeOutXY,1)==size(bioTree{iframe}.node{iNode}.Out,2)
        isRemoveable=true;
    else
        isRemoveable=false;
        return;
    end
end

[distInOut,indexI]=pdist2(nodeInXY,nodeOutXY,'euclidean','Smallest',1);
while size(distInOut(~isnan(distInOut)),2)>=1
    [~,idxIn]=min(distInOut);
    linkIn=[linkIn,indexI(idxIn)];
    linkOut=[linkOut,idxIn];
    nodeInXY(indexI(idxIn))=NaN;
    nodeOutXY(idxIn)=NaN;
    [distInOut,indexI]=pdist2(nodeInXY,nodeOutXY,'euclidean','Smallest',1);
end
end
function bioTree=confuseLinker(bioTree,iframe,iNode,linkIn,linkOut)
for ilink=1:size(linkIn,2)
    if bioTree{iframe}.node{iNode}.In{linkIn(ilink)}.isNode==false && bioTree{iframe}.node{iNode}.Out{linkOut(ilink)}.is2Node==false
        rootInfo=bioTree{iframe}.node{iNode}.In{linkIn(ilink)}.rootInfo;
        leafInfo=bioTree{iframe}.node{iNode}.Out{linkOut(ilink)}.leafInfo;
        bioTree{rootInfo(1)}.root{rootInfo(2)}.is2Node=false;
        bioTree{rootInfo(1)}.root{rootInfo(2)}.nodeInfo=[];
        bioTree{rootInfo(1)}.root{rootInfo(2)}.leafInfo=leafInfo;
        bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.pixelIdxList(end+1:end+size( bioTree{iframe}.node{iNode}.Out{linkOut(ilink)}.traceInfo.pixelIdxList,2))= bioTree{iframe}.node{iNode}.Out{linkOut(ilink)}.traceInfo.pixelIdxList(1:end);
        bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.measurment(end+1:end+size( bioTree{iframe}.node{iNode}.Out{linkOut(ilink)}.traceInfo.pixelIdxList,2))= bioTree{iframe}.node{iNode}.Out{linkOut(ilink)}.traceInfo.measurment(1:end);
        bioTree{leafInfo(1)}.leavies{leafInfo(2)}.is2Node=false;
        bioTree{leafInfo(1)}.leavies{leafInfo(2)}.nodeInfo=[];
        bioTree{leafInfo(1)}.leavies{leafInfo(2)}.rootInfo=rootInfo;
        continue;
    end
    if bioTree{iframe}.node{iNode}.In{linkIn(ilink)}.isNode==true && bioTree{iframe}.node{iNode}.Out{linkOut(ilink)}.is2Node==false
        nodeInfo=bioTree{iframe}.node{iNode}.In{linkIn(ilink)}.nodeInfo;
        leafInfo=bioTree{iframe}.node{iNode}.Out{linkOut(ilink)}.leafInfo;
        bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.is2Node=false;
        bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.nodeInfo=[];
        bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.leafInfo=leafInfo;
        bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.pixelIdxList(end+1:end+size( bioTree{iframe}.node{iNode}.Out{linkOut(ilink)}.traceInfo.pixelIdxList,2))= bioTree{iframe}.node{iNode}.Out{linkOut(ilink)}.traceInfo.pixelIdxList(1:end);
        bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.measurment(end+1:end+size( bioTree{iframe}.node{iNode}.Out{linkOut(ilink)}.traceInfo.pixelIdxList,2))= bioTree{iframe}.node{iNode}.Out{linkOut(ilink)}.traceInfo.measurment(1:end);
        bioTree{leafInfo(1)}.leavies{leafInfo(2)}.nodeInfo=nodeInfo;
        continue;
    end
    if bioTree{iframe}.node{iNode}.In{linkIn(ilink)}.isNode==false && bioTree{iframe}.node{iNode}.Out{linkOut(ilink)}.is2Node==true
        rootInfo=bioTree{iframe}.node{iNode}.In{linkIn(ilink)}.rootInfo;
        nodeInfo=bioTree{iframe}.node{iNode}.Out{linkOut(ilink)}.nodeInfo;
        bioTree{rootInfo(1)}.root{rootInfo(2)}.nodeInfo=nodeInfo;
        bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.pixelIdxList(end+1:end+size( bioTree{iframe}.node{iNode}.Out{linkOut(ilink)}.traceInfo.pixelIdxList,2))= bioTree{iframe}.node{iNode}.Out{linkOut(ilink)}.traceInfo.pixelIdxList(1:end);
        bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.measurment(end+1:end+size( bioTree{iframe}.node{iNode}.Out{linkOut(ilink)}.traceInfo.pixelIdxList,2))= bioTree{iframe}.node{iNode}.Out{linkOut(ilink)}.traceInfo.measurment(1:end);
        bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{nodeInfo(3)}.isNode=false;
        bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{nodeInfo(3)}.nodeInfo=[];
        bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{nodeInfo(3)}.rootInfo=rootInfo;
        continue;
    end
    if bioTree{iframe}.node{iNode}.In{linkIn(ilink)}.isNode==true && bioTree{iframe}.node{iNode}.Out{linkOut(ilink)}.is2Node==true
        nodeInfo1=bioTree{iframe}.node{iNode}.In{linkIn(ilink)}.nodeInfo;
        nodeInfo2=bioTree{iframe}.node{iNode}.Out{linkOut(ilink)}.nodeInfo;
        bioTree{nodeInfo1(1)}.node{nodeInfo1(2)}.Out{nodeInfo1(3)}.nodeInfo=nodeInfo2;
        bioTree{nodeInfo1(1)}.node{nodeInfo1(2)}.Out{nodeInfo1(3)}.traceInfo.pixelIdxList(end+1:end+size( bioTree{iframe}.node{iNode}.Out{linkOut(ilink)}.traceInfo.pixelIdxList,2))= bioTree{iframe}.node{iNode}.Out{linkOut(ilink)}.traceInfo.pixelIdxList(1:end);
        bioTree{nodeInfo1(1)}.node{nodeInfo1(2)}.Out{nodeInfo1(3)}.traceInfo.measurment(end+1:end+size( bioTree{iframe}.node{iNode}.Out{linkOut(ilink)}.traceInfo.pixelIdxList,2))= bioTree{iframe}.node{iNode}.Out{linkOut(ilink)}.traceInfo.measurment(1:end);
        bioTree{nodeInfo2(1)}.node{nodeInfo2(2)}.In{nodeInfo2(3)}.nodeInfo=nodeInfo1;
        continue;
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