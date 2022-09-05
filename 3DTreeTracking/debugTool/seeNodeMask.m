function [pixelIdxListIn,pixelIdxListOut]=seeNodeMask(bioTree,nodeInfo) 
for iIn=1:size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In,2)
    if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.isNode==true
        nodeInfo_pre=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.nodeInfo;
        pixelIdxListIn{iIn}=bioTree{nodeInfo_pre(1)}.node{nodeInfo_pre(2)}.Out{nodeInfo_pre(3)}.traceInfo.pixelIdxList{end};
    end
    if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.isNode==false
        rootInfo=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.rootInfo;
       pixelIdxListIn{iIn}=bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.pixelIdxList{end};
    end
end
for iOut=1:size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out,2)
       pixelIdxListOut{iOut}=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.traceInfo.pixelIdxList;
end
end