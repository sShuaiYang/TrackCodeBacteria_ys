function bioTree=mergeNode(bioTree)
for iFrame=1:size(bioTree,2)
    for iNode=1:size(bioTree{iFrame}.node,2)
        [needMerge,iOutNum]=findProperNode(bioTree{iFrame}.node{iNode});
        if needMerge==1
            for iPage=1:size(bioTree{iFrame}.node{iNode}.Out{1}.traceInfo.pixelIdxList,2)
                for iOut=2:size(bioTree{iFrame}.node{iNode}.Out,2)
                    bioTree{iFrame}.node{iNode}.Out{1}.traceInfo.pixelIdxList{iPage}=[bioTree{iFrame}.node{iNode}.Out{1}.traceInfo.pixelIdxList{iPage};bioTree{iFrame}.node{iNode}.Out{iOut}.traceInfo.pixelIdxList{iPage}];
                end
            end
            bioTree{iFrame}.node{iNode}.Out{1}.nodeInfo(end)=1;
            bioTree{iFrame}.node{iNode}.Out(2:end)=[];
            minOutNum=min(iOutNum);
            iOutNum=iOutNum(iOutNum~=minOutNum);
            bioTree{iFrame}.node{iNode}.Out{1}.nodeInfo(3)=minOutNum;
            newNodeInfo=bioTree{iFrame}.node{iNode}.Out{1}.nodeInfo;
            bioTree{newNodeInfo(1)}.node{newNodeInfo(2)}.In{newNodeInfo(3)}.nodeInfo(end)=1;
            bioTree{newNodeInfo(1)}.node{newNodeInfo(2)}.In(iOutNum)=[];
            for iIn=minOutNum+1:size(bioTree{newNodeInfo(1)}.node{newNodeInfo(2)}.In,2)
                if bioTree{newNodeInfo(1)}.node{newNodeInfo(2)}.In{iIn}.isNode==0
                    rootInfo=bioTree{newNodeInfo(1)}.node{newNodeInfo(2)}.In{iIn}.rootInfo;
                    bioTree{rootInfo(1)}.root{rootInfo(2)}.nodeInfo(3)=iIn;
                end
                if bioTree{newNodeInfo(1)}.node{newNodeInfo(2)}.In{iIn}.isNode==1
                    nodeInfo=bioTree{newNodeInfo(1)}.node{newNodeInfo(2)}.In{iIn}.nodeInfo;
                    bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.nodeInfo(3)=iIn;
                end
            end
        end
    end
end
end
function [needMerge,iOutNum]=findProperNode(bioTreeNode)
needMerge=0;
iOutNum=[];
if size(bioTreeNode.In,2)>=4
    outNum=size(bioTreeNode.Out,2);
    if outNum==1;
        return
    end
    for iOut=1:size(bioTreeNode.Out,2)
        if bioTreeNode.Out{iOut}.is2Node==0
            return
        else
            if iOut==1
                nodeInfo=bioTreeNode.Out{iOut}.nodeInfo;
                iOutNum=nodeInfo(3);
            else
                newNodeInfo=bioTreeNode.Out{iOut}.nodeInfo;
                iOutNum=[iOutNum,newNodeInfo(3)];
                if ~isequal(newNodeInfo(1:2),nodeInfo(1:2))
                    return
                end
            end
        end
    end
    needMerge=1;
end
end