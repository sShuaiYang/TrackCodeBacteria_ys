function bioTree=changeBranchNum(bioTree)
focusBranchNum=[10,14,10000];
for iframe=1:size(bioTree,2)
    for iNode=1:size(bioTree{iframe}.node,2)
        nodeInfo=[iframe,iNode];
        bioTree=needChangeNode(bioTree,nodeInfo,focusBranchNum);
    end
end
end
function bioTree=needChangeNode(bioTree,nodeInfo,focusBranchNum)
inIndex=[];
for iIn=1:size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In,2)
    if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.isNode==1
        preNode=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.nodeInfo;
        inIndex=[inIndex;bioTree{preNode(1)}.node{preNode(2)}.branchIndex];
    else
        preRoot=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.rootInfo;
        inIndex=[inIndex;bioTree{preRoot(1)}.root{preRoot(2)}.branchIndex];
    end
end
inIndex=unique(inIndex);
if numel(inIndex)>=2 || (numel(inIndex)==1 && inIndex==10000)
    if any(ismember(inIndex,focusBranchNum))
        outIndex=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.branchIndex;
        if ismember(outIndex,focusBranchNum)
            bioTree{nodeInfo(1)}.node{nodeInfo(2)}.branchIndex=10000;
        end
    end
end
end
