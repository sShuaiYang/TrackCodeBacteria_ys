function ratio=findtheRightRatio(bioTree,branchList)
rightNum=0;
allBranch=size(branchList,1);
for i=1:allBranch
    nodeInfo=branchList(i,:);
    allNode=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.allNode;
    allRoot=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.allRoot;
    if size(allRoot,1)==1
        isHyperNode=allNode(:,5);
        if max(isHyperNode)==0
            rightNum=rightNum+1;
        end
    end
end
ratio=rightNum/allBranch;
end
            