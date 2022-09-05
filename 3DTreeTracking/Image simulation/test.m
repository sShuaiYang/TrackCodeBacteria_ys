[oriTree,frameNum,imageSize,dirFile]=bacteriaSimulation(4000);
bioTree=makeBioTreeWithImageStack(dirFile);
%%≤Ω÷Ë
bioTree=autoTidyBioTree(bioTree);
for i=1:4
    bioTree=type1NodeReduction(bioTree,1);
    bioTree=type_OutNodeReduction(bioTree);
    bioTree=type2NodeReduction(bioTree,1);
    bioTree=type3NodeReduction(bioTree,15);
    bioTree=type4NodeReduction(bioTree);
    bioTree=autoTidyBioTree(bioTree);
end
bioTree1=bioTree;
for i=1:4
    bioTree=type_InNodeReduction(bioTree);
    bioTree=type2NodeReduction(bioTree,1);
    bioTree=type3NodeReduction(bioTree,15);
    bioTree=type4NodeReduction(bioTree);
    bioTree=autoTidyBioTree(bioTree);
end
bioTree2=bioTree;
for i=1:4
    bioTree=type1NodeReduction(bioTree,2);
    bioTree=type_InNodeReduction(bioTree);
    bioTree=type2NodeReduction(bioTree,2);
    bioTree=type3NodeReduction(bioTree,30);
    bioTree=type4NodeReduction(bioTree);
    bioTree=autoTidyBioTree(bioTree);
end
bioTree=mergeNode(bioTree);
for i=1:4
    bioTree=type1NodeReduction(bioTree,2);
    bioTree=type_InNodeReduction(bioTree);
    bioTree=type2NodeReduction(bioTree,2);
    bioTree=type3NodeReduction(bioTree,30);
    bioTree=type4NodeReduction(bioTree);
end
[bioTree,branchList,unCertainList,emptyNodeList]=myBiograph_new2(bioTree);
bioTree=bioTreeMeasure(bioTree,0,bioTree{1}.imageSize(1),bioTree{1}.imageSize(2));
[bioTree,branchList,allList]=divisionFinder(bioTree,branchList);
makeBranchDemo(bioTree,branchList,allList,1,3,7000);
bioTree1=generateBioTree(oriTree,frameNum,imageSize);
[bioTree1,branchList1,unCertainList,emptyNodeList]=myBiograph_new2(bioTree1);
[bioTree1,branchList1,allList1]=divisionFinder(bioTree1,branchList1);
for i=1:4
    bioTree=type1NodeReduction(bioTree);
    bioTree=type_OutNodeReduction(bioTree);
    bioTree=type1NodeReduction_deep(bioTree);
    bioTree=type2NodeReduction(bioTree);
    bioTree=type3NodeReduction(bioTree);
    bioTree=type4NodeReduction(bioTree);
end

% protocal
bioTree=batchTreeTrackingJob();
bioTree=autoTidyBioTree(bioTree);
for i=1:4
    bioTree=type1NodeReduction(bioTree,1);
    bioTree=type_OutNodeReduction(bioTree);
    bioTree=type2NodeReduction(bioTree,1);
    bioTree=type3NodeReduction(bioTree,15);
    bioTree=type4NodeReduction(bioTree);
    bioTree=autoTidyBioTree(bioTree);
end
bioTree1=bioTree;
for i=1:4
    bioTree=type_InNodeReduction(bioTree);
    bioTree=type2NodeReduction(bioTree,1);
    bioTree=type3NodeReduction(bioTree,15);
    bioTree=type4NodeReduction(bioTree);
    bioTree=autoTidyBioTree(bioTree);
end
bioTree2=bioTree;
for i=1:4
    bioTree=type1NodeReduction(bioTree,2);
    bioTree=type_InNodeReduction(bioTree);
    bioTree=type2NodeReduction(bioTree,2);
    bioTree=type3NodeReduction(bioTree,30);
    bioTree=type4NodeReduction(bioTree);
    bioTree=autoTidyBioTree(bioTree);
end
bioTree=mergeNode(bioTree);
for i=1:4
    bioTree=type1NodeReduction(bioTree,2);
    bioTree=type_InNodeReduction(bioTree);
    bioTree=type2NodeReduction(bioTree,2);
    bioTree=type3NodeReduction(bioTree,30);
    bioTree=type4NodeReduction(bioTree);
end
[bioTree,branchList,unCertainList,emptyNodeList]=myBiograph_new2(bioTree);
bioTree=bioTreeMeasure(bioTree,0,bioTree{1}.imageSize(1),bioTree{1}.imageSize(2));
[bioTree,branchList,allList]=divisionFinder(bioTree,branchList);
makeBranchDemo(bioTree,branchList,allList,1,3,8000);
