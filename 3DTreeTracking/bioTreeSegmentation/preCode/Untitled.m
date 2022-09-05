 maskImage=bioTree{1,163}.node{1,1}.Out{1,1}.traceInfo.pixelIdxList{1,1};
eachInfo{1}=bioTree{1,1}.root{1,34}.traceInfo.pixelIdxList{1,end};
eachInfo{2}=bioTree{1,1}.root{1,44}.traceInfo.pixelIdxList{1,end};



bioTree9=type1NodeReduction(bioTree9);
bioTree9=type_OutNodeReduction(bioTree9);
bioTree9=type1NodeReduction_deep(bioTree9);
bioTree9=type2NodeReduction(bioTree9);
bioTree9=type_InNodeReduction(bioTree9);
bioTree9=type3NodeReduction(bioTree9);
bioTree9=type4NodeReduction(bioTree9);

bioTree=type1NodeReduction(bioTree);
bioTree=type_OutNodeReduction(bioTree);
bioTree=type1NodeReduction_deep(bioTree);
bioTree=type2NodeReduction(bioTree);
bioTree=type_InNodeReduction(bioTree);
bioTree=type3NodeReduction(bioTree);
bioTree=type4NodeReduction(bioTree);

for i=1:size(pixelIdxListIn,2)
     [regionNum,xyMin,regionImage]=findRegionNum(pixelIdxListIn{i},imageSize);
end
for i=1:size(pixelIdxListOut,2)
    [regionNum,xyMin,regionImage]=findRegionNum(pixelIdxListOut{i}{1},imageSize);
end
bioTree1=type1NodeReduction(bioTree);
bioTree2=type_InNodeReduction(bioTree1);
bioTree3=type_OutNodeReduction(bioTree2);
bioTree4=type1NodeReduction_deep(bioTree3);
bioTree5=type2NodeReduction(bioTree4);
bioTree6=type3NodeReduction(bioTree5);
bioTree7=type4NodeReduction(bioTree6);

p=false(imageSize);
for i=1:numel(pixelIdxListIn)
p(pixelIdxListIn{i})=1;
end
figure;imshow(p)
p=zeros(imageSize);
for i=1:numel(pixelIdxListOut)
p(pixelIdxListOut{i}{1})=1;
end
figure;imshow(p)

p=[];
for i=1:size(result.Out{2}.traceInfo.pixelIdxList,2)
    image=false(imageSize);
    image(result.Out{2}.traceInfo.pixelIdxList{i})=1;
    p(:,:,i)=image;
end
close all
for i=1:inputNum
    figure;imshow(in{i}.BWImage)
end
for i=1:outputNum
    figure;imshow(out{i}.BWImage)
end

imageSize=[1050,1000];
p=false(imageSize);
for i=1:numel(eachInfo)
p(eachInfo{i})=1;
end
imshow(p)


bioTree9=bioTreeMeasure(bioTree9,0,1666,1593);
[bioTree9,branchList,unCertainList,emptyNodeList]=myBiograph_new2(bioTree9);
[bioTree9,branchList,allList]=divisionFinder(bioTree9,branchList);
glueMask=false(1666,1593,5124);
makeBranchDemo(bioTree9,branchList,glueMask,allList,1,1,5124);


bioTree1=type1NodeReduction(bioTree);
bioTree2=type2NodeReduction(bioTree1);
bioTree3=type_InNodeReduction(bioTree2);
bioTree4=type1NodeReduction_deep(bioTree3);
bioTree5=type3NodeReduction(bioTree4);
bioTree6=type4NodeReduction(bioTree5);
bioTree9=type1NodeReduction(bioTree6);
bioTree9=type2NodeReduction(bioTree9);
bioTree9=type_InNodeReduction(bioTree9);
bioTree9=type1NodeReduction_deep(bioTree9);
bioTree9=type3NodeReduction(bioTree9);
bioTree9=type4NodeReduction(bioTree9);
bioTree9=type1NodeReduction(bioTree9);
bioTree9=type2NodeReduction(bioTree9);
bioTree9=type_InNodeReduction(bioTree9);
bioTree9=type1NodeReduction_deep(bioTree9);
bioTree9=type3NodeReduction(bioTree9);
bioTree9=type4NodeReduction(bioTree9);
bioTree9=type1NodeReduction(bioTree9);
bioTree9=type2NodeReduction(bioTree9);
bioTree9=type_InNodeReduction(bioTree9);
bioTree9=type1NodeReduction_deep(bioTree9);
bioTree9=type3NodeReduction(bioTree9);
bioTree9=type4NodeReduction(bioTree9);
bioTree9=bioTreeMeasure(bioTree9,0,1666,1593);
[bioTree9,branchList,unCertainList,emptyNodeList]=myBiograph_new2(bioTree9);
[bioTree9,branchList,allList]=divisionFinder(bioTree9,branchList);
glueMask=false(1666,1593,2563);
makeBranchDemo(bioTree9,branchList,glueMask,allList,1,2,5124);