function branchMap=bioTreeBranchDrawing(bioTree,branchNum)
branchMap=zeros(bioTree{1}.imageSize(1),bioTree{1}.imageSize(2),3,'uint8');
map1=branchMap(:,:,1);   % r
map2=branchMap(:,:,2);   % g
map3=branchMap(:,:,3);   % b 
[bioTree,branchList,~,~]=myBiograph_new2(bioTree);
[bioTree,branchList,allList]=divisionFinder(bioTree,branchList);
iBranchInfo=branchList(branchNum,:);
allNode=bioTree{iBranchInfo(1)}.node{iBranchInfo(2)}.allNode;
allRoot=bioTree{iBranchInfo(1)}.node{iBranchInfo(2)}.allRoot;
allLeaf=bioTree{iBranchInfo(1)}.node{iBranchInfo(2)}.allLeaf;
for iRoot=1:size(allRoot,1)
    rootInfo=allRoot(iRoot,:);
    rootMask=bioTree{rootInfo(1)}.root{rootInfo(2)}.rootPixelDetail;
    map1(rootMask)=255;
end
for iNode=1:size(allNode,1)
    nodeInfo=allNode(iNode,:);
    [~,pixelIdxListOut]=seeNodeMask(bioTree,nodeInfo);
    for iOut=1:2
        map2(pixelIdxListOut{iOut}{1})=255;
        map1(pixelIdxListOut{iOut}{1})=255;
    end
end
for iLeaf=1:size(allLeaf,1)
    leafInfo=allLeaf(iLeaf,:);
    leafMask=bioTree{leafInfo(1)}.leavies{leafInfo(2)}.leafPixelDetail;
    map2(leafMask)=255;
end
map1=imfill(map1,'holes');
map2=imfill(map2,'holes');
map3=imfill(map3,'holes');
branchMap=cat(3,map1,map2,map3);
% branchMap=branchMap(1:590,1:1700,:,:);
imshow(branchMap)
hold on
colorAll=colormap(jet(size(allRoot,1)+2*size(allNode,1)));
% colorAll=zeros((size(allRoot,1)+2*size(allNode,1)),3);
colorAll(:,3)=1;
iColor=0;
% plot root trace
for iRoot=1:size(allRoot,1)
    iColor=iColor+1;
    rootInfo=allRoot(iRoot,:);
    rootMeasurment=bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.measurment;
    centroidList=getCentroidFromMeasurment(rootMeasurment);
    plot(centroidList(:,1),centroidList(:,2),'color',colorAll(iColor,:),'linewidth',2);
    text(centroidList(ceil(end/2),1),centroidList(ceil(end/2),2),num2str(size(centroidList,1)),'Color','w')
end
% plot node trace
for iNode=1:size(allNode,1)
    nodeInfo=allNode(iNode,:);
    for iOut=1:2
        iColor=iColor+1;
        nodeMeasurment=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.traceInfo.measurment;
        centroidList=getCentroidFromMeasurment(nodeMeasurment);
        plot(centroidList(:,1),centroidList(:,2),'color',colorAll(iColor,:),'linewidth',2);
%         text(centroidList(ceil(end/2),1),centroidList(ceil(end/2),2),num2str(size(centroidList,1)),'Color','w')
    end
end
end
function centroidList=getCentroidFromMeasurment(measurment)
centroidList=[];
for iTrace=1:size(measurment,2)
    centroidList=[centroidList;measurment{iTrace}(1).Centroid];
end
end
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