function clusterTree=getClusterTree(bacteriaFrameInfo,bioTree)
clusterTree=[];
for iFrame=1:size(bacteriaFrameInfo,2)
    averageLength=mean(bacteriaFrameInfo{iFrame}.lengthInfo);
    linkageMatrix=linkage(bacteriaFrameInfo{iFrame}.centroidInfo,'single','euclidean');
    c=cluster(linkageMatrix,'cutoff',averageLength*1.2,'criterion','distance');
    branchList=bacteriaFrameInfo{iFrame}.bacteriaInfo(:,5);
    NBCMatrix=cat(2,(1:numel(branchList))',branchList,c);
    [~,branchOrder]=sort(branchList);
    NBCMatrix=NBCMatrix(branchOrder,:);
    clusterInfo=getCluster(c);
    clusterTree{iFrame}.NBCMatrix=NBCMatrix;
    clusterTree{iFrame}.clusterInfo=getOneCluster(clusterInfo,bacteriaFrameInfo{iFrame},bioTree);
end
end
function clusterInfo=getCluster(c)
clusterInfo=[];
numInfo=(1:size(c,1))';
iPro=0;
for iCluster=1:max(c);
    iClusterInfo=numInfo(c==iCluster);
    if numel(iClusterInfo)>=2
        iPro=iPro+1;
        clusterInfo{iPro}=iClusterInfo';
    end
end
end
function clusterTree=getOneCluster(clusterInfo,bacteriaFrameInfo,bioTree)
for iCluster=1:size(clusterInfo,2)
    bacteriaList=clusterInfo{iCluster};
    clusterTree{iCluster}.bacteriaInfo=bacteriaFrameInfo.bacteriaInfo(bacteriaList,:);
    for ibacteria=1:numel(bacteriaList)
        bacteriaInfo=clusterTree{iCluster}.bacteriaInfo(ibacteria,:);
        if bacteriaInfo(3)==0
            clusterTree{iCluster}.pixelIdxList{ibacteria}=bioTree{bacteriaInfo(1)}.root{bacteriaInfo(2)}.traceInfo.pixelIdxList{bacteriaInfo(4)};
        end
        if bacteriaInfo(3)~=0
            clusterTree{iCluster}.pixelIdxList{ibacteria}=bioTree{bacteriaInfo(1)}.node{bacteriaInfo(2)}.Out{bacteriaInfo(3)}.traceInfo.pixelIdxList{bacteriaInfo(4)};
        end
    end
end
end
function image=plotImage(testImage,testTree)
clusterNum=size(testTree.clusterInfo,2);
mycolor=colormap(jet(clusterNum));
image1=testImage;
image2=testImage;
image3=testImage;
for iCluster=1:clusterNum
    pixelAll=[];
    for i=1:numel(testTree.clusterInfo{iCluster}.pixelIdxList)
        pixelAll=[pixelAll;testTree.clusterInfo{iCluster}.pixelIdxList{i}];
    end
    image1(pixelAll)=mycolor(iCluster,1)*255;
    image2(pixelAll)=mycolor(iCluster,2)*255;
    image3(pixelAll)=mycolor(iCluster,3)*255;
end
image=cat(3,image1,image2,image3);
end
    