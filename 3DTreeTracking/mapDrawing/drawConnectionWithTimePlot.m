function drawConnectionWithTimePlot(bioTree,bacteriaFrameInfo)
timeSeries=1:200:numel(bioTree);
aimBranch=3;
color=colormap(jet(numel(timeSeries)));
for i=timeSeries
    bacteriaInfo=bacteriaFrameInfo{i}.bacteriaInfo;
    centroidInfo=bacteriaFrameInfo{i}.centroidInfo;
    orderNum=1:size(bacteriaInfo,1);
    orderNum=orderNum(bacteriaInfo(:,5)==aimBranch);
    distMatrix=pdist2(centroidInfo,centroidInfo);
    for iMatrix=1:size(distMatrix,1)
        distMatrix(iMatrix,iMatrix)=10000;
    end
    distMatrix(distMatrix<=50*3)=1;
    distMatrix(distMatrix~=1)=0;
    distMatrix=logical(distMatrix);
    allCentroid=[];
    for iOrder=1:numel(orderNum)
        linkBac=distMatrix(orderNum(iOrder),:);
        linkCentroid=centroidInfo(linkBac,:);
        orderCentroid=centroidInfo(orderNum(iOrder),:);
        linkCentroid(:,1)=linkCentroid(:,1)+orderCentroid(1);
        linkCentroid(:,2)=linkCentroid(:,2)+orderCentroid(2);
        centroid=linkCentroid/2;
        allCentroid=[allCentroid;centroid]; 
    end
    hold on;plot(allCentroid(:,1)*0.064,allCentroid(:,2)*0.064,'color',color(ceil(i/200),:),'lineStyle','none','marker','o','markerSize',5)
%     hold on;plot(allCentroid(:,1)*0.064,allCentroid(:,2)*0.064,'color',color(ceil(i/200),:),'lineStyle','none','marker','.','markerSize',10)
end
end
        
    
    