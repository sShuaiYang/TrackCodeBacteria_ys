function clusterTree=shortestDisNetwork(bacteriaFrameInfo)
clusterTree=[];
for iFrame=1:size(bacteriaFrameInfo,2)
    averageLength=mean(bacteriaFrameInfo{iFrame}.lengthInfo);
    centroidInfo=bacteriaFrameInfo{iFrame}.centroidInfo;
    D=pdist2(centroidInfo,centroidInfo,'euclidean');
    D=(D<=(averageLength*1.2));
    for i=1:size(D,1)
        D(i,i)=0;
    end
    branchList=bacteriaFrameInfo{iFrame}.bacteriaInfo(:,5);
    [branchList,branchOrder]=sort(branchList);
    clusterTree{iFrame}.branchOrder=cat(2,branchList,branchOrder);
    clusterTree{iFrame}.distMatrix=D;
    clusterCoeficient=caculateNetwork(clusterTree{iFrame}.branchOrder,D);
    clusterTree{iFrame}.clusterCoeficient=mean(clusterCoeficient);
end
end
function clusterCoeficient=caculateNetwork(branchInfo,distMatrix)
clusterCoeficient=[];
numList=1:size(distMatrix,1);
for i=1:size(distMatrix,1)
    distI=distMatrix(i,:);
    linkNum=numList(distI==1);
    if numel(linkNum)==0
        clusterCoeficient=[clusterCoeficient;0];
    end
    if numel(linkNum)==1
        clusterCoeficient=[clusterCoeficient;0];
    end
    if numel(linkNum)>1
    linkMatrix=distMatrix(linkNum,linkNum);
    edgeNum=numel(linkMatrix(linkMatrix==1));
    clusterCoeficient=[clusterCoeficient;edgeNum/numel(linkNum)/(numel(linkNum)-1)];
    end
end
end
function distDistribution=getDistDistribution(distMatrix)
distDistribution=zeros(100,1);
for i=1:size(distMatrix,1)
    linkD=distMatrix(i,:);
    linkNum=numel(linkD(linkD==1));
    if isempty(distDistribution(linkNum+1,1))
        distDistribution(linkNum+1,1)=1;
    else
        distDistribution(linkNum+1,1)=distDistribution(linkNum+1,1)+1;
    end
end
end