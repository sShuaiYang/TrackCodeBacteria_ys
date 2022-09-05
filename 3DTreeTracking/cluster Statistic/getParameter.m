function clusterTree=getParameter(clusterTree)
for iFrame=1:size(clusterTree,2)
clusterCoeficient=caculateNetwork(clusterTree{iFrame}.branchOrder,clusterTree{iFrame}.distMatrix);
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