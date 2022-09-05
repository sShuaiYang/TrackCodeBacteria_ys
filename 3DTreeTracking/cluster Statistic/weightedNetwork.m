function clusterTree=weightedNetwork(bacteriaFrameInfo,bioTree,step,maxLength)
for iFrame=1:size(bacteriaFrameInfo,2)
    disp(iFrame)
    averageLength=mean(bacteriaFrameInfo{iFrame}.lengthInfo);
    centroidInfo=bacteriaFrameInfo{iFrame}.centroidInfo;
    bacteriaList=bacteriaFrameInfo{iFrame}.bacteriaInfo;
    D=pdist2(centroidInfo,centroidInfo,'euclidean');
    D=double((D<=(averageLength*1.5)));
    for i=1:size(D,1)
        D(i,i)=0;
    end
    currentTree.distMatrix=D;
    currentTree.bacteriaList=bacteriaList;
    if iFrame~=1
        currentTree=makeNewDistMatrix(preTree,currentTree,bioTree,iFrame,maxLength);
    end
    branchList=currentTree.bacteriaList(:,5);
    [branchList,branchOrder]=sort(branchList);
    currentTree.branchOrder=cat(2,branchList,branchOrder);
    clusterTree{iFrame}.branchOrder=cat(2,branchList,branchOrder);
%     clusterCoeficient=caculateNetwork(currentTree.branchOrder,currentTree.distMatrix);
%     clusterTree{iFrame}.clusterCoeficient=mean(clusterCoeficient);
%     clusterTree{iFrame}.averageDegree=numel(currentTree.distMatrix(currentTree.distMatrix==1))/size(currentTree.distMatrix,1);
    if iFrame>=200 && (mod(iFrame,step)==1 || iFrame==size(bacteriaFrameInfo,2))
        clusterTree{iFrame}.distMatrix=sparse(currentTree.distMatrix);
        clusterTree{iFrame}.bacteriaList=currentTree.bacteriaList;
%         if max(max(full(clusterTree{iFrame}.distMatrix)))==1
%             [maxArea,matrixSize,averageDistance]=getMaxCluster(full(clusterTree{iFrame}.distMatrix));
%         else
%             maxArea=0;
%             matrixSize=0;
%             averageDistance=inf;
%         end
%         clusterTree{iFrame}.networkSizeInfo=[maxArea,matrixSize,maxArea/matrixSize];
%         clusterTree{iFrame}.averageDistance=averageDistance;
    end
    preTree=currentTree;
end
end
function newClusterTree=makeNewDistMatrix(preTree,currentTree,bioTree,iFrame,maxLength)
% maxLength=600;
newClusterTree=preTree;
correspondList=zeros(size(currentTree.bacteriaList,1),1);
currentOrderNum=1:size(preTree.bacteriaList,1);
for i=1:size(currentTree.bacteriaList,1)
    bacteriaInfo=currentTree.bacteriaList(i,:);
    if bacteriaInfo(4)~=1
        ibacteria=currentOrderNum(bacteriaInfo(1)==preTree.bacteriaList(:,1) & bacteriaInfo(2)==preTree.bacteriaList(:,2) & bacteriaInfo(3)==preTree.bacteriaList(:,3));
        newClusterTree.bacteriaList(ibacteria,:)=currentTree.bacteriaList(i,:);
        correspondList(i)=ibacteria;
    end
    if bacteriaInfo(4)==1
        if bacteriaInfo(3)==0
            if bioTree{bacteriaInfo(1)}.root{bacteriaInfo(2)}.is2Node==1
                nextNode=bioTree{bacteriaInfo(1)}.root{bacteriaInfo(2)}.nodeInfo;
                if nextNode(1)==iFrame
                    continue
                end
            end
            newClusterTree.bacteriaList(end+1,:)=currentTree.bacteriaList(i,:);
            newClusterTree.distMatrix(size(newClusterTree.bacteriaList,1),:)=0;
            newClusterTree.distMatrix(:,size(newClusterTree.bacteriaList,1))=0;
            correspondList(i)=size(newClusterTree.bacteriaList,1);
        end
        if bacteriaInfo(3)~=0
            nodeInfo=bacteriaInfo(1:2);
            if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{1}.isNode==0
                preRoot=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{1}.rootInfo;
                if preRoot(1)==iFrame
                    totalNum=0;
                    preParent=[];
                    for j=1:size(currentTree.bacteriaList,1)
                        if currentTree.bacteriaList(j,3)~=0
                            if isequal(bacteriaInfo(1:2),currentTree.bacteriaList(j,1:2))
                                totalNum=totalNum+1;
                                preParent=[preParent;j];
                                if bacteriaInfo(3)==currentTree.bacteriaList(j,3)
                                    preParent(end)=[];
                                    break
                                end
                            end
                        end
                    end
                    newClusterTree.bacteriaList(end+1,:)=currentTree.bacteriaList(i,:);
                    newClusterTree.distMatrix(size(newClusterTree.bacteriaList,1),:)=0;
                    newClusterTree.distMatrix(:,size(newClusterTree.bacteriaList,1))=0;
                    newOrder=size(newClusterTree.bacteriaList,1);
                    correspondList(i)=newOrder;
                    for iParent=1:numel(preParent)
                        if currentTree.distMatrix(i,preParent(iParent))==1
                            newClusterTree.distMatrix(newOrder,correspondList(preParent(iParent)))=1;
                            newClusterTree.distMatrix(correspondList(preParent(iParent)),newOrder)=1;
                        end
                    end
                    continue
                end
            end
            totalNum=0;
            inNum=0;
            preParent=[];
            for j=1:size(currentTree.bacteriaList,1)
                if currentTree.bacteriaList(j,3)~=0
                    if isequal(bacteriaInfo(1:2),currentTree.bacteriaList(j,1:2))
                        totalNum=totalNum+1;
                        preParent=[preParent;j];
                        if bacteriaInfo(3)==currentTree.bacteriaList(j,3)
                            inNum=totalNum;
                            preParent(end)=[];
                            break
                        end
                    end
                end
            end
            if inNum==1
                nodeInfo=bacteriaInfo(1:2);
                if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{1}.isNode==0
                    preRoot=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{1}.rootInfo;
                    preRoot=[preRoot,0];
                    ibacteria=currentOrderNum(preRoot(1)==preTree.bacteriaList(:,1) & preRoot(2)==preTree.bacteriaList(:,2) & preRoot(3)==preTree.bacteriaList(:,3));
                    newClusterTree.bacteriaList(ibacteria,:)=currentTree.bacteriaList(i,:);
                    newClusterTree.distMatrix(ibacteria,:)=newClusterTree.distMatrix(ibacteria,:);
                    newClusterTree.distMatrix(:,ibacteria)=newClusterTree.distMatrix(:,ibacteria);
                    correspondList(i)=ibacteria;
                end
                if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{1}.isNode==1
                    preNode=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{1}.nodeInfo;
                    ibacteria=currentOrderNum(preNode(1)==preTree.bacteriaList(:,1) & preNode(2)==preTree.bacteriaList(:,2) & preNode(3)==preTree.bacteriaList(:,3));
                    newClusterTree.bacteriaList(ibacteria,:)=currentTree.bacteriaList(i,:);
                    newClusterTree.distMatrix(ibacteria,:)=newClusterTree.distMatrix(ibacteria,:);
                    newClusterTree.distMatrix(:,ibacteria)=newClusterTree.distMatrix(:,ibacteria);
                    correspondList(i)=ibacteria;
                end
            end
            if inNum>1
                newClusterTree.bacteriaList(end+1,:)=currentTree.bacteriaList(i,:);
                newClusterTree.distMatrix(size(newClusterTree.bacteriaList,1),:)=0;
                newClusterTree.distMatrix(:,size(newClusterTree.bacteriaList,1))=0;
                nodeInfo=bacteriaInfo(1:2);
                if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{1}.isNode==0
                    preRoot=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{1}.rootInfo;
                    preRoot=[preRoot,0];
                    ibacteria=currentOrderNum(preRoot(1)==preTree.bacteriaList(:,1) & preRoot(2)==preTree.bacteriaList(:,2) & preRoot(3)==preTree.bacteriaList(:,3));
                end
                if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{1}.isNode==1
                    preNode=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{1}.nodeInfo;
                    ibacteria=currentOrderNum(preNode(1)==preTree.bacteriaList(:,1) & preNode(2)==preTree.bacteriaList(:,2) & preNode(3)==preTree.bacteriaList(:,3));
                end
                linkInfo=preTree.distMatrix(ibacteria,:);
                newClusterTree.distMatrix(1:size(linkInfo,2),end)=linkInfo';
                newClusterTree.distMatrix(end,1:size(linkInfo,2))=linkInfo;
                newOrder=size(newClusterTree.bacteriaList,1);
                correspondList(i)=newOrder;
                for iParent=1:numel(preParent)
                    if currentTree.distMatrix(i,preParent(iParent))==1
                        newClusterTree.distMatrix(newOrder,correspondList(preParent(iParent)))=100000;
                        newClusterTree.distMatrix(correspondList(preParent(iParent)),newOrder)=100000;
                    end
                end
            end
        end
    end
end
newClusterTree.distMatrix(correspondList,correspondList)=currentTree.distMatrix+newClusterTree.distMatrix(correspondList,correspondList);
% newClusterTree.distMatrix(newClusterTree.distMatrix>1)=1;
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
function [maxArea,matrixSize,averageDistance]=getMaxCluster(distMatrix)
gObj=biograph(distMatrix);
matrixSize=size(distMatrix,1);
[~,C]=conncomp(gObj,'directed','false');
c=sort(C);
cc=regionprops(c,'FilledArea');
maxArea=0;
for i=1:size(cc,1)
    if cc(i).FilledArea>maxArea
        maxArea=cc(i).FilledArea;
        maxNum=i;
    end
end
if maxArea>=100
    newDist=allshortestpaths(gObj);
    averageDistance=mean(newDist(newDist~=inf));
else
    averageDistance=+inf;
end
end