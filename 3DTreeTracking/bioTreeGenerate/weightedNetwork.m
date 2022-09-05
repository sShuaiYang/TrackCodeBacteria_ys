function clusterTree=weightedNetwork(bacteriaFrameInfo,bioTree,step)
for iFrame=1:size(bacteriaFrameInfo,2)
    disp(iFrame)
    averageLength=mean(bacteriaFrameInfo{iFrame}.lengthInfo);
    centroidInfo=bacteriaFrameInfo{iFrame}.centroidInfo;
    bacteriaList=bacteriaFrameInfo{iFrame}.bacteriaInfo;
    D=pdist2(centroidInfo,centroidInfo,'euclidean');
    D=(D<=(averageLength*1.2));
    for i=1:size(D,1)
        D(i,i)=0;
    end
    currentTree.distMatrix=D;
    currentTree.bacteriaList=bacteriaList;
    if iFrame~=1
        currentTree=makeNewDistMatrix(preTree,currentTree,bioTree,iFrame);
    end
    branchList=currentTree.bacteriaList(:,5);
    [branchList,branchOrder]=sort(branchList);
    currentTree.branchOrder=cat(2,branchList,branchOrder);
    clusterTree{iFrame}.branchOrder=cat(2,branchList,branchOrder);
    clusterCoeficient=caculateNetwork(currentTree.branchOrder,currentTree.distMatrix);
    clusterTree{iFrame}.clusterCoeficient=mean(clusterCoeficient);
    clusterTree{iFrame}.averageDegree=numel(currentTree.distMatrix(currentTree.distMatrix==1))/size(currentTree.distMatrix,1);
    if iFrame>=200 && (mod(iFrame,step)==1 || iFrame==size(bacteriaFrameInfo,2))
        clusterTree{iFrame}.distMatrix=sparse(currentTree.distMatrix);
        clusterTree{iFrame}.bacteriaList=currentTree.bacteriaList;
        if max(max(full(clusterTree{iFrame}.distMatrix)))==1
            [maxArea,matrixSize,averageDistance]=getMaxCluster(full(clusterTree{iFrame}.distMatrix));
        else
            maxArea=0;
            matrixSize=0;
            averageDistance=inf;
        end
        clusterTree{iFrame}.networkSizeInfo=[maxArea,matrixSize,maxArea/matrixSize];
        clusterTree{iFrame}.averageDistance=averageDistance;
    end
    preTree=currentTree;
end
end
% function clusterTree=accumulateNetwork(bacteriaFrameInfo,bioTree)
% for iFrame=1:size(bacteriaFrameInfo,2)
%     disp(iFrame)
%     averageLength=mean(bacteriaFrameInfo{iFrame}.lengthInfo);
%     centroidInfo=bacteriaFrameInfo{iFrame}.centroidInfo;
%     bacteriaList=bacteriaFrameInfo{iFrame}.bacteriaInfo;
%     D=pdist2(centroidInfo,centroidInfo,'euclidean');
%     D=(D<=(averageLength*1.2));
%     for i=1:size(D,1)
%         D(i,i)=0;
%     end   
%     clusterTree{iFrame}.distMatrix=D;
%     clusterTree{iFrame}.bacteriaList=bacteriaList;
%     if iFrame~=1
%         clusterTree{iFrame}=makeNewDistMatrix(clusterTree{iFrame-1},clusterTree{iFrame},bioTree);
%     end
%     branchList=clusterTree{iFrame}.bacteriaList(:,5);
%     [branchList,branchOrder]=sort(branchList); 
%     clusterTree{iFrame}.branchOrder=cat(2,branchList,branchOrder);
%     clusterCoeficient=caculateNetwork(clusterTree{iFrame}.branchOrder,clusterTree{iFrame}.distMatrix);
%     clusterTree{iFrame}.clusterCoeficient=mean(clusterCoeficient);
%     clusterTree{iFrame}.averageDegree=numel(clusterTree{iFrame}.distMatrix(clusterTree{iFrame}.distMatrix==1))/size(clusterTree{iFrame}.distMatrix,1)/2;    
% end
% end
function newClusterTree=makeNewDistMatrix(preTree,currentTree,bioTree,iFrame)
newClusterTree=preTree;
correspondList=[];
for i=1:size(currentTree.bacteriaList,1)
    bacteriaInfo=currentTree.bacteriaList(i,:);
    if bacteriaInfo(4)~=1
        for ibacteria=1:size(preTree.bacteriaList,1)
            if isequal(bacteriaInfo(1:3),preTree.bacteriaList(ibacteria,1:3))
                newClusterTree.bacteriaList(ibacteria,:)=currentTree.bacteriaList(i,:);
                correspondList=[correspondList,ibacteria];
                break
            end
        end
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
            correspondList=[correspondList,size(newClusterTree.bacteriaList,1)];
        end
        if bacteriaInfo(3)~=0
            nodeInfo=bacteriaInfo(1:2);
            if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{1}.isNode==0
                preRoot=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{1}.rootInfo;
                if preRoot(1)==iFrame
                    newClusterTree.bacteriaList(end+1,:)=currentTree.bacteriaList(i,:);
                    newClusterTree.distMatrix(size(newClusterTree.bacteriaList,1),:)=0;
                    newClusterTree.distMatrix(:,size(newClusterTree.bacteriaList,1))=0;
                    correspondList=[correspondList,size(newClusterTree.bacteriaList,1)];
                    continue
                end
            end
            totalNum=0;
            inNum=0;
            for j=1:size(currentTree.bacteriaList,1)
                if currentTree.bacteriaList(j,3)~=0
                    if isequal(bacteriaInfo(1:2),currentTree.bacteriaList(j,1:2))
                        totalNum=totalNum+1;
                        if bacteriaInfo(3)==currentTree.bacteriaList(j,3)
                            inNum=totalNum;
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
                    for ibacteria=1:size(preTree.bacteriaList,1)
                        if isequal(preRoot,preTree.bacteriaList(ibacteria,1:3))
                            newClusterTree.bacteriaList(ibacteria,:)=currentTree.bacteriaList(i,:);
                            correspondList=[correspondList,ibacteria];
                            break
                        end
                    end
                end
                if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{1}.isNode==1
                    preNode=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{1}.nodeInfo;
                    for ibacteria=1:size(preTree.bacteriaList,1)
                        if isequal(preNode,preTree.bacteriaList(ibacteria,1:3))
                            newClusterTree.bacteriaList(ibacteria,:)=currentTree.bacteriaList(i,:);
                            correspondList=[correspondList,ibacteria];
                            break
                        end
                    end
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
                    for ibacteria=1:size(preTree.bacteriaList,1)
                        if isequal(preRoot,preTree.bacteriaList(ibacteria,1:3))
                            break
                        end
                    end
                end
                if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{1}.isNode==1
                    preNode=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{1}.nodeInfo;
                    for ibacteria=1:size(preTree.bacteriaList,1)
                        if isequal(preNode,preTree.bacteriaList(ibacteria,1:3))
                            break
                        end
                    end
                end
                linkInfo=preTree.distMatrix(ibacteria,:);
                newClusterTree.distMatrix(1:size(linkInfo,2),end)=linkInfo';
                newClusterTree.distMatrix(end,1:size(linkInfo,2))=linkInfo;
                correspondList=[correspondList,size(newClusterTree.bacteriaList,1)];
            end
        end
    end
end
newClusterTree.distMatrix(correspondList,correspondList)=currentTree.distMatrix | newClusterTree.distMatrix(correspondList,correspondList);
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