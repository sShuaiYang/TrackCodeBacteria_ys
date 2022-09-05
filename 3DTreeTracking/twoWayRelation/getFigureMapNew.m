% 在原来的基础上改进了作图，可以作两种图，100%正确的就全部绘出，非100%正确的就框起
function getFigureMapNew(bioTree,bacteriaFrameInfo)
% 初始化全部的node 并进行标记谁是一定正确的
sizeNeed=1200*13.5;
bioTree=bioTreeCut(bioTree,sizeNeed);
bacteriaFrameInfo=bacteriaFrameInfo(1:sizeNeed);

minDivisionTime=400;  %% original 800
[bioTree,branchList,~,~]=myBiograph_new2(bioTree);
coreBranch=branchList(:,3)==1;
branchList=branchList(coreBranch,:);
for ibac=1:size(branchList,1)
    iList=branchList(ibac,:);
    allNode=bioTree{iList(1)}.node{iList(2)}.allNode;
    for iNode=1:size(allNode,1)
        allNode(iNode,4)=morphCorrection(bioTree,allNode(iNode,:));
        bioTree{allNode(iNode,1)}.node{allNode(iNode,2)}.divisionNode=allNode(iNode,4);
    end
    allNode=allNode(allNode(:,4)==1,:);
    for iNode=1:size(allNode,1)
        nodeInfo=allNode(iNode,1:2);
        if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{1}.isNode==0
            continue
        end
        preNode=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{1}.nodeInfo;
        if bioTree{preNode(1)}.node{preNode(2)}.divisionNode==0
            continue
        end
        if nodeInfo(1)-preNode(1)<=minDivisionTime
            allNode(iNode,4)=0;
            bioTree{nodeInfo(1)}.node{nodeInfo(2)}.divisionNode=0;
        end
    end
    bioTree{iList(1)}.node{iList(2)}.allNode=allNode; %将是不是divisionNode的信息放在allNode的第四列
end

% aveLength=40;
finalInfo=bacteriaFrameInfo{end};
aveLength=mean(finalInfo.lengthInfo(finalInfo.lengthInfo<=80));
multi=3;
% [bioTree,branchList,~,~]=myBiograph_new2(bioTree);
isCoreBranch=branchList(:,3);
coreBranchNum=numel(isCoreBranch(isCoreBranch==1));
branchNum=coreBranchNum;
colorBar=colormap(jet(coreBranchNum+1));
close all
figure1=figure;
treeSize=numel(bioTree);
if coreBranchNum<=7
    branchNum=coreBranchNum;
else
    branchNum=7;
end
colorBar=colormap(jet(branchNum+1));

% down map
% num=1;
% axes('Parent',figure1,'Units','inches','Position',[7 0.5 7 1.5],'xLim',[0,numel(bioTree)/1200],'box','on');
% leafAll=[];
% for i=1:branchNum
%     if i<=coreBranchNum
%         color=colorBar(i,:);
%     else
%         color=[0,0,1];
%     end
%     [linkMatrix,centroidInfo,leafList,leafNum,allList,getLeaf,allLeaf,linkTwoPoint]=clusterLikeOneBranchPhytree(bioTree,i,branchList(i,:),branchList);
%     if getLeaf==1
%         [positionInfo,positionNum]=generatePos(leafNum,centroidInfo,numel(bioTree));
%         leafYCoo=num+positionInfo;
%         num=num+positionNum;
%         leafAll=[leafAll;leafList];
%     else
%         leafYCoo=num+(1:size(allLeaf,1));
%         num=num+size(allLeaf,1);
%         leafAll=[leafAll;allLeaf];
%     end
%     hold on;clusterLikeTreePlot(linkMatrix,centroidInfo,leafYCoo,color,allLeaf,getLeaf,allList,leafNum,linkTwoPoint);
% end
% set(gca,'yLim',[1,num]);
% view(gca,[90 -90]);

% leaf map
num=1;
axes('Parent',figure1,'Units','inches','Position',[4.646 2 2.354 7],'xLim',[0,numel(bioTree)/1200],'box','on','YTickLabel',{});
leafAll=[];
for i=1:branchNum
    if i<=coreBranchNum
        color=colorBar(i,:);
    else
        color=[0,0,1];
    end;
    [linkMatrix,centroidInfo,leafList,leafNum,allList,getLeaf,allLeaf,linkTwoPoint]=clusterLikeOneBranchPhytree(bioTree,i,branchList(i,:),branchList);
    if getLeaf==1
        [positionInfo,positionNum]=generatePos(leafNum,centroidInfo,numel(bioTree));
        leafYCoo=num+positionInfo;
        num=num+positionNum;
        leafAll=[leafAll;leafList];
    else
        leafYCoo=num+(1:size(allLeaf,1));
        num=num+size(allLeaf,1);
        leafAll=[leafAll;allLeaf];
    end
    hold on;clusterLikeTreePlot(linkMatrix,centroidInfo,leafYCoo,color,allLeaf,getLeaf,allList,leafNum,linkTwoPoint);
end
set(gca,'yLim',[1,num]);

% middle map
axes('Parent',figure1,'Units','inches','Position',[7 2 7 7],'Color',[0,0,0],'box','on','XTickLabel',{},'YTickLabel',{});
leafNum=numel(leafAll(:,1));
aimOrder=leafAll(:,1)==size(bioTree,2);
distMatrix=zeros(leafNum);
orderNum=1:leafNum;
orderNum=orderNum(aimOrder);
centroidInfo=zeros(orderNum,2);
for i=1:numel(orderNum)
    leafInfo=leafAll(orderNum(i),:);
    centroidInfo(i,:)=bioTree{leafInfo(1)}.leavies{leafInfo(2)}.leafMeasurment.Centroid;
    branchIndex(i)=bioTree{leafInfo(1)}.leavies{leafInfo(2)}.branchIndex;
end
aimDist=pdist2(centroidInfo,centroidInfo);
aimDist(aimDist<=multi*aveLength)=1;
% for i=1:numel(orderNum)
%     aimDist(i,i)=0;
% end
aimDist(aimDist~=1)=0;
% distMatrix(aimOrder,aimOrder)=aimDist;
distMatrix=aimDist;
% imshow(distMatrix)
pos1=[];
pos2=[];
for i=1:size(distMatrix,1)
    for j=1:size(distMatrix,1)
        if distMatrix(i,j)==1
            if branchIndex(i)==branchIndex(j)
                pos1=[pos1;i,j];
            else
                pos2=[pos2;i,j];
            end
        end
    end
end
hold on
if ~isempty(pos1)
    plot(pos1(:,1),pos1(:,2),'Marker','.','LineStyle','none','Color',[1,0,0],'MarkerSize',6);
end
if ~isempty(pos2)
    plot(pos2(:,1),pos2(:,2),'Marker','.','LineStyle','none','Color',[0,1,0],'MarkerSize',6)
end
set(gca,'XLim',[0,num],'YLim',[0,num])
% axis xy
end
function [positionInfo,positionNum]=generatePos(leafNum,centroidInfo,treeSize)
leafInfo=centroidInfo(end-leafNum+1:end);
is2End=leafInfo==treeSize;
positionNum=numel(is2End(is2End==1));
positionInfo=zeros(leafNum,1);
positionInfo(is2End)=1:positionNum;
orderNum=1:leafNum;
orderNum=orderNum(is2End);
for i=1:positionNum
    if i==1
        if orderNum(i)~=1
            positionInfo(1:orderNum(i)-1)=1/orderNum(i)*(1:(orderNum(i)-1));
        end
    end
    if i==positionNum
        if orderNum(i)~=leafNum;
            positionInfo(orderNum(i)+1:end)=i+1/(leafNum-orderNum(i)+1)*(1:(leafNum-orderNum(i)));
        end
        continue
    end
    if orderNum(i)+1~=orderNum(i+1)
        positionInfo(orderNum(i)+1:orderNum(i+1)-1)=i+1/(orderNum(i+1)-orderNum(i))*(1:(orderNum(i+1)-orderNum(i)-1));
    end
end
if positionNum==0;
    positionInfo(1:leafNum)=1/(leafNum+1)*(1:leafNum);
end
end
function rightNode=morphCorrection(bioTree,nodeInfo)
    pictureSize=bioTree{1}.imageSize;
    if ~(size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In,2)==1 && size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out,2)==2)
        rightNode=0;
        return
    end
    rightNode=1;
end