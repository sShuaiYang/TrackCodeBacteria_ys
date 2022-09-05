function getFigureMap(bioTree,bacteriaFrameInfo,cutNum)
sizeNeed=1200*10;
% sizeNeed=cutNum;
bioTree=bioTreeCut(bioTree,sizeNeed);
bacteriaFrameInfo=bacteriaFrameInfo(1:sizeNeed);
% aveLength=40;
finalInfo=bacteriaFrameInfo{end};
aveLength=mean(finalInfo.lengthInfo(finalInfo.lengthInfo<=80));
multi=3;
[bioTree,branchList,~,~]=myBiograph_new2(bioTree);
isCoreBranch=branchList(:,3);
coreBranchNum=numel(isCoreBranch(isCoreBranch==1));
branchNum=coreBranchNum;
colorBar=colormap(jet(coreBranchNum+1));
close all
figure1=figure;
treeSize=numel(bioTree);
if coreBranchNum<=20
    branchNum=coreBranchNum;
else
    branchNum=20;
end

% down map
num=1;
axes('Parent',figure1,'Units','inches','Position',[7 0.5 7 1.5],'xLim',[0,numel(bioTree)/1200],'box','on');
leafAll=[];
for i=1:branchNum
    if i<=coreBranchNum
        color=colorBar(i,:);
    else
        color=[0,0,1];
    end
    [linkMatrix,centroidInfo,leafList,leafNum]=generateOneBranchPhytree(bioTree,i,branchList(i,:),'normal');
    [positionInfo,positionNum]=generatePos(size(leafList,1),centroidInfo,numel(bioTree));
    leafYCoo=num+positionInfo;
    num=num+positionNum;
    leafAll=[leafAll;leafList];
    hold on;plotLinkTree(linkMatrix,centroidInfo,leafYCoo,color,treeSize);

%     [linkMatrix,centroidInfo,leafList,leafNum,allList,getLeaf,allLeaf,linkTwoPoint]=clusterLikeOneBranchPhytree(bioTree,i,branchList(i,:),branchList);
%     leafYCoo=num+(0:size(allLeaf,1)-1);
%     num=num+size(allLeaf,1);
%     leafAll=[leafAll;allLeaf];
%     hold on;clusterLikeTreePlot(linkMatrix,centroidInfo,leafYCoo,color,allLeaf,getLeaf,allList,leafNum,linkTwoPoint);
end
set(gca,'yLim',[1,num]);
view(gca,[90 -90]);

% left map
num=1;
axes('Parent',figure1,'Units','inches','Position',[4.646 2 2.354 7],'xLim',[0,numel(bioTree)/1200],'box','on','YTickLabel',{});
leafAll=[];
for i=1:branchNum
    if i<=coreBranchNum
        color=colorBar(i,:);
    else
        color=[0,0,1];
    end;
    [linkMatrix,centroidInfo,leafList,leafNum]=generateOneBranchPhytree(bioTree,i,branchList(i,:),'normal');
    [positionInfo,positionNum]=generatePos(size(leafList,1),centroidInfo,numel(bioTree));
    leafYCoo=num+positionInfo;
    num=num+positionNum;
    leafAll=[leafAll;leafList];
    hold on;plotLinkTree(linkMatrix,centroidInfo,leafYCoo,color,treeSize);

%     [linkMatrix,centroidInfo,leafList,leafNum,allList,getLeaf,allLeaf,linkTwoPoint]=clusterLikeOneBranchPhytree(bioTree,i,branchList(i,:),branchList);
%     leafYCoo=num+(0:size(allLeaf,1)-1);
%     num=num+size(allLeaf,1);
%     leafAll=[leafAll;allLeaf];
%     hold on;clusterLikeTreePlot(linkMatrix,centroidInfo,leafYCoo,color,allLeaf,getLeaf,allList,leafNum,linkTwoPoint);
end
set(gca,'yLim',[1,num]);

% middle map
axes('Parent',figure1,'Units','inches','Position',[7 2 7 7],'Color',[0,0,0],'box','on','XTickLabel',{},'YTickLabel',{});
hold on
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
for i=1:numel(orderNum)
    aimDist(i,i)=0;
end
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
% plot(pos1(:,1),pos1(:,2),'Marker','.','LineStyle','none','Color',[1,0,0],'MarkerSize',4);
% plot(pos2(:,1),pos2(:,2),'Marker','.','LineStyle','none','Color',[0,1,0],'MarkerSize',4);
if ~isempty(pos1)
    plot(pos1(:,1),pos1(:,2),'Marker','.','LineStyle','none','Color',[1,0,0],'MarkerSize',4);
end
if ~isempty(pos2)
    plot(pos2(:,1),pos2(:,2),'Marker','.','LineStyle','none','Color',[0,1,0],'MarkerSize',4);
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
function plotLinkTree(linkMatrix,centroidInfo,leafYCoo,color,treeSize)
% generally linkMatrix means the relationship for [N, L], N means all
% node, L means all leaf, centroidInfo here could be changed for the second
% column. The first column means the frame that the node begin\
% color=[1,0,1];
leafNum=numel(leafYCoo);
nodeNum=size(linkMatrix,1)-leafNum;
centroidInfo(nodeNum+1:end,2)=leafYCoo;
leafOrder=1:leafNum;
for i=1:nodeNum
    iLine=linkMatrix(nodeNum+1-i,nodeNum+1-i:end);
    centroid2=centroidInfo(:,2);
    centroid2=centroid2(nodeNum+1-i:end);
    leafPos=(centroid2(iLine==1));
    centroidInfo(nodeNum+1-i,2)=mean(leafPos(~isnan(leafPos)));
end
% figure;
for i=1:size(linkMatrix,1)
    for j=i+1:size(linkMatrix,1)
        if linkMatrix(i,j)==1
            linkTwoPoints(centroidInfo(i,:),centroidInfo(j,:),color,nodeNum,i,j,treeSize);
        end
    end
end
end
function linkTwoPoints(c1,c2,color,nodeNum,i,j,treeSize)
c1(1)=c1(1)/1200;
c2(1)=c2(1)/1200;
% support c1(x)<c2(x),link c1-c3-c2
if c1(1)>c2(1)
    c=c2;
    c2=c1;
    c1=c;
end
c3=[c1(1),c2(2)];
hold on;
line([c1(1),c3(1),c2(1)],[c1(2),c3(2),c2(2)],'Color',color,'LineWidth',1);
maker='.';
size=16;
if i==1
    plot(c1(1),c1(2),'Marker',maker,'MarkerSize',size,'Color',[1,0,0]);
else
    if j<=nodeNum
        plot(c1(1),c1(2),'Marker',maker,'MarkerSize',size,'Color',[205,155,29]/255);
        plot(c2(1),c2(2),'Marker',maker,'MarkerSize',size,'Color',[205,155,29]/255);
    end
    if j>=nodeNum
        plot(c1(1),c1(2),'Marker',maker,'MarkerSize',size,'Color',[205,155,29]/255);
        if c2(1)<treeSize
            plot(c2(1),c2(2),'Marker',maker,'MarkerSize',size,'Color',[85,107,47]/255);
        end
    end
end
end