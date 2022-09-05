function getFigureMap(bioTree,bacteriaFrameInfo)
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

% down map
ageAll=[];
num=1;
axes('Parent',figure1,'Units','inches','Position',[7 0.5 7 1.5],'xLim',[0,numel(bioTree)]/1200,'box','on');
leafAll=[];
for i=1:branchNum
    if i<=coreBranchNum
        color=colorBar(i,:);
    else
        color=[0,0,1];
    end
    [linkMatrix,centroidInfo,leafList,leafNum,ageAl]=generateOneBranchPhytree(bioTree,i,branchList(i,:),'aging');
    ageAll=[ageAll;ageAl];
    [positionInfo,positionNum]=generatePos(leafNum,centroidInfo,numel(bioTree));
    leafYCoo=num+positionInfo;
    num=num+positionNum;
    leafAll=[leafAll;leafList];
    hold on;plotLinkTree(linkMatrix,centroidInfo,leafYCoo,color,'aging');
end
ageAll=ageAll(leafAll(:,1)==numel(bacteriaFrameInfo));
set(gca,'yLim',[1,num]);
view(gca,[90 -90]);

% leaf map
num=1;
axes('Parent',figure1,'Units','inches','Position',[5.5 2 1.5 7],'xLim',[0,numel(bioTree)]/1200,'box','on');
leafAll=[];
for i=1:branchNum
    if i<=coreBranchNum
        color=colorBar(i,:);
    else
        color=[0,0,1];
    end;
    [linkMatrix,centroidInfo,leafList,leafNum]=generateOneBranchPhytree(bioTree,i,branchList(i,:),'aging');
    [positionInfo,positionNum]=generatePos(leafNum,centroidInfo,numel(bioTree));
    leafYCoo=num+positionInfo;
    num=num+positionNum;
    %     leafYCoo=(num:num+leafNum-1);
    %     num=num+leafNum;
    leafAll=[leafAll;leafList];
    hold on;plotLinkTree(linkMatrix,centroidInfo,leafYCoo,color,'aging');
end
set(gca,'yLim',[1,num]);

% middle map
axes('Parent',figure1,'Units','inches','Position',[7 2 7 7],'Color',[0,0,0],'box','on');
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
pos11=[];
pos12=[];
pos13=[];
pos2=[];
for i=1:size(distMatrix,1)
    for j=1:size(distMatrix,1)
        if distMatrix(i,j)==1
            if branchIndex(i)==branchIndex(j)
                % pos1=[pos1;i,j];
                if ageAll(i)==1 && ageAll(j)==1
                    pos11=[pos11;i,j];
                end
                if ageAll(i)~=1 && ageAll(j)==1 || ageAll(j)~=1 && ageAll(i)==1
                    pos12=[pos12;i,j];
                end
                if ageAll(i)==0 && ageAll(j)==0 || ageAll(j)==0 && ageAll(i)==0
                    pos13=[pos13;i,j];
                end
            else
                pos2=[pos2;i,j];
            end
        end
    end
end
hold on
plot(pos11(:,1),pos11(:,2),'Marker','.','LineStyle','none','Color',[1,0,0],'MarkerSize',4);
plot(pos12(:,1),pos12(:,2),'Marker','.','LineStyle','none','Color',[1,1,0],'MarkerSize',4);
plot(pos13(:,1),pos13(:,2),'Marker','.','LineStyle','none','Color',[0,1,0],'MarkerSize',4);
plot(pos2(:,1),pos2(:,2),'Marker','.','LineStyle','none','Color',[0,0,1],'MarkerSize',4);
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