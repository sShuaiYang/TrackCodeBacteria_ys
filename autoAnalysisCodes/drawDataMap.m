function [xData,yData]=drawDataMap(allData,dimNum,varargin)
%DRAWDATAMAP Summary of this function goes here
%   Detailed explanation goes here
% e.g drawDataMap(tree,'tree',1,'CyOFP')
if strcmp(dimNum,'tree')
    linkMatrix=zeros(size(allData(varargin{1}).linkMatrix));
    linkRow=allData(varargin{1}).linkRow;
    linkColumn=allData(varargin{1}).linkColumn;
    for i=1:numel(linkRow)
        switch varargin{2}
            case 'CyOFP'
                linkMatrix(linkRow(i),linkColumn(i))=mean(allData(varargin{1}).linkInfo(i).CyOFP(allData(varargin{1}).linkInfo(i).CyOFP~=0));
            case 'GFP'
                linkMatrix(linkRow(i),linkColumn(i))=mean(allData(varargin{1}).linkInfo(i).GFP(allData(varargin{1}).linkInfo(i).GFP~=0));
            case 'mScalet'
                linkMatrix(linkRow(i),linkColumn(i))=mean(allData(varargin{1}).linkInfo(i).mScalet(allData(varargin{1}).linkInfo(i).mScalet~=0));
            case 'RFP'
                linkMatrix(linkRow(i),linkColumn(i))=mean(allData(varargin{1}).linkInfo(i).RFP(allData(varargin{1}).linkInfo(i).RFP~=0));
            case 'Age'
                linkMatrix(linkRow(i),linkColumn(i))=mean(allData(varargin{1}).linkInfo(i).linkAge);
            case 'simpleTree'
                linkMatrix(linkRow(i),linkColumn(i))=20;
        end
    end
    plotLinkTree(linkMatrix,allData(varargin{1}).timer,1:allData(varargin{1}).leafNum,allData(varargin{1}).cutNum);
end

% e.g drawDataMap(tree,'generation',1,500,'CyOFP')  tree 1 , time =500min
if strcmp(dimNum,'generation')
    smallTree=allData(varargin{1});
    linkInfo=smallTree.linkInfo;
    linkGeneration{2}=[];
    for iLink=2:numel(linkInfo)
        targetInfo=smallTree.linkColumn(smallTree.linkRow==smallTree.linkRow(iLink));
        linkTagetNum=numel(targetInfo);
        index=find(targetInfo==smallTree.linkColumn(iLink));
        linkGeneration{smallTree.linkColumn(iLink)}=[linkGeneration{smallTree.linkRow(iLink)};index];
        linkInfo(iLink).linkGeneration=linkGeneration{smallTree.linkColumn(iLink)};
    end
    bacInfo=getLeafInfo(linkInfo,varargin{2},varargin{3});
    result=calculateCorrelation(bacInfo);
    yData=abs(result);
    xData=0:numel(yData)-1;
end
end
function iLine=getMaskResultLine(dataAll,inputVari)
switch inputVari
    case 'time'
        iLine=dataAll(:,1);
    case 'maxLen'
        iLine=dataAll(:,2);
    case 'minLen'
        iLine=dataAll(:,3);
    case 'CyOFP'
        iLine=dataAll(:,4);
    case 'GFP'
        iLine=dataAll(:,5);
    case 'mScalet'
        iLine=dataAll(:,6);
    case 'RFP'
        iLine=dataAll(:,7);
end
end
function iLine=getTrackingResultLine(dataAll,inputVari)
switch inputVari
    case 'time'
        iLine=dataAll(:,1);
    case 'maxLen'
        iLine=dataAll(:,2);
    case 'minLen'
        iLine=dataAll(:,3);
    case 'CyOFP'
        iLine=dataAll(:,4);
    case 'GFP'
        iLine=dataAll(:,5);
    case 'mScalet'
        iLine=dataAll(:,6);
    case 'RFP'
        iLine=dataAll(:,7);
    case 'angle'
        iLine=dataAll(:,8);
    case 'area'
        iLine=dataAll(:,9);
    case 'growthRate'
        iLine=dataAll(:,10);
    case 'fpProduce'
        iLine=dataAll(:,11);
end
end
function plotLinkTree(linkMatrix,centroidInfo,leafYCoo,cutNum)
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
    leafPos=(centroid2(iLine~=0));
    centroidInfo(nodeNum+1-i,2)=mean(leafPos(~isnan(leafPos)));
end
% figure;
colorAll=colormap(jet(124));
colorAll=colorAll(25:124,:);
close all
maxLinkNum=max(max(linkMatrix(linkMatrix~=0)));
minLinkNum=min(min(linkMatrix(linkMatrix~=0)));
colorRange=input(['please tap the color range ',num2str(minLinkNum) ,' to ',num2str(maxLinkNum),':']);
minLinkNum=colorRange(1);
maxLinkNum=colorRange(2);

% new add map size
figure1=figure;
axes1 = axes('Parent',figure1,...
    'Position',[0.13 0.11 0.288229166666667 0.74854189336235],'FontSize',20);
% ylim(axes1,[0 118])
% xlim(axes1,[0 12.6])
view(axes1,[90 -90]);
box(axes1,'on');
hold all
% Create xlabel
xlabel('time(h)','VerticalAlignment','bottom','Rotation',90,...
    'HorizontalAlignment','center','FontSize',20);
scrsz=get(0,'ScreenSize');
set (gcf,'Position',[10 10 scrsz(3) scrsz(4)-80]);
set(gcf, 'PaperPositionMode', 'auto');

for i=1:size(linkMatrix,1)
    for j=i+1:size(linkMatrix,1)
        if linkMatrix(i,j)~=0
            if maxLinkNum==minLinkNum
                colorIndex=1;
            else
                colorIndex=round(99*(linkMatrix(i,j)-minLinkNum)/(maxLinkNum-minLinkNum))+1;
            end
            if colorIndex<=1
                colorIndex=1;
            end
            if colorIndex>=100;
                colorIndex=100;
            end
            linkTwoPointsAging(centroidInfo(i,:),centroidInfo(j,:),nodeNum,i,j,cutNum,colorAll(colorIndex,:));
        end
    end
end
end
function linkTwoPointsAging(c1,c2,nodeNum,i,j,cutNum,colorI)
c1(1)=c1(1)/60;
c2(1)=c2(1)/60;
cutNum=cutNum/60;
if c1(1)>c2(1)
    c=c2;
    c2=c1;
    c1=c;
end
hold on;
if c1(1)>cutNum
    return
end
% link line
% c3=[c1(1),c2(2)];
% line([c1(1),c3(1),c2(1)],[c1(2),c3(2),c2(2)],'Color',color,'LineWidth',1);

% link with ellipse
c3=[c2(1),c1(2)];
axisA=abs(c2(1)-c1(1));
axisB=abs(c2(2)-c1(2));
if c2(1)>cutNum
    x=c1(1):0.01:cutNum;
else
    x=c1(1):0.01:c2(1);
end
if c2(2)<c1(2)
    y=c3(2)-axisB*(1-((x-c3(1))/axisA).^2).^(1/2);
else
    y=c3(2)+axisB*(1-((x-c3(1))/axisA).^2).^(1/2);
end
line(x,y,'Color',colorI,'LineWidth',1)

maker='.';
size=16;
if i==1
    plot(c1(1),c1(2),'Marker','^','MarkerSize',12,'Color',[1,0,0]);
else
    if j<=nodeNum
        plot(c1(1),c1(2),'Marker',maker,'MarkerSize',size,'Color',[205,155,29]/255);
        plot(c2(1),c2(2),'Marker',maker,'MarkerSize',size,'Color',[205,155,29]/255);
    end
    if j>=nodeNum
        plot(c1(1),c1(2),'Marker',maker,'MarkerSize',size,'Color',[205,155,29]/255);
        if c2(1)<360
            plot(c2(1),c2(2),'Marker',maker,'MarkerSize',size,'Color',[85,107,47]/255);
        end
    end
end
end
function bacInfo=getLeafInfo(linkInfo,givenTime,fluoProtein)
n=0;
for i=2:numel(linkInfo);
    linkTimer=linkInfo(i).linkTimer;
    diffTime=abs(linkTimer-givenTime);
    indexDiff=find(diffTime<=2);
    if isempty(indexDiff)
        continue
    end
    switch fluoProtein
        case 'CyOFP'
            result=linkInfo(i).CyOFP(indexDiff(2));
        case 'GFP'
            result=linkInfo(i).GFP(indexDiff(2));
        case 'mScalet'
            result=linkInfo(i).mScalet(indexDiff(2));
        case 'RFP'
            result=linkInfo(i).RFP(indexDiff(2));
    end
    if result~=0;
        n=n+1;
        bacInfo{n}.generateMark=linkInfo(i).linkGeneration;
        bacInfo{n}.meanIntensity=result;
    end
end
end
function result=calculateCorrelation(bacInfo)
allTable=[];
for i=1:numel(bacInfo)
    bacI.generateMark=bacInfo{i}.generateMark;
    bacI.meanIntensity=bacInfo{i}.meanIntensity;
    for j=1:numel(bacInfo)
        if j<=i
            bacJ.generateMark=bacInfo{j}.generateMark;
            bacJ.meanIntensity=bacInfo{j}.meanIntensity;
            dis=calDis(bacI.generateMark,bacJ.generateMark);
            allTable=[allTable;dis,bacI.meanIntensity,bacJ.meanIntensity];
        end
    end
end
for i=0:max(allTable(:,1))
    corrResult=corrcoef(allTable(allTable(:,1)==i,2),allTable(allTable(:,1)==i,3));
    result(i+1,1)=corrResult(1,2);
end
end
function dis=calDis(mark1,mark2)
% 分裂后两个子代之间的距离为1
sameGeneration=0;
for i=1:min(numel(mark1),numel(mark2))
    if mark1(i)==mark2(i)
        sameGeneration=i;
    else
        break
    end
end
dis=numel(mark1)-sameGeneration+numel(mark2)-sameGeneration-1;
if dis==-1
    dis=0;
end
end
