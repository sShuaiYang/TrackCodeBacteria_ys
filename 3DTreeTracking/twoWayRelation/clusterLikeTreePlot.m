function clusterLikeTreePlot(linkMatrix,centroidInfo,allLeafCoo,color,allLeaf,getLeaf,allList,leafNum,linkTwoPoint)
% generally linkMatrix means the relationship for [N, L], N means all
% node, L means all leaf, centroidInfo here could be changed for the second
% column. The first column means the frame that the node begin\
% color=[1,0,1];

% 如果是getLeaf=0的时候，node形成的‘叶子’的位置需要重新定义，因为allLeafCoo定义的是最终所有叶子的位置
if getLeaf==0
    if ~isempty(allLeaf)
    minCoo=min(allLeafCoo);
    maxCoo=max(allLeafCoo);
    minLeafCoo=minCoo+double(maxCoo-minCoo)/6;
    maxLeafCoo=maxCoo-double(maxCoo-minCoo)/6;
    leafYCoo=linspace(minLeafCoo,maxLeafCoo,leafNum);
    else
        leafYCoo=allLeafCoo;
    end
else
    leafYCoo=allLeafCoo;
end
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
for i=1:size(linkMatrix,1)
    for j=i+1:size(linkMatrix,1)
        if linkMatrix(i,j)~=0
            linkTwoPoints(centroidInfo(i,:),centroidInfo(j,:),color,i,j,allList);
        end
    end
end
if getLeaf==0
    hold on
    if ~isempty(allLeaf)
        leafNum=size(allLeaf,1);
        xCen=ones(1,leafNum)*allLeaf(1,1);
        plot(xCen/1200,allLeafCoo,'lineStyle','none','Marker','.','MarkerSize',16,'Color',[85,107,47]/255);
        plot([xCen(1)/1200,centroidInfo(linkTwoPoint(1),1)/1200],[allLeafCoo(1),centroidInfo(linkTwoPoint(1),2)],'lineStyle','--','lineWidth',2);
        plot([xCen(1)/1200,centroidInfo(linkTwoPoint(2),1)/1200],[allLeafCoo(end),centroidInfo(linkTwoPoint(2),2)],'lineStyle','--','lineWidth',2);
    end
end
end
function linkTwoPoints(c1,c2,color,i,j,allList)
c1(1)=c1(1)/1200;
c2(1)=c2(1)/1200;
% support c1(x)<c2(x),link c1-c3-c2

% root [1,0,0]
% node [205,155,29]/255
% leaf [85,107,47]/255
    
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
switch allList(i,end)
    case -1
        plot(c1(1),c1(2),'Marker',maker,'MarkerSize',size,'Color',[1,0,0]);
    case 1
        plot(c1(1),c1(2),'Marker',maker,'MarkerSize',size,'Color',[205,155,29]/255);
    case 0
        plot(c1(1),c1(2),'Marker',maker,'MarkerSize',size,'Color',[85,107,47]/255);
end
switch allList(j,end)
    case -1
        plot(c2(1),c2(2),'Marker',marker,'MarkerSize',size,'Color',[1,0,0]);
    case 1
        plot(c2(1),c2(2),'Marker',maker,'MarkerSize',size,'Color',[205,155,29]/255);
    case 0
        plot(c2(1),c2(2),'Marker',maker,'MarkerSize',size,'Color',[85,107,47]/255);
end
end