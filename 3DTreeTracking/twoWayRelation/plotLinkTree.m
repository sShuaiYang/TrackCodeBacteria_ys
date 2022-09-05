function plotLinkTree(linkMatrix,centroidInfo,leafYCoo,color,type,cutNum)
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
color=[1,0,0];
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
            if nargin==5 && strcmp(type,'aging') 
                linkTwoPointsAging(centroidInfo(i,:),centroidInfo(j,:),color,nodeNum,i,j,linkMatrix(i,j));
            else
                linkTwoPoints(centroidInfo(i,:),centroidInfo(j,:),color,nodeNum,i,j,linkMatrix(i,j),cutNum);
            end
        end
    end
end
end
function linkTwoPoints(c1,c2,color,nodeNum,i,j,para,cutNum)
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
line(x,y,'Color',color,'LineWidth',1)


maker='.';
size=16;
if i==1
    plot(c1(1),c1(2),'Marker','^','MarkerSize',12,'Color',[1,0,0]);
else
    if j<=nodeNum
        plot(c1(1),c1(2),'Marker',maker,'MarkerSize',size,'Color',[205,155,29]/255);
        if c2(1)<=cutNum
            plot(c2(1),c2(2),'Marker',maker,'MarkerSize',size,'Color',[205,155,29]/255);
        end
    end
    if j>=nodeNum
        plot(c1(1),c1(2),'Marker',maker,'MarkerSize',size,'Color',[205,155,29]/255);
        if c2(1)<=cutNum
            plot(c2(1),c2(2),'Marker',maker,'MarkerSize',size,'Color',[85,107,47]/255);
        end
    end
end
end
function linkTwoPointsAging(c1,c2,color,nodeNum,i,j,para)
% c1(1)=c1(1)/1200;
% c2(1)=c2(1)/1200;
% support c1(x)<c2(x),link c1-c3-c2

% debug plus para(linkMatrix(i,j))
if para==-inf
    color=[0,0,0];
else
    colorAll=colormap(jet(10));
    color=colorAll(para,:);
end
    
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
