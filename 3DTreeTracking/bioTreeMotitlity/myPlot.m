function  serialNum=myPlot(para1,islog1,xRange,para2,islog2,yRange,group1,varargin)
dirSave=uigetdir();
% dirSave='C:\Users\Administrator\Desktop\ϸ�����ݿ�';
if mod(nargin,2)==1  %%possible type: myPlot('speedH1','speedL1',[1,3,4,5],'time',[2,5],'standPercent',0)
    dirResultSave=strcat(dirSave,'\ͳ�Ʋ���\');
    otherInElement=varargin;
    [xdata,ydata,serialNum]=getXYdataTwoGroup(dirResultSave,group1,para1,para2,otherInElement);
    if islog1==1
        binSpaceX=logspace(xRange(1),xRange(2),91);
    elseif islog1==0
        binSpaceX=linspace(xRange(1),xRange(2),91);
    end
    if islog2==1
        binSpaceY=logspace(yRange(1),yRange(2),91);
    elseif islog2==0
        binSpaceY=linspace(yRange(1),yRange(2),91);
    end
    [countMap,~,~]=get2Dhist(xdata,ydata,binSpaceX,binSpaceY);
    [xCount,xBin]=get1Dhist(xdata,binSpaceX);
    [yCount,yBin]=get1Dhist(ydata,binSpaceY);
    createHistfigure(countMap', yBin, xBin, countMap', xCount, yCount,para1,para2,[islog1,islog2],[min(binSpaceX),max(binSpaceX)],[min(binSpaceY),max(binSpaceY)]);
    %     figure;h=gcf;
    %     axes('Parent',h,'Color',[0,0,0]);
    %     hold on;plot(xdata,ydata,'lineStyle','none','marker','none')
    %     myColor=colormap(jet(size(group1,2)));
    %     for i=1:size(xdata,1)
    %         colorNum=find(group1==serialNum(i,4));
    %         text(xdata(i),ydata(i),num2str(i),'FontSize',6,'color',myColor(colorNum,:));
    %     end
    %
    %     % Create xlabel
    %     xlabel(para1,'FontSize',14);
    %
    %     % Create ylabel
    %     ylabel(para2,'FontSize',14);
    %
    %     % Create colormap
    %     colormap(jet(size(group1,2)))
    %     for iGroup=1:size(group1,2)
    %         yMaskLable{iGroup}=num2str(group1(iGroup));
    %     end
    %     h=colorbar;
    %     currentY=get(h,'YTick');
    %     deltaY=(currentY(end)-currentY(1))/size(group1,2)/2;
    %     colorbar('YTickLabel',yMaskLable,'YTick',linSpace(currentY(1)+deltaY,currentY(end)-deltaY,size(group1,2)));
end

% if mod(nargin,2)==1  %%possible type: myPlot('speedH1','speedL1',[1,3,4,5],'time',[2,5],'standPercent',0)
%     dirResultSave=strcat(dirSave,'\ͳ�Ʋ���\');
%     otherInElement=varargin;
%     [xdata,ydata,serialNum]=getXYdataTwoGroup(dirResultSave,group1,para1,para2,otherInElement);
%     figure;h=gcf;
%     axes('Parent',h,'Color',[0,0,0]);
%     hold on;plot(xdata,ydata,'lineStyle','none','marker','none')
%     for i=1:size(xdata,1)
%         text(xdata(i),ydata(i),num2str(i),'FontSize',6,'color',[0,1,1]);
%     end
%     
%     % Create xlabel
%     xlabel(para1,'FontSize',14);
%     
%     % Create ylabel
%     ylabel(para2,'FontSize',14);
% 
% end

if mod(nargin,2)==0  %%possible type: myPlot('speedH1','speedL1',[1,3,4,5],[2,8],'time',[2,5],'standPercent',0)
    dirResultSave=strcat(dirSave,'\ͳ�Ʋ���\');
    otherInElement=varargin(2:end);
    [xdata1,ydata1,serialNum1]=getXYdataTwoGroup(dirResultSave,group1,para1,para2,otherInElement);
    if islog1==1
        binSpaceX=logspace(xRange(1),xRange(2),91);
    elseif islog1==0
        binSpaceX=linspace(xRange(1),xRange(2),91);
    end
    if islog2==1
        binSpaceY=logspace(yRange(1),yRange(2),91);
    elseif islog2==0
        binSpaceY=linspace(yRange(1),yRange(2),91);
    end
    [countMap1,~,~]=get2Dhist(xdata1,ydata1,binSpaceX,binSpaceY);
    [xCount1,xBin1]=get1Dhist(xdata1,binSpaceX);
    [yCount1,yBin1]=get1Dhist(ydata1,binSpaceY);
    [axes1,axes2,axes3]=createHistfigure(countMap1', yBin1, xBin1, countMap1', xCount1, yCount1,para1,para2,[islog1,islog2],[min(binSpaceX),max(binSpaceX)],[min(binSpaceY),max(binSpaceY)]);
    group2=varargin{1};
    [xdata2,ydata2,serialNum2]=getXYdataTwoGroup(dirResultSave,group2,para1,para2,otherInElement);
    [xCount2,xBin2]=get1Dhist(xdata2,binSpaceX);
    [yCount2,yBin2]=get1Dhist(ydata2,binSpaceY);
    createHistOverlapfigure(xdata2,ydata2, yBin2, xBin2, xCount2, yCount2,[axes1,axes2,axes3]);
    serialNum=serialNum2;   
%     % Create colormap
%     colormap(myColor)
%     for iGroup=1:2
%         yMaskLable{iGroup}=num2str(iGroup);
%     end
%     h=colorbar;
%     currentY=get(h,'YTick');
%     deltaY=(currentY(end)-currentY(1))/2/2;
%     colorbar('YTickLabel',yMaskLable,'YTick',linspace(currentY(1)+deltaY,currentY(end)-deltaY,2));
end
end
function dataLine=getData(para,sheet)
switch para
    case 'speedH1'
        dataLine=sheet(:,1);
    case 'speedH2'
        dataLine=sheet(:,2);
    case 'speedL1'
        dataLine=sheet(:,3);
    case 'speedL2'
        dataLine=sheet(:,4);
    case 'veloH1'
        dataLine=sheet(:,5);
    case 'veloH2'
        dataLine=sheet(:,6);
    case 'veloL1'
        dataLine=sheet(:,7);
    case 'veloL2'
        dataLine=sheet(:,8);
    case 'jump1'
        dataLine=sheet(:,9);
    case 'jump2'
        dataLine=sheet(:,10);
    case 'omiga'
        dataLine=sheet(:,11);
    case 'MSDslope1'
        dataLine=sheet(:,12);
    case 'MSDslope2'
        dataLine=sheet(:,13);
    case 'MSDslopeCentroid'
        dataLine=sheet(:,14);
    case 'aveVelo1'
        dataLine=sheet(:,15);
    case 'aveVelo2'
        dataLine=sheet(:,16);
    case 'aveVeloCentroid'
        dataLine=sheet(:,17);
    case 'growthRate'
        dataLine=sheet(:,18);
    case 'aveLen'
        dataLine=sheet(:,19);
    case 'standPercent'
        dataLine=sheet(:,20);
    case 'time'
        dataLine=sheet(:,21);
    case 'frameNum'
        dataLine=sheet(:,22);
end
end
function [xdata,ydata,serialNum]=getXYdataTwoGroup(dirResultSave,group,para1,para2,otherInElement)
groupData.paraSheet=[];
groupData.serialNum=[];
for iData=1:size(group,2)
    dirFile=strcat(dirResultSave,'parameterResult',num2str(group(iData)));
    load(dirFile)
    groupData.paraSheet=[groupData.paraSheet;parameterAll.paraSheet];
    parameterAll.serialNum(:,4)=group(iData);
    groupData.serialNum=[groupData.serialNum;parameterAll.serialNum];
end
xdata=getData(para1,groupData.paraSheet);
ydata=getData(para2,groupData.paraSheet);
inNum=size(otherInElement,2);
limitNum=inNum/2;
limitMatrix=true(size(groupData.paraSheet,1),1);
for iLimit=1:limitNum
    iLimitData=getData(otherInElement{iLimit*2-1},groupData.paraSheet);
    iLimitData=iLimitData>=otherInElement{iLimit*2}(1) & iLimitData<=otherInElement{iLimit*2}(end);
    limitMatrix=iLimitData & limitMatrix;
end
xdata=xdata(limitMatrix==1);
ydata=ydata(limitMatrix==1);
serialNum=groupData.serialNum;
serialNum(limitMatrix==0,:)=[];
end
function [countMap,xBin,yBin]=get2Dhist(xP1,yP2,binSpaceX,binSpaceY)
countMap=zeros((size(binSpaceX,2)-1)/2,(size(binSpaceY,2)-1)/2);
for i=2:2:size(binSpaceX,2)
    dataCol=yP2(xP1>=binSpaceX(i-1)&xP1<binSpaceX(i+1));
    for j=2:2:size(binSpaceY,2)
        countMap(i/2,j/2)=size(yP2(dataCol>=binSpaceY(j-1)&dataCol<binSpaceY(j+1)),1);
    end
end
xBin=binSpaceX(2:2:end);
yBin=binSpaceY(2:2:end);
end
function [count,Bin]=get1Dhist(data,binSpace)
count=zeros(1,(size(binSpace,2)-1)/2);
for i=2:2:size(binSpace,2)
    count(i/2)=size(data(data>=binSpace(i-1)&data<binSpace(i+1)),1);
end
count=count/sum(count);
Bin=binSpace(2:2:end);
end
function [axes1,axes2,axes3]=createHistfigure(ZData1, YData1, XData1, CData1, xCount1, yCount1,para1,para2,logorNot,xRange,yRange)
%CREATEFIGURE(ZDATA1,YDATA1,XDATA1,CDATA1,XCOUNT1,YCOUNT1)
%  ZDATA1:  surface zdata
%  YDATA1:  surface ydata
%  XDATA1:  surface xdata
%  CDATA1:  surface cdata
%  XCOUNT1:  vector of y data
%  YCOUNT1:  vector of x data

%  Auto-generated by MATLAB on 18-Oct-2012 19:59:03

% Create figure
figure1 = figure;
colormap('pink');

% Create axes
axes1 = axes('Parent',figure1,'XMinorTick','on',...
    'Position',[0.278724333770205 0.256068911511355 0.258191349934469 0.418950665622553]);
% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,xRange);
% Uncomment the following line to preserve the Y-limits of the axes
ylim(axes1,yRange);
hold(axes1,'all');

% Create surface
surface('Parent',axes1,'ZData',ZData1,'YData',YData1,'XData',XData1,...
    'LineStyle','none',...
    'CData',CData1);
%     % Create xlabel
xlabel(para1,'FontSize',14);

% Create ylabel
ylabel(para2,'FontSize',14);
% Create axes
axes2 = axes('Parent',figure1,'XMinorTick','on',...
    'Position',[0.279410222804718 0.705559906029757 0.255757972913936 0.151135473766641]);
% Uncomment the following line to preserve the X-limits of the axes
% xlim(axes2,[0.0001 0.1]);
box(axes2,'on');
hold(axes2,'all');

% Create semilogx
semilogx(XData1,xCount1,'Parent',axes2,'Marker','o','Color',[0 0 1],...
    'DisplayName','xCount vs xBin');

% Create axes
axes3 = axes('Parent',figure1,...
    'Position',[0.556138051550896 0.258418167580266 0.0795107033639144 0.412685982772122]);
% Uncomment the following line to preserve the X-limits of the axes
% xlim(axes3,[0 150]);
box(axes3,'on');
hold(axes3,'all');

% Create plot
plot(yCount1,YData1,'Parent',axes3,'Marker','o','Color',[0 0 1],...
    'DisplayName','yBin vs yCount');
if logorNot(1)==1
    set(axes1,'XScale','log');
    set(axes2,'XScale','log');
end
if logorNot(2)==1
    set(axes1,'YScale','log');
    set(axes3,'YScale','log');
end
end
function createHistOverlapfigure(xdata,ydata,YData1,XData1,xCount1,yCount1,axesInfo)
%CREATEFIGURE(ZDATA1,YDATA1,XDATA1,CDATA1,XCOUNT1,YCOUNT1)
%  ZDATA1:  surface zdata
%  YDATA1:  surface ydata
%  XDATA1:  surface xdata
%  CDATA1:  surface cdata
%  XCOUNT1:  vector of y data
%  YCOUNT1:  vector of x data

%  Auto-generated by MATLAB on 18-Oct-2012 19:59:03

% Create figure
plot(xdata,ydata,'o','Parent',axesInfo(1),'EraseMode','none','MarkerFaceColor',[0 1 1],'Marker','o','MarkerSize',5,...
    'LineStyle','none','Color',[0 1 1]);
% Create semilogx
semilogx(XData1,xCount1,'Parent',axesInfo(2),'Marker','o','Color',[1 0 0],...
    'DisplayName','xCount vs xBin');
% Create plot
plot(yCount1,YData1,'Parent',axesInfo(3),'Marker','o','Color',[1 0 0],...
    'DisplayName','yBin vs yCount');
end