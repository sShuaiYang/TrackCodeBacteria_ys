function serialNum=myPlot(para1,para2,group1,varargin)
dirSave=uigetdir();
if mod(nargin,2)==1  %%possible type: myPlot('speedH1','speedL1',[1,3,4,5],'time',[2,5],'isStand',0)
    dirResultSave=strcat(dirSave,'\统计参数\');
    otherInElement=varargin;
    [xdata,ydata,serialNum]=getXYdataTwoGroup(dirResultSave,group1,para1,para2,otherInElement);
    figure;h=gcf;
    axes('Parent',h,'Color',[0,0,0]);
    hold on;plot(xdata,ydata,'lineStyle','none','marker','none')
    myColor=colormap(jet(size(group1,2)));
    for i=1:size(xdata,1)
        colorNum=find(group1==serialNum(i,4));
        text(xdata(i),ydata(i),num2str(i),'FontSize',6,'color',myColor(colorNum,:));
    end
    
    % Create xlabel
    xlabel(para1,'FontSize',14);
    
    % Create ylabel
    ylabel(para2,'FontSize',14);
    
    % Create colormap
    colormap(jet(size(group1,2)))
    for iGroup=1:size(group1,2)
        yMaskLable{iGroup}=num2str(group1(iGroup));
    end
    h=colorbar;
    currentY=get(h,'YTick');
    deltaY=(currentY(end)-currentY(1))/size(group1,2)/2;
    colorbar('YTickLabel',yMaskLable,'YTick',linSpace(currentY(1)+deltaY,currentY(end)-deltaY,size(group1,2)));
end

if mod(nargin,2)==0  %%possible type: myPlot('speedH1','speedL1',[1,3,4,5],[2,8],'time',[2,5],'isStand',0)
    dirResultSave=strcat(dirSave,'\统计参数\');
    otherInElement=varargin(2:end);
    [xdata1,ydata1,serialNum1]=getXYdataTwoGroup(dirResultSave,group1,para1,para2,otherInElement);
    group2=varargin{1};
    [xdata2,ydata2,serialNum2]=getXYdataTwoGroup(dirResultSave,group2,para1,para2,otherInElement);
    serialNum=[serialNum1;serialNum2];
    xdata=[xdata1;xdata2];
    ydata=[ydata1;ydata2];
    figure;h=gcf;
    axes('Parent',h,'Color',[0,0,0]);
    hold on;plot(xdata,ydata,'lineStyle','none','marker','none')
    myColor(1,:)=[0,1,1];
    myColor(2,:)=[1,1,0];
    for i=1:size(xdata,1)
        if i<=size(xdata1,1)
            text(xdata(i),ydata(i),num2str(i),'FontSize',6,'color',myColor(1,:));
        else
            text(xdata(i),ydata(i),num2str(i),'FontSize',6,'color',myColor(2,:));
        end
    end
    
    % Create xlabel
    xlabel(para1,'FontSize',14);
    
    % Create ylabel
    ylabel(para2,'FontSize',14);
    
    % Create colormap
    colormap(myColor)
    for iGroup=1:2
        yMaskLable{iGroup}=num2str(iGroup);
    end
    h=colorbar;
    currentY=get(h,'YTick');
    deltaY=(currentY(end)-currentY(1))/2/2;
    colorbar('YTickLabel',yMaskLable,'YTick',linspace(currentY(1)+deltaY,currentY(end)-deltaY,2));
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
    case 'aveVelo1'
        dataLine=sheet(:,14);
    case 'aveVelo2'
        dataLine=sheet(:,15);
    case 'growthRate'
        dataLine=sheet(:,16);
    case 'aveLen'
        dataLine=sheet(:,17);
    case 'standPercent'
        dataLine=sheet(:,18);
    case 'time'
        dataLine=sheet(:,19);
    case 'frameNum'
        dataLine=sheet(:,20);
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