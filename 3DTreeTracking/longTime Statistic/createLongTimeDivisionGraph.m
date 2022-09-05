function createLongTimeDivisionGraph(bioTree,dirFile)
dirSave=strcat(dirFile,'\longTimeResult');
% cd(dirSave)
% divisionSeries(bioTree,dirSave);
getGrowthRate(bioTree);
end
%% find the maxLength[40,100] and minLength[20,55] of one bacteria and the length vs time statistic
function divisionInfo=divisionSeries(bioTree,dirSave)
divisionInfo.maxLength=[];
divisionInfo.minLength=[];
[bioTree,branchList,~,~]=myBiograph_new2(bioTree);
[bioTree,~,allList]=divisionFinder(bioTree,branchList);
allList.allNode(allList.allNode(:,4)==0,:)=[];
statisticDivisionTimeInfo(allList,bioTree);   % 这个图要经过某种变换
for iList=1:size(allList.allNode,1)
    divisionNode=allList.allNode(iList,:);
    if bioTree{divisionNode(1)}.node{divisionNode(2)}.In{1}.isNode==1
        nodeInfo=bioTree{divisionNode(1)}.node{divisionNode(2)}.In{1}.nodeInfo;
        divisionInfo.maxLength=[divisionInfo.maxLength;[divisionNode(1),bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.measurment{end}.MajorAxisLength]];
    else if bioTree{divisionNode(1)}.node{divisionNode(2)}.In{1}.isNode==0
            rootInfo=bioTree{divisionNode(1)}.node{divisionNode(2)}.In{1}.rootInfo;
            divisionInfo.maxLength=[divisionInfo.maxLength;[divisionNode(1),bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.measurment{end}.MajorAxisLength]];
        end
    end
    for iOut=1:size(bioTree{divisionNode(1)}.node{divisionNode(2)}.Out,2)
        divisionInfo.minLength=[divisionInfo.minLength;[divisionNode(1),bioTree{divisionNode(1)}.node{divisionNode(2)}.Out{iOut}.traceInfo.measurment{1}.MajorAxisLength]];
    end
end
h1=createDivisionLengthMinfigure(divisionInfo.minLength(:,1),divisionInfo.minLength(:,2));
saveas(h1,'minLength.fig');
h2=createDivisionLengthMaxfigure(divisionInfo.maxLength(:,1),divisionInfo.maxLength(:,2));
saveas(h2,'maxLength.fig');
binSpaceX=linspace(1,size(bioTree,2),61);
binSpaceY=linspace(20,100,61);
frameLengthInfo=getEachBacteriaLengthInFrame(bioTree);
[countMap1,~,~]=get2Dhist(frameLengthInfo(:,1),frameLengthInfo(:,2),binSpaceX,binSpaceY);
[xCount1,xBin1]=get1DhistTime(frameLengthInfo(:,1),binSpaceX);
[yCount1,yBin1]=get1DhistforAll(frameLengthInfo(:,2),binSpaceY);
h3=createHistfigure(countMap1', yBin1, xBin1, countMap1', xCount1, yCount1,'Time','Length',[min(binSpaceX),max(binSpaceX)],[min(binSpaceY),max(binSpaceY)]);
saveas(h3,'aveLength.fig');
end
function statisticDivisionTimeInfo(allList,bioTree)
% this function is used to get the division period interval and division
% point vs time
divisionTimePoint=allList.allNode(:,1);
divisionInterval=[];
for iList=1:size(allList.allNode,1)
    nodeInfo=allList.allNode(iList,:);
    for i=1:2
       if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{i}.is2Node==1;
           nextNode=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{i}.nodeInfo;
           if size(bioTree{nextNode(1)}.node{nextNode(2)}.In,2)==1 && size(bioTree{nextNode(1)}.node{nextNode(2)}.Out,2)==2
               divisionInterval=[divisionInterval;nextNode(1)-nodeInfo(1)];
           end
       end
    end
end 
[count,Bin]=get1DhistforAll(divisionInterval,linspace(2000/5,8000/5,21));
sumAll=0;
time=0:1:360000;
for i=1:numel(count)
    sumAll=sumAll+count(i)*cos(1/(Bin(i)*5)*time);
end
end
function h=createDivisionLengthMinfigure(X1, Y1)
% Create figure
figure1=figure;
scrsz=get(0,'ScreenSize');
set(gcf,'PaperPositionMode','auto');
set(gcf,'Position',[1 1 scrsz(3) scrsz(4)]);

% Create axes
axes1 = axes('Parent',figure1,...
    'Position',[0.164635416666667 0.211563731931669 0.279114583333333 0.492763561580123]);
% Uncomment the following line to preserve the Y-limits of the axes
ylim(axes1,[20 55]);
box(axes1,'on');
hold(axes1,'all');

% Create plot
plot(X1,Y1,'Marker','o','LineStyle','none');

% Create axes
axes2 = axes('Parent',figure1,...
    'Position',[0.485937500000001 0.208935611038108 0.268749999999999 0.495351511169514]);
box(axes2,'on');
hist(Y1,linspace(20,55,14));
xlim(axes2,[20,55]);
view(axes2,[90 -90]);

% Create textbox
annotation(figure1,'textbox',...
    [0.392666666666667 0.784210526315789 0.228166666666667 0.0447368421052632],...
    'String',{'Length After Division'},...
    'FontSize',20,...
    'FitBoxToText','off',...
    'LineStyle','none');
h=gca;
end
function h=createDivisionLengthMaxfigure(X1, Y1)
% Create figure
figure1=figure;
scrsz=get(0,'ScreenSize');
set(gcf,'PaperPositionMode','auto');
set(gcf,'Position',[1 1 scrsz(3) scrsz(4)]);

% Create axes
axes1 = axes('Parent',figure1,...
    'Position',[0.164635416666667 0.211563731931669 0.279114583333333 0.492763561580123]);
% Uncomment the following line to preserve the Y-limits of the axes
ylim(axes1,[40 100]);
box(axes1,'on');
hold(axes1,'all');

% Create plot
plot(X1,Y1,'Marker','o','LineStyle','none');

% Create axes
axes2 = axes('Parent',figure1,...
    'Position',[0.485937500000001 0.208935611038108 0.268749999999999 0.495351511169514]);
box(axes2,'on');
hist(Y1,linspace(40,100,12));
xlim(axes2,[40,100]);
view(axes2,[90 -90]);

% Create textbox
annotation(figure1,'textbox',...
    [0.392666666666667 0.784210526315789 0.228166666666667 0.0447368421052632],...
    'String',{'Length before Division'},...
    'FontSize',20,...
    'FitBoxToText','off',...
    'LineStyle','none');
h=gca;
end
function frameLengthInfoFinal=getEachBacteriaLengthInFrame(bioTree)
for iframe=1:size(bioTree,2)
    disp(iframe)
    frameLengthInfo{iframe}=[];
    for smallFrame=1:iframe
        if ~isempty(bioTree{smallFrame}.root)
            if ~isempty(bioTree{smallFrame}.root)
                for iRoot=1:size(bioTree{smallFrame}.root,2)
                    traceInfo=bioTree{smallFrame}.root{iRoot}.traceInfo.measurment;
                    if size(traceInfo,2)>=iframe+1-smallFrame
                        frameLengthInfo{iframe}=[frameLengthInfo{iframe};[iframe,traceInfo{iframe+1-smallFrame}(1).MajorAxisLength]];
                    end
                end
            end
        end
        if ~isempty(bioTree{smallFrame}.node)
            for iNode=1:size(bioTree{smallFrame}.node,2)
                for iOut=1:size(bioTree{smallFrame}.node{iNode}.Out,2);
                    traceInfo=bioTree{smallFrame}.node{iNode}.Out{iOut}.traceInfo.measurment;
                    if size(traceInfo,2)>=iframe+1-smallFrame
                        frameLengthInfo{iframe}=[frameLengthInfo{iframe};[iframe,traceInfo{iframe+1-smallFrame}(1).MajorAxisLength]];
                    end
                end
            end
        end
    end
end
frameLengthInfoFinal=[];
for iframe=1:size(bioTree,2)
    frameLengthInfoFinal=[frameLengthInfoFinal;frameLengthInfo{iframe}];
end
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
function [count,Bin]=get1DhistforAll(data,binSpace)
count=zeros(1,(size(binSpace,2)-1)/2);
for i=2:2:size(binSpace,2)
    count(i/2)=size(data(data>=binSpace(i-1)&data<binSpace(i+1)),1);
end
count=count/sum(count);
Bin=binSpace(2:2:end);
end
function h=createHistfigure(ZData1, YData1, XData1, CData1, xCount1, yCount1,para1,para2,xRange,yRange)
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
xlim(axes2,xRange);
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
h=gca;
end

%% growth rate distribution(each bacteria in each frame has its growth rate)
function [growthRateInfo,bioTree]=getGrowthRate(bioTree)
growthRateInfo=[];
aveGrowthRateVSgrowthRate=[];
for iframe=1:size(bioTree,2)
    for iRoot=1:size(bioTree{iframe}.root,2)
        lengthInfo=[];
        for iTrace=1:size(bioTree{iframe}.root{iRoot}.traceInfo.measurment,2)
            lengthInfo=[lengthInfo;bioTree{iframe}.root{iRoot}.traceInfo.measurment{iTrace}(1).MajorAxisLength];
        end
        [thr,sorh,keepapp] = ddencmp('den','wv',lengthInfo);
        newLength=wdencmp('gbl',lengthInfo,'db5',5,thr,sorh,keepapp);
        diffLength=diff(newLength);
        diffLength(diffLength<0)=NaN;
        growthRateInfo=[growthRateInfo;diffLength];
        if ~isempty(diffLength)
            p=polyfit(1:size(lengthInfo,1),lengthInfo',1);
            aveGrowthRate=p(1);
            aveGrowthRateInfo(1:size(diffLength,1),1)=aveGrowthRate;
            aveGrowthRateInfo(:,2)=diffLength;
            aveGrowthRateVSgrowthRate=[aveGrowthRateVSgrowthRate;aveGrowthRateInfo];
            aveGrowthRateInfo=[];
        end
    end
    for iNode=1:size(bioTree{iframe}.node,2)
        for iOut=1:size(bioTree{iframe}.node{iNode}.Out,2)
            lengthInfo=[];
            for iTrace=1:size(bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.measurment,2)
                lengthInfo=[lengthInfo;bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.measurment{iTrace}(1).MajorAxisLength];
            end
            [thr,sorh,keepapp] = ddencmp('den','wv',lengthInfo);
            newLength=wdencmp('gbl',lengthInfo,'db5',5,thr,sorh,keepapp);
            diffLength=diff(newLength);
            diffLength(diffLength<0)=NaN;
            growthRateInfo=[growthRateInfo;diffLength];
            if ~isempty(diffLength)
                p=polyfit(1:size(lengthInfo,1),lengthInfo',1);
                aveGrowthRate=p(1);
                aveGrowthRateInfo(1:size(diffLength,1),1)=aveGrowthRate;
                aveGrowthRateInfo(:,2)=diffLength;
                aveGrowthRateVSgrowthRate=[aveGrowthRateVSgrowthRate;aveGrowthRateInfo];
                aveGrowthRateInfo=[];
            end
        end
    end
end
growthRateInfo(isnan(growthRateInfo))=[];
clear growthRateInfo
aveGrowthRateVSgrowthRate(isnan(aveGrowthRateVSgrowthRate(:,2)),:)=[];
% binspace=logspace(-8,2,61);
% [count n]=get1Dhist(growthRateInfo,binspace);
% h1=createGrowthRateDistributionfigure(n,count);
% saveas(h1,'growthRate.fig');
binSpaceX=logspace(-2,-1,91);
binSpaceY=logspace(-8,0,91);
[countMap,~,~]=get2Dhist(aveGrowthRateVSgrowthRate(:,1),aveGrowthRateVSgrowthRate(:,2),binSpaceX,binSpaceY);
[xCount,xBin]=get1Dhist(aveGrowthRateVSgrowthRate(:,1),binSpaceX);
[yCount,yBin]=get1Dhist(aveGrowthRateVSgrowthRate(:,2),binSpaceY);
countMap=log10(countMap);
countMap(isnan(countMap))=0;
countMap(abs(countMap)==inf)=0;
h1=createHistfigure1(countMap',yBin,xBin, countMap', xCount, log10(yCount),[1,1],[min(binSpaceX),max(binSpaceX)],[min(binSpaceY),max(binSpaceY)]);
saveas(h1,'growthRateVsAveRate');
end
function [count,n]=get1DhistTime(velocityData,binSpace)
nBin=binSpace;
count=[];
n=[];
for i=2:2:size(nBin,2)
    count=[count,numel(velocityData(velocityData>=nBin(i-1)&velocityData<nBin(i+1)))/(nBin(i+1)-nBin(i-1))];
    n=[n,nBin(i)];
end
end
function [count,n]=get1Dhist(velocityData,binSpace)
nBin=binSpace;
count=[];
n=[];
for i=2:2:size(nBin,2)
    count=[count,numel(velocityData(velocityData>=nBin(i-1)&velocityData<nBin(i+1)))];
    n=[n,nBin(i)];
end
end
function h1=createGrowthRateDistributionfigure(X1,Y1)
% Create figure
figure1 = figure;
scrsz=get(0,'ScreenSize');
set(gcf,'PaperPositionMode','auto');
set(gcf,'Position',[1 1 scrsz(3) scrsz(4)]);

% Create axes
axes1 = axes('Parent',figure1,'YScale','log','YMinorTick','on',...
    'XScale','log',...
    'XMinorTick','on',...
    'Position',[0.128958333333333 0.11 0.775 0.815]);
box(axes1,'on');
hold(axes1,'all');

% Create loglog
loglog(X1,Y1,'Marker','o');

% Create title
title('growth rate distribution','FontSize',20);
axis equal
h1=gca;
end
function h=createHistfigure1(ZData1, YData1, XData1, CData1, xCount1, yCount1,logorNot,xRange,yRange)
%CREATEFIGURE(ZDATA1,YDATA1,XDATA1,CDATA1,XCOUNT1,YCOUNT1)
%  ZDATA1:  surface zdata
%  YDATA1:  surface ydata
%  XDATA1:  surface xdata
%  CDATA1:  surface cdata
%  XCOUNT1:  vector of y data
%  YCOUNT1:  vector of x data

%  Auto-generated by MATLAB on 18-Oct-2012 19:59:03

% Create figure
scrsz=get(0,'ScreenSize');
figure1 = figure('Renderer','zbuffer',...
    'PaperSize',[20.98404194812 29.67743169791],...
    'InvertHardcopy','off',...
    'Color',[1 1 1]);
set (gcf,'Position',[10 10 scrsz(3) scrsz(4)-80]);
set(gcf, 'PaperPositionMode', 'auto');
colormap('pink');

% Create axes
axes1 = axes('Parent',figure1,'XMinorTick','on',...
    'Position',[0.324036833770205 0.256068911511355 0.258191349934469 0.418950665622553]);
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
xlabel('aveGrowthRate','FontSize',14);

% Create ylabel
ylabel('growthRatePerFrame','FontSize',14);
% Create axes
axes2 = axes('Parent',figure1,'XMinorTick','on',...
    'Position',[0.325764389471385 0.705559906029757 0.255757972913936 0.151135473766641]);
% Uncomment the following line to preserve the X-limits of the axes
% xlim(axes2,[0.0001 0.1]);
box(axes2,'on');
hold(axes2,'all');

% Create semilogx
semilogx1=semilogx(XData1,xCount1,'Parent',axes2,'Marker','o','Color',[0 0 1],...
    'DisplayName','xCount vs xBin');
set(semilogx1(1),'DisplayName','空白');

% Create axes
axes3 = axes('Parent',figure1,...
    'Position',[0.599367218217563 0.258418167580266 0.0795107033639144 0.412685982772122]);

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
h=gcf;
end