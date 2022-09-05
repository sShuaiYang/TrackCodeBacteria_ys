function allData=detailsForBacteriaGrowthandMoving(allData,maxFrame,dirFile)
cd(dirFile)
allData=division2DivisionBacteria(allData);
   %this function has two parts, length statistic and velocity statistic,
   %    the first one uses a coarse smooth method for length and then get the mean value
   %  of every 30frames(1.5min),finallly get the growthRate in different growth stage
   %    the other one has a partly smooth function for x,y position and then
   %  get the absolute speed in different growth stage
division2DetachingBacteria(allData,maxFrame);
allData=attachingBacteria(allData);
end
function allData=division2DivisionBacteria(allData)
allData=lengthStatistic(allData);
% allData=velocityVsTime(allData);
end
function allData=velocityVsTime(allData)
velocityAll=[];
for iCell=1:size(allData,2)
    if allData{iCell}.isDivision==1 && allData{iCell}.is2Divi==1
        velocityInfo=[];
        positionInfo=allData{iCell}.bacteriaInfo(:,2:3);
        positionInfo(:,1)=sectionDenoise(positionInfo(:,1),[0.2,10],7);
        positionInfo(:,2)=sectionDenoise(positionInfo(:,2),[0.2,10],7);
        velocityInfo(:,1)=diff(positionInfo(:,1))';
        velocityInfo(:,2)=diff(positionInfo(:,2))';
        velocityInfo(:,3)=(velocityInfo(:,1).^2+velocityInfo(:,2).^2).^0.5;
        allData{iCell}.velocityInfo=velocityInfo;
        velocityAll=[velocityAll;cat(2,((1:numel(velocityInfo(:,3)))/numel(velocityInfo(:,3)))',velocityInfo(:,3))];        
    end  
end
h=createAtoBHistGram(linspace(0,1,81),logspace(-6,2,81),velocityAll(:,1),velocityAll(:,2),[0,1],'Time','velocityCen');
saveas(h,'division2Div_velocityCenDist.fig');
end
function allData=lengthStatistic(allData)
lengthAll=[];
for iCell=1:size(allData,2)
    if allData{iCell}.isDivision==1 && allData{iCell}.is2Divi==1
        lengthInfo=allData{iCell}.bacteriaInfo(:,1);
%         plot(lengthInfo)
        lengthAfterSmooth=mySmoothForLength(lengthInfo);
        if ~isempty(lengthAfterSmooth)
            allData{iCell}.LengthData(:,2)=lengthAfterSmooth;
            %         hold on;plot(allData{iCell}.LengthData(:,2),'r')
            allData{iCell}.maxLength=max(allData{iCell}.LengthData(:,2));
            allData{iCell}.LengthData(:,2)=allData{iCell}.LengthData(:,2)./allData{iCell}.maxLength;
            allData{iCell}.LengthData(:,1)=(1:size(allData{iCell}.LengthData(:,1),1))./size(allData{iCell}.LengthData(:,1),1);
            if allData{iCell}.LengthData(1,2)<0.4 || allData{iCell}.LengthData(1,2)>0.7
                allData{iCell}.LengthData=[];
            else
                lengthAll=[lengthAll;cat(2,allData{iCell}.LengthData(1:end-1,1),diff(allData{iCell}.LengthData(:,2)))];
            end
        end
        close all
    end
end
h=createAtoBHistGram(linspace(0,1,81),logspace(-7,-2,81),lengthAll(:,1),lengthAll(:,2),[0,1],'Time','GrowthRate');
saveas(h,'division2Div_growthRateDistribution.fig');
end

%% partly denoise
function newSignal=sectionDenoise(originalSignal,threshold,strength)
diffSignal=diff(originalSignal);
diffSignal=[0;diffSignal];
position=find(abs(diffSignal)>=threshold(1) & abs(diffSignal)<=threshold(2));
position=[1;position;numel(diffSignal)];
for i=1:numel(position)-2
    Data=originalSignal(position(i):position(i+1)-1,1);
    if numel(Data)>=3
        [~,sorh,keepapp] = ddencmp('den','wv',Data);
        thr=1/0.01;
        newSignal(position(i):position(i+1)-1,1)=wdencmp('gbl',Data,'db7',strength,thr,sorh,keepapp);
    else
        newSignal(position(i):position(i+1)-1,1)=Data;
    end
end
Data=originalSignal(position(end-1):position(end),1);
[~,sorh,keepapp] = ddencmp('den','wv',Data);
thr=1/0.01;
newSignal(position(end-1):position(end),1)=wdencmp('gbl',Data,'db7',strength,thr,sorh,keepapp);
end
%% A to B hist function
function h=createAtoBHistGram(binSpaceX,binSpaceY,funA,funB,logorNot,para1,para2)
[countMap1,~,~]=get2Dhist(funA,funB,binSpaceX,binSpaceY);
[xCount1,xBin1]=get1DhistforAll(funA,binSpaceX);
[yCount1,yBin1]=get1DhistforAll(funB,binSpaceY);
h=createHistfigure(log10(countMap1)', yBin1, xBin1, log10(countMap1)', xCount1, log10(yCount1),para1,para2,[min(binSpaceX),max(binSpaceX)],[min(binSpaceY),max(binSpaceY)],logorNot);
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
% count=count/sum(count);
Bin=binSpace(2:2:end);
end
function h=createHistfigure(ZData1, YData1, XData1, CData1, xCount1, yCount1,para1,para2,xRange,yRange,logorNot)
% Create figure
figure1 = figure();
colormap('jet');

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

if logorNot(1)==1
    set(axes1,'XScale','log');
    set(axes2,'XScale','log');
end
if logorNot(2)==1
    set(axes1,'YScale','log');
    set(axes3,'YScale','log');
end
h=gca;
end

%% useful function for length smooth
function smoothLength=mySmoothForLength(lengthInfo)
[~,sorh,keepapp] = ddencmp('den','wv',lengthInfo);
thr=0.5;
lengthInfo=wdencmp('gbl',lengthInfo,'db7',7,thr,sorh,keepapp);
if numel(lengthInfo)>30
    for i=1:numel(lengthInfo)-30
        smoothLength(i)=mean(lengthInfo(i:i+29));
    end
else
    smoothLength=[];
end
% for i=1:numel(lengthInfo)
% [~,sorh,keepapp] = ddencmp('den','wv',lengthInfo);
% thr=1.3;
% smoothLength=wdencmp('gbl',lengthInfo,'db7',7,thr,sorh,keepapp);
% diffForSmoLen=abs(diff(smoothLength));
% diffForSmoLen=[0;diffForSmoLen];
% isTooBig=diffForSmoLen>0.1;
% bigInfo=regionprops(isTooBig,'PixelIdxList');
% for i=1:size(bigInfo,1);
%     pixelIdx=bigInfo(i).PixelIdxList;
%     if pixelIdx(end)==numel(lengthInfo)
%         smoothLength(pixelIdx)=[];
%         continue
%     end
%     yy=spline([pixelIdx(1)-1,pixelIdx(end)+1],[smoothLength(pixelIdx(1)-1),smoothLength(pixelIdx(end)+1)],pixelIdx(1)-1:pixelIdx(end)+1);
%     smoothLength(pixelIdx)=yy(2:end-1);
% end
% j=2;
% i=1;
% while j<numel(smoothLength)
%     if smoothLength(j)>=smoothLength(i)
%         if j==i+1
%         end
%         if j~=i+1;
%             yy=spline([i,j],[smoothLength(i),smoothLength(j)],i:1:j);
%             smoothLength(i:j)=yy;
%         end
%     i=j;
%     j=j+1;
%     else
%         if smoothLength(j)<smoothLength(i)
%             j=j+1;
%         end
%     end
% end
end

%% statistic for division to detaching bacteria
function division2DetachingBacteria(allData,maxFrame)
divi2Info=[];
step=1000;
interval=[];
for i=1:size(allData,2)
    if allData{i}.isDivision==1
        if (allData{i}.is2Divi==1 || allData{i}.isDetaching==1) && allData{i}.frameNum(2)~=maxFrame+1
            divi2Info=[divi2Info;[allData{i}.frameNum(1),allData{i}.is2Divi]];
            if allData{i}.is2Divi==1
                interval=[interval;allData{i}.frameNum(end)-allData{i}.frameNum(1)];
            end
        end
    end
end
divi2InfoSta=[];
regionInfo=linspace(1,max(divi2Info(:,1)),fix(max(divi2Info(:,1))/step));
for i=1:numel(regionInfo)-1;
    pickInfo=divi2Info(divi2Info(:,1)>regionInfo(i) & divi2Info(:,1)<regionInfo(i+1),2);
    divi2InfoSta=[divi2InfoSta;(regionInfo(i)+regionInfo(i+1))/2,numel(pickInfo),numel(pickInfo(pickInfo==1))/numel(pickInfo)];
end
interval(interval<=400)=[];
h=createfigure1(divi2InfoSta(:,1),divi2InfoSta(:,3),divi2InfoSta(:,2),interval);
saveas(h,'div2DetachRateVsTime.fig')
end

function h=createfigure1(X1, Y1, Y2, C1)
%CREATEFIGURE(X1,Y1,Y2)
%  X1:  vector of x data
%  Y1:  vector of y data
%  Y2:  vector of y data

%  Auto-generated by MATLAB on 24-Nov-2012 16:54:40

% Create figure
figure1 = figure();

% Create axes
axes1 = axes('Parent',figure1,...
    'Position',[0.19375 0.493024453711505 0.277624316893035 0.427070636848726]);
box(axes1,'on');
hold(axes1,'all');

% Create plot
plot(X1,Y1,'Parent',axes1,'MarkerSize',15,'Marker','x');

% Create xlabel
xlabel('Time£¨frame£©','FontSize',14);

% Create ylabel
ylabel('division ratio','FontSize',14);

% Create axes
axes2 = axes('Parent',figure1,...
    'Position',[0.522916666666667 0.489706730404671 0.281770833333333 0.430095791959495],...
    'CLim',[1 2]);
box(axes2,'on');
hold(axes2,'all');

% Create plot
plot(X1,Y2,'Parent',axes2,'MarkerSize',15,'Marker','*');

% Create xlabel
xlabel('Time£¨frame£©','FontSize',14);

% Create ylabel
ylabel('case Num','FontSize',14);

% Create axes
axes3 = axes('Parent',figure1,...
    'Position',[0.194270833333333 0.0769351230425055 0.613020833333333 0.360425055928412]);
box(axes3,'on');
hist(C1)

% Create xlabel
xlabel('division interval','FontSize',14);

% Create ylabel
ylabel('case num','FontSize',14);
h=gca;
end

%% statistic for attaching bacteria that stay longer than a proper time
function allData=attachingBacteria(allData)
allData=attaching_velocityVsTime(allData,50);
allData=attaching_travellingDistance(allData);
end
function allData=attaching_velocityVsTime(allData,threShold)
velocityAll=[];
for iCell=1:size(allData,2)
    if allData{iCell}.isAttaching==1
        velocityInfo=[];
        if isempty(allData{iCell}.bacteriaInfo)
            continue
        end
        positionInfo=allData{iCell}.bacteriaInfo(:,2:3);
        if size(positionInfo,1)>=2
            positionInfo(:,1)=sectionDenoise(positionInfo(:,1),[0.2,10],7);
            positionInfo(:,2)=sectionDenoise(positionInfo(:,2),[0.2,10],7);
            velocityInfo(:,1)=diff(positionInfo(:,1))';
            velocityInfo(:,2)=diff(positionInfo(:,2))';
            velocityInfo(:,3)=(velocityInfo(:,1).^2+velocityInfo(:,2).^2).^0.5;
            allData{iCell}.velocityInfo=velocityInfo;
            allData{iCell}.positionInfo=positionInfo;
        end
        if size(positionInfo,1)>=threShold
             velocityAll=[velocityAll;cat(2,(1:numel(allData{iCell}.velocityInfo(:,3)))',allData{iCell}.velocityInfo(:,3))];
        end
    end  
end
h=createVelo2TimeHistGram(linspace(0,max(velocityAll(:,1)),81),logspace(-6,2,81),velocityAll(:,1),velocityAll(:,2),[0,1],'Time','velocityCen');
saveas(h,'attaching_VelocityVsTime');
end
function h=createVelo2TimeHistGram(binSpaceX,binSpaceY,funA,funB,logorNot,para1,para2)
[countMap1,~,~]=get2DhistVelo2Time(funA,funB,binSpaceX,binSpaceY);
countMap1=log10(countMap1);
countMap1(abs(countMap1)==inf)=min(min(countMap1(abs(countMap1)~=inf)));
[xCount1,xBin1]=get1DhistforAll(funA,binSpaceX);
[yCount1,yBin1]=get1DhistforAll(funB,binSpaceY);
h=createHistfigure(countMap1', yBin1, xBin1, countMap1', xCount1, log10(yCount1),para1,para2,[min(binSpaceX),max(binSpaceX)],[min(binSpaceY),max(binSpaceY)],logorNot);
end
function [countMap,xBin,yBin]=get2DhistVelo2Time(xP1,yP2,binSpaceX,binSpaceY)
countMap=zeros((size(binSpaceX,2)-1)/2,(size(binSpaceY,2)-1)/2);
for i=2:2:size(binSpaceX,2)
    dataCol=yP2(xP1>=binSpaceX(i-1)&xP1<binSpaceX(i+1));
    for j=2:2:size(binSpaceY,2)
        countMap(i/2,j/2)=size(yP2(dataCol>=binSpaceY(j-1)&dataCol<binSpaceY(j+1)),1);
    end
end
for i=1:size(countMap,1)
    countMap(i,:)=countMap(i,:)/sum(countMap(i,:));
end
xBin=binSpaceX(2:2:end);
yBin=binSpaceY(2:2:end);
end
function allData=attaching_travellingDistance(allData)
distanceInfo=[];
for iCell=1:size(allData,2)
    if allData{iCell}.isAttaching==1
        if isempty(allData{iCell}.bacteriaInfo)
            continue
        end
        positionInfo=allData{iCell}.bacteriaInfo(:,2:3);
        if size(positionInfo,1)>=2
            positionInfo=allData{iCell}.positionInfo;
            positionInfo(:,3)=((positionInfo(:,1)-positionInfo(1,1)).^2+(positionInfo(:,2)-positionInfo(1,2)).^2).^0.5;
            distanceInfo=[distanceInfo;cat(2,(1:numel(positionInfo(:,3)))',positionInfo(:,3))];
        end
    end
end
h=createVelo2TimeHistGram(linspace(0,max(distanceInfo(:,1)),81),linspace(0,max(distanceInfo(:,2)),81),distanceInfo(:,1),distanceInfo(:,2),[0,0],'Time','distance');
saveas(h,'attaching_tracevelDistance.fig');
end
