function plotResultsAndSave(saveGraphic,dirSave,timeResult,histResult,clusterDataList,iData)
% in this function, we create six kinds of figures to show the motility
% characteristic of bacterias
saveFile1=strcat(dirSave,'\','histResult');
saveFile2=strcat(dirSave,'\','p1hist');
saveFile3=strcat(dirSave,'\','p2hist');
saveFile5=strcat(dirSave,'\','timeResult');
mkdir(saveGraphic);
saveGraphicFile=saveGraphic;
save(saveFile1,'histResult');
save(saveFile5,'timeResult');

%figure1 has two parts. First part shows us that how long a bacteria can
%move with velocity Vi,Second part tells us that if a bacteria moves with
%it the velocity of its head is V1,how long it moves while the tail of the
%bacteria is V2
h1=createP1VelocityHistgramfigure(histResult.p1p2MapP1sum, histResult.p1(:,1), histResult.p1p2MapP1sum, histResult.p1(:,2),clusterDataList,iData);
saveas(h1,saveFile2,'fig');
saveas(h1,strcat(saveGraphicFile,'\','p1hist'),'pdf');

h2=createP2VelocityHistgramfigure(histResult.p1p2MapP2sum, histResult.p2(:,1), histResult.p1p2MapP2sum, histResult.p2(:,2),clusterDataList,iData);
saveas(h2,saveFile3,'fig');
saveas(h2,strcat(saveGraphicFile,'\','p2hist'),'pdf');

%Orientation means the angle between V1 and the orientation of the major
%axis.Figure2 shows us how long a bacteria can move if its orientation is
%sita.
for iThreshold=1:size(histResult.p1VelocityOrientation,2)
    threshold=histResult.p1VelocityOrientation{iThreshold}.threshold;
    if ~isempty(threshold) && (threshold>0.1&&threshold<1)
     saveFile4=strcat(dirSave,'\','p1OrientationHist',num2str(iThreshold));
    dataList=histResult.p1VelocityOrientation{iThreshold}.histOrientation(:,3:end);
    theta=histResult.p1VelocityOrientation{iThreshold}.histOrientation(:,1);
     h3=myPolarPlot(theta,dataList,threshold);
     saveas(h3,saveFile4,'fig');
     saveas(h3,strcat(saveGraphicFile,'\','p1OrientationHist',num2str(iThreshold)),'pdf')
    end
end

for iThreshold=1:size(histResult.p2VelocityOrientation,2)
    threshold=histResult.p2VelocityOrientation{iThreshold}.threshold;
    if ~isempty(threshold) && (threshold>0.1&&threshold<1)
     saveFile4=strcat(dirSave,'\','p2OrientationHist',num2str(iThreshold));
    dataList=histResult.p2VelocityOrientation{iThreshold}.histOrientation(:,3:end);
    theta=histResult.p2VelocityOrientation{iThreshold}.histOrientation(:,1);
     h3=myPolarPlot(theta,dataList,threshold);
     saveas(h3,saveFile4,'fig');
     saveas(h3,strcat(saveGraphicFile,'\','p2OrientationHist',num2str(iThreshold)),'pdf')
    end
end


% % % % % Here Lei Ni change the orientation figures to two threshold specified figures 
% % % % allThreshold=[0.15,0.6];
% % % % for iThreshold=1:2
% % % % threshold=allThreshold(iThreshold);
% % % % if histResult.p1VelocityOrientation{end}.threshold>allThreshold(iThreshold)
% % % % saveFile4=strcat(dirSave,'\','p1OrientationHist',num2str(iThreshold));
% % % % dataList=histResult.p1VelocityOrientation{iThreshold}.histOrientation(:,3:end);
% % % % theta=histResult.p1VelocityOrientation{iThreshold}.histOrientation(:,1);
% % % % h3=myPolarPlot(theta,dataList,threshold);
% % % % saveas(h3,saveFile4,'fig');
% % % % saveas(h3,strcat(saveGraphicFile,'\','p1OrientationHist',num2str(iThreshold)),'pdf')
% % % % end
% % % % end
% % % % 
% % % % for iThreshold=1:2
% % % % threshold=allThreshold(iThreshold);
% % % % if histResult.p2VelocityOrientation{end}.threshold>allThreshold(iThreshold)
% % % % saveFile4=strcat(dirSave,'\','p2OrientationHist',num2str(iThreshold));
% % % % dataList=histResult.p2VelocityOrientation{iThreshold}.histOrientation(:,3:end);
% % % % theta=histResult.p2VelocityOrientation{iThreshold}.histOrientation(:,1);
% % % % h3=myPolarPlot(theta,dataList,threshold);
% % % % saveas(h3,saveFile4,'fig');
% % % % saveas(h3,strcat(saveGraphicFile,'\','p1OrientationHist',num2str(iThreshold)),'pdf');
% % % % end
% % % % end



%MSD means mean squared displacement.From the diffMSD figure we can guess
%the movement of the bacteria.
h4=createMSDfigure(timeResult.Msd(:,1), timeResult.Msd(:,2:4),timeResult.MsdDiff(:,1),timeResult.MsdDiff(:,2:4),clusterDataList,iData);
saveFile6=strcat(dirSave,'\','MSD');
saveas(h4,saveFile6,'fig');
saveas(h4,strcat(saveGraphicFile,'\','MSD'),'pdf');

% screened by Lei Ni
% % % % %correlation and power spectrum
% % % % h5=createCorrelationfigure(timeResult.Correlation(:,1), timeResult.Correlation(:,2:4),timeResult.PowerSpectrum(:,1),timeResult.PowerSpectrum(:,2:3),clusterDataList,iData);
% % % % saveFile7=strcat(dirSave,'\','Correlation');
% % % % saveas(h5,saveFile7,'fig');
% % % % saveas(h5,strcat(saveGraphicFile,'\','Correlation'),'pdf');

if numel(find(clusterDataList==0))~=2
%this figure show how the length of bacteria changes with time 
h6=createLenfigure([timeResult.beforeWdencmpLen,timeResult.afterWdencmpLen],timeResult.diffLen,clusterDataList,iData);
saveFile8=strcat(dirSave,'\','LenChange');
saveas(h6,saveFile8,'fig');
saveas(h6,strcat(saveGraphicFile,'\','LenChange'),'pdf');

%this figure has four parts. In first part,we get the mean value of all
%orientation and make the orientation distribution. In second part, we
%caculate the angle velocity and show their distribution. Then we think
%about the correlation between a row of angle velocity and plot its power
%spectrum.
% Lei Ni screen this figure
% % % % h7=createOritationFigure(timeResult.OritationInfo.Distribution,timeResult.OritationInfo.correlation(:,1),timeResult.OritationInfo.correlation(:,2),timeResult.OritationInfo.power(:,1),timeResult.OritationInfo.power(:,2),timeResult.OritationInfo.angleVelocity,clusterDataList,iData);
% % % % saveFile9=strcat(dirSave,'\','OrientationStatistic');
% % % % saveas(h7,saveFile9,'fig');
% % % % saveas(h7,strcat(saveGraphicFile,'\','OrientationStatistic'),'pdf');
end
end

%% create V1 Histgram
function h=createP1VelocityHistgramfigure(ZData1, YData1, CData1, histResult1,dataList,iData)
%CREATEFIGURE(ZDATA1,YDATA1,CDATA1,HISTRESULT1)
%  ZDATA1:  surface zdata
%  YDATA1:  surface ydata
%  CDATA1:  surface cdata
%  HISTRESULT1:  vector of y data

%  Auto-generated by MATLAB on 16-Jan-2012 10:34:32

%find peaks
[peaks,legs]=findpeaks(histResult1);
peakAll=[];
for i=1:numel(peaks)
    peakAll=strcat(peakAll,'   *',num2str(YData1(legs(i))));
end
% Create figure
iData=strcat('#',num2str(dataList(1)),'-',num2str(dataList(2)),'  number',num2str(iData));
scrsz=get(0,'ScreenSize');
figure1=figure;
set (gcf,'Position',[10 10 scrsz(3) scrsz(4)-80]);
set(gcf, 'PaperPositionMode', 'auto');


% Create axes
axes1 = axes('Parent',figure1,'YTick',[0.001 0.01 0.1 1 10 100],...
    'YScale','log',...
    'YMinorTick','on',...
    'XTick',[0.001 0.01 0.1 1 10 100],...
    'XScale','log',...
    'XMinorTick','on',...
    'Position',[0.399591836734694 0.222028352918767 0.174285714285714 0.304228660924041]);
% Uncomment the following line to preserve the X-limits of the axes
% xlim(axes1,[0.001 100]);
% Uncomment the following line to preserve the Y-limits of the axes
% ylim(axes1,[0.001 100]);
box(axes1,'on');
hold(axes1,'all');

% Create surface
surface('Parent',axes1,'ZData',ZData1,'YData',YData1,'XData',YData1,...
    'CData',CData1);

% Create xlabel
xlabel('|p1V| pixels/frame','FontSize',16);

% Create ylabel
ylabel('|p2V| pixels/frame','FontSize',16);

% Create axes
axes2 = axes('Parent',figure1,'XTick',[0.01 0.1 1 10 100],...
    'XScale','log',...
    'XMinorTick','on',...
    'Position',[0.396875 0.575751761942052 0.177410714285714 0.166613155833986]);
box(axes2,'on');
hold(axes2,'all');

% Create semilogx
semilogx(YData1,histResult1,'Parent',axes2,'Marker','o','LineWidth',2,...
    'DisplayName','histResult.p1(:,2) vs histResult.p1(:,1)');

% Create ylabel
ylabel('Displacement pixels','FontSize',14);

% Create colorbar
colorbar('peer',axes1,...
    [0.585170068027211 0.218872357086923 0.0131519274376417 0.306969459671104]);

% Create textbox
annotation(figure1,'textbox',...
    [0.563673469387757 0.535413469068128 0.0701754385964912 0.0234925606891151],...
    'String',{'Displacement pixels'},...
    'FitBoxToText','off',...
    'LineStyle','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.409387755102042 0.780519185591229 0.172514619883041 0.0438527799530149],...
    'String',{'P1 Velocity Histogram'},...
    'FontSize',24,...
    'FitBoxToText','off',...
    'LineStyle','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.589122797572435 0.613985496347545 0.0291063690942314 0.0735607675906183],...
    'String',{iData},...
    'FontSize',30,...
    'FitBoxToText','off',...
    'LineStyle','none');

annotation(figure1,'textbox',...
    [0.369256450351837 0.867096474297436 0.280469898358092 0.0245202558635369],...
    'String',{peakAll},...
    'FontSize',14,...
    'FitBoxToText','off',...
    'LineStyle','none');
hold on
stem(YData1(legs),peaks,'--g');
h=gca;
end

%% create V2 Histgram
function h=createP2VelocityHistgramfigure(ZData1, YData1, CData1, histResult1,dataList,iData)
%CREATEFIGURE(ZDATA1,YDATA1,CDATA1,HISTRESULT1)
%  ZDATA1:  surface zdata
%  YDATA1:  surface ydata
%  CDATA1:  surface cdata
%  HISTRESULT1:  vector of y data

%  Auto-generated by MATLAB on 16-Jan-2012 10:34:32
[peaks,legs]=findpeaks(histResult1);
peakAll=[];
for i=1:numel(peaks)
    peakAll=strcat(peakAll,'   *',num2str(YData1(legs(i))));
end

% Create figure
iData=strcat('#',num2str(dataList(1)),'-',num2str(dataList(2)),'  number',num2str(iData));
scrsz=get(0,'ScreenSize');
figure1=figure;
set(gcf, 'PaperPositionMode', 'auto');
set (gcf,'Position',[10 10 scrsz(3) scrsz(4)-80]);

% Create axes
axes1 = axes('Parent',figure1,'YTick',[0.001 0.01 0.1 1 10 100],...
    'YScale','log',...
    'YMinorTick','on',...
    'XTick',[0.001 0.01 0.1 1 10 100],...
    'XScale','log',...
    'XMinorTick','on',...
    'Position',[0.399591836734694 0.222028352918767 0.174285714285714 0.304228660924041]);
% Uncomment the following line to preserve the X-limits of the axes
% xlim(axes1,[0.001 100]);
% Uncomment the following line to preserve the Y-limits of the axes
% ylim(axes1,[0.001 100]);
box(axes1,'on');
hold(axes1,'all');

% Create surface
surface('Parent',axes1,'ZData',ZData1,'YData',YData1,'XData',YData1,...
    'CData',CData1);

% Create xlabel
xlabel('|p1V| pixels/frame','FontSize',16);

% Create ylabel
ylabel('|p2V| pixels/frame','FontSize',16);

% Create axes
axes2 = axes('Parent',figure1,'XTick',[0.001 0.01 0.1 1 10 100],...
    'XScale','log',...
    'XMinorTick','on',...
    'Position',[0.396875 0.575751761942052 0.177410714285714 0.166613155833986]);
box(axes2,'on');
hold(axes2,'all');

% Create semilogx
semilogx(YData1,histResult1,'Parent',axes2,'Marker','o','LineWidth',2,...
    'DisplayName','histResult.p1(:,2) vs histResult.p1(:,1)');

% Create ylabel
ylabel('Displacement pixels','FontSize',14);

% Create colorbar
colorbar('peer',axes1,...
    [0.585170068027211 0.218872357086923 0.0131519274376417 0.306969459671104]);

% Create textbox
annotation(figure1,'textbox',...
    [0.563673469387757 0.535413469068128 0.0701754385964912 0.0234925606891151],...
    'String',{'Displacement pixels'},...
    'FitBoxToText','off',...
    'LineStyle','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.409387755102042 0.780519185591229 0.172514619883041 0.0438527799530149],...
    'String',{'P2 Velocity Histogram'},...
    'FontSize',24,...
    'FitBoxToText','off',...
    'LineStyle','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.589122797572435 0.613985496347545 0.0291063690942314 0.0735607675906183],...
    'String',{iData},...
    'FontSize',30,...
    'FitBoxToText','off',...
    'LineStyle','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.369256450351837 0.867096474297436 0.280469898358092 0.0245202558635369],...
    'String',{peakAll},...
    'FontSize',14,...
    'FitBoxToText','off',...
    'LineStyle','none');
h=gca;
hold on
stem(YData1(legs),peaks,'--g');
end

%% create Polar figure
function h=myPolarPlot(theta,dataList,threshold)
% figure1=figure;
% h=gca;
% myColor=[1,0,0;0,1,0;0,0,1];
% for iData=1:size(dataList,2)
%     x=dataList(:,iData).*cos(theta);
%     y=dataList(:,iData).*sin(theta);
%      plot(x,y, ...
%                 'Marker','.','MarkerSize',6,'Color',myColor(iData,:),'LineStyle','-');hold on
% end
% text1=strcat(' Velocity Orientation  threshold=', num2str(threshold));
% annotation(figure1,'textbox',...
%     [0.409387755102042 0.780519185591229 0.172514619883041 0.0438527799530149],...
%     'String',{text1},...
%     'FontSize',24,...
%     'LineStyle','none');
theta(end+1)=theta(end)+theta(2)-theta(1);
dataList(end+1,:)=dataList(1,:);
[theta,dataList]=makeSpline(theta,dataList);
figure;
set(gcf, 'PaperPositionMode', 'auto');
h1=polar(theta,dataList(:,1));
set(h1,'LineWidth',3.5,...
    'Color',[1 0 0],...
    'DisplayName','Alldata');
hold on
if ~isempty(threshold)
    h2=polar(theta,dataList(:,2));
    set(h2,'LineWidth',2.5,'Color',[0 1 0],...
        'DisplayName','beforeThereshold');
    h3=polar(theta,dataList(:,3));
    set(h3,'LineWidth',2,'Color',[0 0 1],...
        'DisplayName','afterThereshold');
end
axis1=gca;
figure1=gcf;

set(axis1,'view',[270,90]);

% Create legend;
axis1=gca;
legend1 = legend(axis1,'show');
set(legend1,...
    'Position',[0.714224897110486 0.838975683666983 0.278571428571429 0.144444444444444]);

% Create textbox
if ~isempty(threshold)
annotation(figure1,'textbox',...
    [0.0420475575668196 0.168489260464597 0.0834444444444445 0.107200000000001],...
    'String',{'Velocity','Orientation',strcat('threShold=',num2str(threshold))},...
    'FontSize',20,...
    'FitBoxToText','off',...
    'LineStyle','none');
else
    annotation(figure1,'textbox',...
    [0.0420475575668196 0.168489260464597 0.0834444444444445 0.107200000000001],...
    'String',{'Velocity','Orientation','No Threshold'},...
    'FontSize',30,...
    'FitBoxToText','off',...
    'LineStyle','none');
end
h=gcf;
end
function [theta,dataList]=makeSpline(theta0,dataList0)
theta=theta0(1):(theta0(end)-theta0(1))/numel(theta0)/5:theta0(end);
theta=theta';
for i=1:size(dataList0,2)
    dataList(:,i)=spline(theta0,dataList0(:,i),theta);
end
end

%% create MSD figure
function h=createMSDfigure(X1, YMatrix1, X2, YMatrix2,dataList,iData)
%CREATEFIGURE1(TIMERESULT1,YMATRIX1,TIMERESULT2,YMATRIX2)
%  TIMERESULT1:  vector of x data
%  YMATRIX1:  matrix of y data
%  TIMERESULT2:  vector of x data
%  YMATRIX2:  matrix of y data

%  Auto-generated by MATLAB on 18-Feb-2012 20:39:50

% % Create figure
% figure1 = figure('Position',[1 1 scrsz(3) scrsz(4)-100]);
iData=strcat('#',num2str(dataList(1)),'-',num2str(dataList(2)),'  number',num2str(iData));
scrsz=get(0,'ScreenSize');
figure1=figure;
set(gcf, 'PaperPositionMode', 'auto');
set (gcf,'Position',[1 1 scrsz(3) scrsz(4)-100]); 

% Create axes
axes1 = axes('Parent',figure1,'YScale','log','YMinorTick','on',...
    'XScale','log',...
    'XMinorTick','on',...
    'Position',[0.340104166666667 0.563706975327452 0.325 0.33181541273225]);
box(axes1,'on');
hold(axes1,'all');

% Create multiple lines using matrix input to loglog
loglog1 = loglog(X1,YMatrix1,'Parent',axes1);
set(loglog1(1),'Marker','o','DisplayName','p1');
set(loglog1(2),'Marker','square','DisplayName','p2');
set(loglog1(3),'Marker','*','DisplayName','Centroid');

% Create ylabel
ylabel('MSD / pixels^2','FontSize',20);

% Create axes
axes2 = axes('Parent',figure1,...
    'YTickLabel',{'0','0.25','0.5','0.75','1','1.25','1.5','1.75','2','2.25','2.5'},...
    'YTick',[0 0.25 0.5 0.75 1 1.25 1.5 1.75 2 2.25 2.5],...
    'YGrid','on',...
    'XScale','log',...
    'XMinorTick','on',...
    'Position',[0.343229166666667 0.193024672555591 0.323958333333333 0.328571428571433]);
% Uncomment the following line to preserve the Y-limits of the axes
% ylim(axes2,[0 2.5]);
box(axes2,'on');
hold(axes2,'all');
axis([0 100 0 2.5])
axis 'auto x'

% Create xlabel
xlabel('time delay / frame','FontSize',20);

% Create ylabel
ylabel('slope','FontSize',20);

% Create multiple lines using matrix input to semilogx
semilogx1 = semilogx(X2,YMatrix2,'Parent',axes2);
set(semilogx1(1),'Marker','o','DisplayName','timeResult.MsdDiff(:,2)');
set(semilogx1(2),'Marker','square','Color',[0 0.498039215803146 0],...
    'DisplayName','timeResult.MsdDiff(:,3)');
set(semilogx1(3),'Marker','*','DisplayName','timeResult.MsdDiff(:,4)');

% Create legend
legend1 = legend(axes1,'show');
set(legend1,'EdgeColor',[1 1 1],'YColor',[1 1 1],'XColor',[1 1 1],...
    'Position',[0.360169923554111 0.777097877897138 0.0585800764458893 0.0896398619749296]);

% Create textbox
annotation(figure1,'textbox',...
    [0.617205168776371 0.518359081461742 0.0709520319786809 0.117136659436009],...
    'String',{iData},...
    'FontSize',30,...
    'FitBoxToText','off',...
    'LineStyle','none');
h=gcf;
end

%% create Correlation figure
function h=createCorrelationfigure(X1, YMatrix1, X2, YMatrix2,dataList,iData)
iData=strcat('#',num2str(dataList(1)),'-',num2str(dataList(2)),'  number',num2str(iData));
%CREATEFIGURE(TIMERESULT1,YMATRIX1,TIMERESULT2,YMATRIX2)
%  TIMERESULT1:  vector of x data
%  YMATRIX1:  matrix of y data
%  TIMERESULT2:  vector of x data
%  YMATRIX2:  matrix of y data

%  Auto-generated by MATLAB on 19-Feb-2012 23:05:18

% Create figure
scrsz=get(0,'ScreenSize');
figure1=figure;
set(gcf, 'PaperPositionMode', 'auto');
set (gcf,'Position',[1 1 scrsz(3) scrsz(4)-100]); 

% Create axes
axes1 = axes('Parent',figure1,'XScale','log','XMinorTick','on',...
    'Position',[0.352604166666667 0.559541758655969 0.3265625 0.387302110817942]);
box(axes1,'on');
hold(axes1,'all');

% Create multiple lines using matrix input to semilogx
semilogx1 = semilogx(X1,YMatrix1,'Parent',axes1);
set(semilogx1(1),'DisplayName','<vp1,vp1>');
set(semilogx1(2),'DisplayName','<vp2,vp2>');
set(semilogx1(3),'DisplayName','<vp1,vp2>');

% Create xlabel
xlabel('TimeDelay/Frames','FontSize',20);

% Create ylabel
ylabel('Correlation','FontSize',20);

% Create legend
legend1 = legend(axes1,'show');
set(legend1,'EdgeColor',[1 1 1],'YColor',[1 1 1],'XColor',[1 1 1],...
    'Position',[0.603724558069551 0.83410669265801 0.0616898148148148 0.0802829980181316]);

% Create axes
axes2 = axes('Parent',figure1,'XScale','log','XMinorTick','on',...
    'Position',[0.352604166666667 0.0914175531998879 0.328645833333333 0.387302110817942]);
box(axes2,'on');
hold(axes2,'all');

% Create multiple lines using matrix input to semilogx
semilogx2 = semilogx(X2,YMatrix2,'Parent',axes2);
set(semilogx2(1),'DisplayName','vp1.PowerSpectrum');
set(semilogx2(2),'DisplayName','vp2.PowerSpectrum');

% Create xlabel
xlabel('/Frame','FontSize',20);

% Create ylabel
ylabel('Power','FontSize',20);

% Create legend
legend2 = legend(axes2,'show');
set(legend2,'EdgeColor',[1 1 1],'YColor',[1 1 1],'XColor',[1 1 1],...
    'Position',[0.580124683734031 0.374917104564533 0.0913773148148148 0.0745792891359728]);

% Create textbox
annotation(figure1,'textbox',...
    [0.398286119812059 0.800161435346404 0.075742364917776 0.113006396588486],...
    'String',{iData},...
    'FontSize',30,...
    'FitBoxToText','off',...
    'LineStyle','none');
h=gcf;
end

%% here create the bacteria length figure
function h=createLenfigure(YMatrix1, Y1,dataList,iData)
%CREATEFIGURE(YMATRIX1,Y1)
%  YMATRIX1:  matrix of y data
%  Y1:  vector of y data

%  Auto-generated by MATLAB on 05-Mar-2012 22:22:06

% Create figure
iData=strcat('#',num2str(dataList(1)),'-',num2str(dataList(2)),'  number',num2str(iData));
scrsz=get(0,'ScreenSize');
figure1=figure;
set (gcf,'Position',[10 10 scrsz(3) scrsz(4)-80]);
set(gcf, 'PaperPositionMode', 'auto');

% Create axes
axes1 = axes('Parent',figure1,...
    'Position',[0.336664040616246 0.566938786311604 0.34001663165266 0.4075]);
box(axes1,'on');
hold(axes1,'all');

% Create multiple lines using matrix input to plot
plot1 = plot(YMatrix1,'Parent',axes1);
set(plot1(1),'DisplayName','beforeSmooth');
set(plot1(2),'Color',[1 0 0],'DisplayName','afterSmooth');

% Create ylabel
ylabel('bacteria length','FontSize',20);

% Create legend
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.349281230560334 0.876536477441454 0.0735443561507648 0.0652186681573255]);

% Create axes
axes2 = axes('Parent',figure1,'YScale','log','YMinorTick','on',...
    'Position',[0.338230917366946 0.117013081201402 0.341058298319328 0.406103851444292]);
box(axes2,'on');
hold(axes2,'all');

% Create semilogy
semilogy(Y1,'Parent',axes2);

% Create xlabel
xlabel('Time/frame','FontSize',20);

% Create ylabel
ylabel('length change','FontSize',20);

% Create textbox
annotation(figure1,'textbox',...
    [0.357348587300379 0.440167828755512 0.0449249676584734 0.0276073619631911],...
    'String',{iData},...
    'FontSize',30,...
    'FitBoxToText','off',...
    'LineStyle','none');
h=gca;
end

%% here create the Oritation figure
function h=createOritationFigure(Distribution, X1, Y1, X2, Y2,angleVelocity,dataList,iData)
% Create figure
iData=strcat('#',num2str(dataList(1)),'-',num2str(dataList(2)),'  number',num2str(iData));
scrsz=get(0,'ScreenSize');
figure1=figure;
set (gcf,'Position',[10 10 scrsz(3) scrsz(4)-80]);
set(gcf, 'PaperPositionMode', 'auto');

% Create axes
hist(Distribution,16);
h=gca;
set(h, 'Position',[0.342708333333333 0.707889125799574 0.161095177663425 0.271855010660981],...
    'CLim',[1 2]);
box(h,'on');

% Create ylabel
ylabel('Number','FontSize',20);

% Create xlabel
xlabel('Orientation Distribution','FontSize',20);

% Create axes
axes2 = axes('Parent',figure1,'XScale','log','XMinorTick','on',...
    'Position',[0.344270833333333 0.414712153518124 0.352604166666667 0.238928982835504]);
box(axes2,'on');
hold(axes2,'all');

% Create semilogx
semilogx(X1,Y1,'Parent',axes2);

% Create xlabel
xlabel('Orientation Velocity Correlation','FontSize',20);

% Create axes
axes3 = axes('Parent',figure1,'XScale','log','XMinorTick','on',...
    'Position',[0.3453125 0.0490405117270789 0.352604166666667 0.28997867803838]);
box(axes3,'on');
hold(axes3,'all');

% Create semilogx
semilogx(X2,Y2,'Parent',axes3);

% Create xlabel
xlabel({'/Frame',''},'FontSize',20);

% Create ylabel
ylabel('Orientation Velocity Spectrum','FontSize',20);

% Create axes
axes4 = axes('Parent',figure1,...
    'Position',[0.521875 0.706140511727079 0.174479166666667 0.272233901918977],...
    'CLim',[1 2]);
hist(angleVelocity,16);
h=gca;
box(h,'on');

% Create xlabel
xlabel('Orientation Velocity Distribution','FontSize',20);

% Create textbox
annotation(figure1,'textbox',...
    [0.603177735650846 0.220338129571035 0.0947389310158203 0.0930947062498607],...
    'String',{iData},...
    'FontSize',30,...
    'FitBoxToText','off',...
    'LineStyle','none');
h=gca;
end
