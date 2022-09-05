function JumpTimeHist(NewAllData,saveFile,threShold)
HistforEach(NewAllData,saveFile,threShold,2,3,1)
HistforEach(NewAllData,saveFile,threShold,4,5,2)
end
function HistforEach(NewAllData,saveFile,threShold,c1,c2,num)
jumpTimePeriodAll=[];
jumpAmptitudeAll=[];
meanSpeed=[];
allDisplacement=[]; 
sectionDisplacement=[];
meanVelocity=[];
angle=[];
n=0;
for j=1:numel(NewAllData)
    positionInfo=NewAllData{1,j}.denosieData(:,c1:c2);
    positionInfoBefore=positionInfo(1:end-1,:);
    positionInfoAfter=positionInfo(2:end,:);
    positionChange=[];
    n=n+size(positionInfo,1);
    positionChange(:,1)=((positionInfoAfter(:,1)-positionInfoBefore(:,1)).^2+(positionInfoAfter(:,2)-positionInfoBefore(:,2)).^2).^0.5;
    jumpTimePoint=find(positionChange>=threShold);
 
    %find continues jump and statistic
    Matrix1=false(size(positionChange));
    Matrix1(jumpTimePoint)=1;
    Matrix2=logical(1-Matrix1);
    cc1=regionprops(Matrix1,positionChange,'PixelValues','MeanIntensity','PixelIdxList');
    cc2=regionprops(Matrix2,positionChange,'Area','MeanIntensity','PixelValues','PixelIdxList');
    if Matrix1(1)==1
        for i=1:numel(cc1)-1
            jumpTimePeriodAll=[jumpTimePeriodAll;cc2(i,1).Area];
            meanSpeed=[meanSpeed;cc2(i,1).MeanIntensity];
            allDisplacement=[allDisplacement;sum(cc2(i,1).PixelValues)];
            sectionDisplacement=[sectionDisplacement;sqrt((positionInfo(cc2(i).PixelIdxList(1),1)-positionInfo(cc2(i).PixelIdxList(end)+1,1)).^2+(positionInfo(cc2(i).PixelIdxList(1),2)-positionInfo(cc2(i).PixelIdxList(end)+1,2)).^2)];
            meanVelocity=[meanVelocity;sectionDisplacement(end)/numel(cc2(i).PixelIdxList)];
            vector2=[positionInfo(cc2(i).PixelIdxList(end)+1,1)-positionInfo(cc2(i).PixelIdxList(1),1),positionInfo(cc2(i).PixelIdxList(end)+1,2)-positionInfo(cc2(i).PixelIdxList(1),2)];
            vector1=[NewAllData{1,j}.denosieData(cc2(i).PixelIdxList(1),2)-NewAllData{1,j}.denosieData(cc2(i).PixelIdxList(1),4),NewAllData{1,j}.denosieData(cc2(i).PixelIdxList(1),3)-NewAllData{1,j}.denosieData(cc2(i).PixelIdxList(1),5)];
            angle=[angle;angleBetweenTwoVector(vector1,vector2)];
        end
    else
        for i=1:numel(cc1)-1
            jumpTimePeriodAll=[jumpTimePeriodAll;cc2(i+1,1).Area];
            meanSpeed=[meanSpeed;cc2(i+1,1).MeanIntensity];
            allDisplacement=[allDisplacement;sum(cc2(i+1,1).PixelValues)];
            sectionDisplacement=[sectionDisplacement;sqrt((positionInfo(cc2(i+1).PixelIdxList(1),1)-positionInfo(cc2(i+1).PixelIdxList(end)+1,1)).^2+(positionInfo(cc2(i+1).PixelIdxList(1),2)-positionInfo(cc2(i+1).PixelIdxList(end)+1,2)).^2)];
            meanVelocity=[meanVelocity;sectionDisplacement(end)/numel(cc2(i+1).PixelIdxList)];
            vector2=[positionInfo(cc2(i+1).PixelIdxList(end)+1,1)-positionInfo(cc2(i+1).PixelIdxList(1),1),positionInfo(cc2(i+1).PixelIdxList(end)+1,2)-positionInfo(cc2(i+1).PixelIdxList(1),2)];
            vector1=[NewAllData{1,j}.denosieData(cc2(i+1).PixelIdxList(1),2)-NewAllData{1,j}.denosieData(cc2(i+1).PixelIdxList(1),4),NewAllData{1,j}.denosieData(cc2(i+1).PixelIdxList(1),3)-NewAllData{1,j}.denosieData(cc2(i+1).PixelIdxList(1),5)];
            angle=[angle;angleBetweenTwoVector(vector1,vector2)];
        end
    end
    for i=1:numel(cc1)
        jumpAmptitudeAll=[jumpAmptitudeAll;sum(cc1(i,1).PixelValues)];
        meanSpeed=[meanSpeed;cc1(i,1).MeanIntensity];
        allDisplacement=[allDisplacement;sum(cc1(i,1).PixelValues)];
        sectionDisplacement=[sectionDisplacement;sqrt((positionInfo(cc1(i).PixelIdxList(1),1)-positionInfo(cc1(i).PixelIdxList(end)+1,1)).^2+(positionInfo(cc1(i).PixelIdxList(1),2)-positionInfo(cc1(i).PixelIdxList(end)+1,2)).^2)];
        meanVelocity=[meanVelocity;sectionDisplacement(end)/numel(cc1(i).PixelIdxList)];
        vector2=[positionInfo(cc1(i).PixelIdxList(end)+1,1)-positionInfo(cc1(i).PixelIdxList(1),1),positionInfo(cc1(i).PixelIdxList(end)+1,2)-positionInfo(cc1(i).PixelIdxList(1),2)];
        vector1=[NewAllData{1,j}.denosieData(cc1(i).PixelIdxList(1),2)-NewAllData{1,j}.denosieData(cc1(i).PixelIdxList(1),4),NewAllData{1,j}.denosieData(cc1(i).PixelIdxList(1),3)-NewAllData{1,j}.denosieData(cc1(i).PixelIdxList(1),5)];
        angle=[angle;angleBetweenTwoVector(vector1,vector2)];
    end   
end
% disp(n)
% disp(numel(meanSpeed))
if numel(jumpTimePeriodAll)>=100
    [theta,dataList]=makeDataList(angle,sectionDisplacement,meanVelocity,threShold);
    h=myPolarPlot(theta',dataList,threShold);
    saveas(h,strcat(saveFile,'\caseOrientation',num2str(num)),'fig')
saveas(h,strcat(saveFile,'\Graphic','\caseOrientation',num2str(num)),'pdf')
% for i=1:ceil((max(jumpTimePeriodAll)-min(jumpTimePeriodAll))/40)
%     jumpTimeHist(i,2)=numel(find(jumpTimePeriodAll>=40*(i-1)+min(jumpTimePeriodAll)&jumpTimePeriodAll<40*(i)+min(jumpTimePeriodAll)))./n;
%     jumpTimeHist(i,1)=40*(i-1)+5;
% end
[jumpTimeHist(:,2),jumpTimeHist(:,1)]=hist(jumpTimePeriodAll,logspace(0,3,10));
jumpTimeHist(:,2)=jumpTimeHist(:,2)./n;
jumpTimeHist2=jumpTimeHist(:,2);
jumpTimeHist(jumpTimeHist==0)=min(jumpTimeHist2(jumpTimeHist2~=0));
p1=polyfit(log(jumpTimeHist(:,1)),log(jumpTimeHist(:,2)),1);
jumpTimeLinear=exp(polyval(p1,log(jumpTimeHist(:,1))));

% for i=1:ceil((max(jumpAmptitudeAll)-min(jumpAmptitudeAll))/threShold*2)
%     jumpAmptitudeHist(i,2)=numel(find(jumpAmptitudeAll>=threShold/2*(i-1)+min(jumpAmptitudeAll)&jumpAmptitudeAll<threShold/2*(i)+min(jumpAmptitudeAll)))./n;
%     jumpAmptitudeHist(i,1)=threShold/2*(i-1)+threShold/6;
% end
[jumpAmptitudeHist(:,2),jumpAmptitudeHist(:,1)]=hist(jumpAmptitudeAll,logspace(-0.09691,1.176,6));
jumpAmptitudeHist(:,2)=jumpAmptitudeHist(:,2)./n;
jumpAmptitudeHist2=jumpAmptitudeHist(:,2);
jumpAmptitudeHist(jumpAmptitudeHist==0)=min(jumpAmptitudeHist2(jumpAmptitudeHist2~=0));
p2=polyfit(log(jumpAmptitudeHist(:,1)),log(jumpAmptitudeHist(:,2)),1);
jumpAmptitudeLinear=exp(polyval(p2,log(jumpAmptitudeHist(:,1))));

createfigure(jumpTimeHist(:,1),[jumpTimeHist(:,2),jumpTimeLinear],jumpAmptitudeHist(:,1),[jumpAmptitudeHist(:,2),jumpAmptitudeLinear],p1,p2)
h=gca;
saveas(h,strcat(saveFile,'\jumpHist',num2str(num)),'fig')
saveas(h,strcat(saveFile,'\Graphic','\jumpHist',num2str(num)),'pdf')
createhist(meanVelocity,n)
h=gca;
saveas(h,strcat(saveFile,'\velocityhista-',num2str(num)),'fig')
saveas(h,strcat(saveFile,'\Graphic','\velocityHista-',num2str(num)),'pdf')
createhist2(allDisplacement,meanVelocity,n)
h=gca;
saveas(h,strcat(saveFile,'\velocityhistb-',num2str(num)),'fig')
saveas(h,strcat(saveFile,'\Graphic','\velocityHistb-',num2str(num)),'pdf')
createhist3(allDisplacement,meanVelocity,n)
h=gca;
saveas(h,strcat(saveFile,'\velocityhistc-',num2str(num)),'fig')
saveas(h,strcat(saveFile,'\Graphic','\velocityHistc-',num2str(num)),'pdf')
createhist4(sectionDisplacement,meanVelocity,n)
h=gca;
saveas(h,strcat(saveFile,'\velocityhistd-',num2str(num)),'fig')
saveas(h,strcat(saveFile,'\Graphic','\velocityHistd-',num2str(num)),'pdf')
end
end

function createfigure(X1, YMatrix1, X2, YMatrix2,p1,p2)
%CREATEFIGURE(X1,YMATRIX1,X2,YMATRIX2)
%  X1:  vector of x data
%  YMATRIX1:  matrix of y data
%  X2:  vector of x data
%  YMATRIX2:  matrix of y data

%  Auto-generated by MATLAB on 25-Mar-2012 22:02:35

% Create figure
scrsz=get(0,'ScreenSize');
figure1=figure;
set (gcf,'Position',[10 10 scrsz(3) scrsz(4)-80]);
set(gcf, 'PaperPositionMode', 'auto');

% Create axes
axes1 = axes('Parent',figure1,'YScale','log','YMinorTick','on',...
    'XScale','log',...
    'XMinorTick','on',...
    'Position',[0.327244258872651 0.555555555555556 0.347755741127349 0.395834623261141]);
box(axes1,'on');
hold(axes1,'all');

% Create multiple lines using matrix input to loglog
loglog1 = loglog(X1,YMatrix1,'Parent',axes1);
set(loglog1(1),'MarkerSize',8,'Marker','o','LineWidth',2,'LineStyle','none',...
    'DisplayName','jumpTime');
set(loglog1(2),'LineWidth',3,'Color',[1 0 0],'DisplayName','linear fit');

% Create xlabel
xlabel('JumpTime(frames)','FontSize',20);

% Create ylabel
ylabel('frequency','FontSize',20);

% Create legend
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.60573144977646 0.877887016782087 0.0600986227824463 0.0466836449019125]);

% Create axes
axes2 = axes('Parent',figure1,'YScale','log','YMinorTick','on',...
    'XScale','log',...
    'XMinorTick','on',...
    'Position',[0.327766179540709 0.105505812561027 0.349837987125958 0.375076687116564]);
box(axes2,'on');
hold(axes2,'all');

% Create multiple lines using matrix input to loglog
loglog2 = loglog(X2,YMatrix2,'Parent',axes2);
set(loglog2(1),'MarkerSize',8,'Marker','o','LineWidth',2,'LineStyle','none',...
    'DisplayName','Amptitude');
set(loglog2(2),'LineWidth',3,'Color',[1 0 0],'DisplayName','linear fit');

% Create xlabel
xlabel('Amplitude(pixels)','FontSize',20);

% Create ylabel
ylabel('frequency','FontSize',20);

% Create legend
legend2 = legend(axes2,'show');
set(legend2,...
    'Position',[0.606628601118592 0.397057769265071 0.0614495798319328 0.0454717122720601]);

% Create textbox
annotation(figure1,'textbox',...
    [0.343960767839167 0.558071875556538 0.0874886452623336 0.0567484662576701],...
    'String',{strcat('slope=',num2str(p1(1)))},...
    'FontSize',20,...
    'FitBoxToText','off',...
    'LineStyle','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.347886659259504 0.122162337229639 0.0741761942051684 0.0506134969325154],...
    'String',{strcat('slope=',num2str(p2(1)))},...
    'FontSize',20,...
    'FitBoxToText','off',...
    'LineStyle','none');
end
function createhist(meanSpeed,n)
[num,xout]=hist(meanSpeed,logspace(-3,2,50));
num=num/n;
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1,'XScale','log','XMinorTick','on');
box(axes1,'on');
hold(axes1,'all');

% Create semilogx
semilogx(xout,num);

% Create xlabel
xlabel('section mean Velocity','FontSize',20);

% Create ylabel
ylabel('number of the case/all frame','FontSize',20);
end
function createhist2(allDisplacement,meanSpeed,n)
t=logspace(-3,2,50);
for i=1:numel(t)-1
    p=find(meanSpeed>t(i)&meanSpeed<=t(i+1));
    num(i)=sum(allDisplacement(p));
end
disp(strcat('allDisplacement/n is =',num2str(sum(allDisplacement)/n)))
num=num/sum(allDisplacement);
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1,'XScale','log','XMinorTick','on');
box(axes1,'on');
hold(axes1,'all');

% Create semilogx
semilogx(t(2:end),num);

% Create xlabel
xlabel('section mean velocity','FontSize',20);

% Create ylabel
ylabel('distance travelled during one section/all distance','FontSize',20);
end
function createhist3(allDisplacement,meanSpeed,n)
t=logspace(-3,2,50);
for i=1:numel(t)-1
    p=find(meanSpeed>t(i)&meanSpeed<=t(i+1));
    num(i)=sum(allDisplacement(p));
end
num=num/n;
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1,'XScale','log','XMinorTick','on');
box(axes1,'on');
hold(axes1,'all');

% Create semilogx
semilogx(t(2:end),num);

% Create xlabel
xlabel('section mean velocity','FontSize',20);

% Create ylabel
ylabel('distance travelled during one section/all frame','FontSize',20);
end
function createhist4(sectionDisplacement,meanSpeed,n)
t=logspace(-3,2,50);
for i=1:numel(t)-1
    p=find(meanSpeed>t(i)&meanSpeed<=t(i+1));
    num(i)=sum(sectionDisplacement(p));
end
num=num/n;
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1,'XScale','log','XMinorTick','on');
box(axes1,'on');
hold(axes1,'all');

% Create semilogx
semilogx(t(2:end),num);

% Create xlabel
xlabel('section mean velocity','FontSize',20);

% Create ylabel
ylabel('section displacement/all frame','FontSize',20);
end
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
function angle=angleBetweenTwoVector(vector1,vector2) %return the angle ([-180,180]) betwwen the vector1 and vector2

angle1=atan(vector1(:,2)./abs(vector1(:,1))).*(180/pi);
angle1(vector1(:,1)<0&vector1(:,2)>=0)= 180-angle1(vector1(:,1)<0&vector1(:,2)>=0);
angle1(vector1(:,1)<0&vector1(:,2)<0)= -180-angle1(vector1(:,1)<0&vector1(:,2)<0);

angle2=atan(vector2(:,2)./abs(vector2(:,1))).*(180/pi);
angle2(vector2(:,1)<0&vector2(:,2)>=0)= 180-angle2(vector2(:,1)<0&vector2(:,2)>=0);
angle2(vector2(:,1)<0&vector2(:,2)<0)= -180-angle2(vector2(:,1)<0&vector2(:,2)<0);

angle=angle2- angle1;
angle(angle>180)=angle(angle>180)-360;
angle(angle<-180)=angle(angle<-180)+360;
end
function [theta,dataList]=makeDataList(angle,sectionDisplacement,meanVelocity,threShold)
theta=linspace(-180,180,40);
for i=1:numel(theta)-1
    p=find(angle>=theta(i)&angle<theta(i+1));
    dataList(i,1)=sum(sectionDisplacement(p));
end
p=find(meanVelocity<threShold);
angle1=angle(p);
sectionDisplacement1=sectionDisplacement(p);
for i=1:numel(theta)-1
    p=find(angle1>=theta(i)&angle1<theta(i+1));
    dataList(i,2)=sum(sectionDisplacement1(p));
end
p=find(meanVelocity>=threShold);
angle2=angle(p);
sectionDisplacement2=sectionDisplacement(p);
for i=1:numel(theta)-1
    p=find(angle2>=theta(i)&angle2<theta(i+1));
    dataList(i,3)=sum(sectionDisplacement2(p));
end
theta=theta+360/80;
theta=theta/180*pi;
theta(end)=[];
theta=theta';
end