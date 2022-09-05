% function [h1,result]=getMixBranchSurfaceAndPlot(bacteriaFrameInfo)
function [h1,result]=getMixBranchSurfaceAndPlot(clusterTree)
% generate direct link matrix
% aveLen=40;
colorChoose1=colormap(jet(4));
multi=3;
% aimFrame=numel(bacteriaFrameInfo);
n=0;
for aimFrame=1200*[3,6,9,12]
    n=n+1;
% finalInfo=bacteriaFrameInfo{aimFrame};
% % aveLen=50;
% aveLen=mean(finalInfo.lengthInfo(finalInfo.lengthInfo<=80));
% distMatrix=pdist2(finalInfo.centroidInfo,finalInfo.centroidInfo);
% distMatrix(distMatrix<=aveLen*multi & distMatrix~=0)=1;
% distMatrix(distMatrix~=1)=0;
% distMatrix=logical(distMatrix);
% for i=1:size(distMatrix,1)
%     distMatrix(i,i)=10000;
% end
distMatrix=clusterTree{aimFrame+1}.distMatrix;
% branch mix
maxColorNum=500;
colorChoose=zeros(maxColorNum,3);
colorChoose(:,1)=0.5;
colorChoose(1:300,:)=colormap(jet(300));
index=find((colorChoose(:,1)==0 & colorChoose(:,2)==0 & colorChoose(:,3)==1)==1);
colorChoose(1:index,1)=linspace(1,0,index);
colorChoose(1:index,2)=linspace(1,0,index);
colorChoose(1:index,3)=1;

colorChoose(1,:)=[1,1,1];
degree=[];
% branchIndex=bacteriaFrameInfo{aimFrame}.bacteriaInfo(:,5);
branchIndex=clusterTree{aimFrame+1}.bacteriaList(:,5);
for i=1:size(distMatrix,2)
    iLine=distMatrix(i,:);
    degree(i,1)=numel(iLine(iLine==1));
    attachBranch=branchIndex(iLine);
    bacBranch=branchIndex(i);
    attachBranch=unique(attachBranch);
    attachBranch(attachBranch==bacBranch)=[];
    attachNum(i,1)=numel(attachBranch);
end
[a,b]=hist(degree);
hold on;plot(b,a)
% % % attachNum=attachNum./degree;
% % % attachNum(degree==0)=0;
% % % mixDegree=cat(2,attachNum,degree);
% % % mixResult=twoDimGaussianForMix(mixDegree);
% % % hold on;plot(mixResult.row(:,1),mixResult.row(:,2),'color',colorChoose1(n,:,:))
end
% result.mixResult=mixResult;
% h1=createSurfacefigure(mixResult.y',mixResult.x',mixResult.Image/max(max(mixResult.Image))*maxColorNum,colorChoose,mixResult.row(:,1),mixResult.row(:,2),mixResult.column(:,1),mixResult.column(:,2));
end
%% for mix
function result=twoDimGaussianForMix(mixDegree)
xMin=0;
yMin=0;
xMax=1;
% xMax=1;
yMax=70;
sigma1=0.01;
interval1=0.001;
sigma2=0.5;
interval2=0.05;
rowPos=xMin:interval1:xMax;
rowPos=rowPos';
columnPos=yMin:interval2:yMax;
resultImage=zeros(numel(rowPos),numel(columnPos));
for i=1:size(mixDegree,1)
%     disp(i)
    rowGauss=1/(2*pi)^0.5/sigma1*exp(-(rowPos-mixDegree(i,1)).^2/2/sigma1^2);
    columnGauss=1/(2*pi)^0.5/sigma2*exp(-(columnPos-mixDegree(i,2)).^2/2/sigma2^2);
    resultImage=resultImage+rowGauss*columnGauss;
end
resultImage=resultImage/size(mixDegree,1);
result.Image=resultImage;
result.numPoint=size(mixDegree,1);
result.x=rowPos';
result.y=columnPos;

neighbour=mixDegree(:,2);
ratio=mixDegree(:,1);
t=0;
for i=0.0125:0.025:1
    t=t+1;
    rowVectorNew(t)=sum(neighbour(logical(ratio<i+0.0125 & ratio>= i-0.0125)));
    if i==0.9875
        rowVectorNew(t)=sum(neighbour(logical(ratio<=i+0.0125 & ratio>=i-0.0125)));
    end
end
rowVectorNew=rowVectorNew/sum(rowVectorNew);
result.row=cat(2,(0.0125:0.025:1)',rowVectorNew');
t=0;
for i=1:2:70
    t=t+1;
    columnVectorNew(t)=sum(ratio(logical(neighbour<i+1 & neighbour>=i-1)));
    if i==69
        columnVectorNew(t)=sum(ratio(logical(neighbour<=i+1 & neighbour>=i-1)));
    end
end
columnVectorNew=columnVectorNew.*(1:2:70);
columnVectorNew=columnVectorNew/sum(columnVectorNew);
result.column=cat(2,(1:2:70)',columnVectorNew');
end


function h=createSurfacefigure(xdata,ydata,zdata,colorChoose,X1,Y1,X2,Y2)
% Create figure
figure1 = figure;
% Create axes
axes1 = axes('Parent',figure1,'YTick',zeros(1,0),'XTick',zeros(1,0),...
    'Position',[0.35 0.31 0.30 0.60],...
    'FontSize',20,...
    'FontName','Times New Roman');
xlim(axes1,[0 70]);
box(axes1,'on');
% ylim(axes1,[0 15]);
ylim(axes1,[0 1]);
hold(axes1,'all');
% Create surface
surface('Parent',axes1,'ZData',zdata,'YData',ydata,'XData',xdata,'CData',zdata,...
    'LineStyle','none');
colorChoose=colorChoose(1:ceil(max(max(zdata))),:);
colormap(colorChoose)
% Create title
title(strcat('mixNumSurface'),'FontSize',20,'FontName','Times New Roman');

% Create xlabel
h=gcf;

% Create axes
axes2 = axes('Parent',figure1,...
    'Position',[0.25 0.31 0.10 0.60],...
    'FontSize',14);
view(axes2,[-90 90]);
box(axes2,'on');
hold(axes2,'all');

% Create xlabel
xlabel('ratio of different branch','VerticalAlignment','bottom',...
    'Rotation',90,...
    'HorizontalAlignment','center',...
    'FontSize',20);

% Create bar
bar(X1,Y1,'FaceColor',[0.47843137383461 0.062745101749897 0.894117653369904],...
    'EdgeColor','none',...
    'LineWidth',2,...
    'Parent',axes2);

% Create axes
axes3 = axes('Parent',figure1,...
    'Position',[0.35 0.16 0.30 0.15],...
    'FontSize',14);
box(axes3,'on');
hold(axes3,'all');

% Create xlabel
xlabel('neighbour bacteria number','FontSize',20);

% Create bar
bar(X2,Y2,'FaceColor',[0.47843137383461 0.062745101749897 0.894117653369904],...
    'EdgeColor','none',...
    'LineWidth',2,...
    'Parent',axes3);

end
function h=createSurfacefigure1(xdata,ydata,zdata,colorChoose,X1,Y1,X2,Y2)
% Create figure
figure1 = figure;
% Create axes
axes1 = axes('Parent',figure1,'YTick',zeros(1,0),'XTick',zeros(1,0),...
    'Position',[0.35 0.311111111111111 0.2921875 0.616346270100348],...
    'FontSize',20,...
    'FontName','Times New Roman');
xlim(axes1,[0 70]);
box(axes1,'on');
% ylim(axes1,[0 15]);
ylim(axes1,[0 1]);
hold(axes1,'all');
% Create surface
surface('Parent',axes1,'ZData',zdata,'YData',ydata,'XData',xdata,'CData',zdata,...
    'LineStyle','none');
colorChoose=colorChoose(1:ceil(max(max(zdata))),:);
colormap(colorChoose)
% Create title
title(strcat('mixNumSurface'),'FontSize',20,'FontName','Times New Roman');

% Create xlabel
h=gcf;

% Create axes
axes2 = axes('Parent',figure1,...
    'Position',[0.2578125 0.310941782436898 0.0969067879499213 0.616152885680621],...
    'FontSize',14);
view(axes2,[-90 90]);
box(axes2,'on');
hold(axes2,'all');

% Create xlabel
xlabel('ratio of different branch','VerticalAlignment','bottom',...
    'Rotation',90,...
    'HorizontalAlignment','center',...
    'FontSize',20);

% Create bar
bar(X1,Y1,'FaceColor',[0.47843137383461 0.062745101749897 0.894117653369904],...
    'EdgeColor','none',...
    'LineWidth',2,...
    'Parent',axes2);

% Create axes
axes3 = axes('Parent',figure1,...
    'Position',[0.354687500000001 0.127312295973885 0.292708333333333 0.183764921702293],...
    'FontSize',14);
box(axes3,'on');
hold(axes3,'all');

% Create xlabel
xlabel('neighbour bacteria number','FontSize',20);

% Create bar
bar(X2,Y2,'FaceColor',[0.47843137383461 0.062745101749897 0.894117653369904],...
    'EdgeColor','none',...
    'LineWidth',2,...
    'Parent',axes3);

end