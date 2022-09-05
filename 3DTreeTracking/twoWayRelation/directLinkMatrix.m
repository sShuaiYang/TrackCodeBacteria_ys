function [h1,result]=directLinkMatrix(bacteriaFrameInfo,bioTree)
% generate direct link matrix
% aveLen=40;
multi=3;
aimFrame=numel(bacteriaFrameInfo);
finalInfo=bacteriaFrameInfo{aimFrame};
% aveLen=50;
aveLen=mean(finalInfo.lengthInfo(finalInfo.lengthInfo<=80));
distMatrix=pdist2(finalInfo.centroidInfo,finalInfo.centroidInfo);
% distMatrix=distMatrix/aveLen;
% distMatrix(distMatrix>multi)=0;
% distMatrix=1./distMatrix;
% distMatrix(abs(distMatrix)==inf)=0;
% distMatrix(distMatrix>=1)=1;

distMatrix(distMatrix<=aveLen*multi & distMatrix~=0)=1;
distMatrix(distMatrix~=1)=0;
distMatrix=logical(distMatrix);

%degree distribution
% degree=[];
% for i=1:size(distMatrix,1)
% iLine=distMatrix(i,:);
% degree(i)=sum(iLine(iLine~=0));
% end
% easyHist(degree,color);

% cluster coefficient
% degree=[];
% linkNum=[];
% cc=[];
% for i=1:size(distMatrix,2)
%     iLine=distMatrix(i,:);
%     linkNum(i)=numel(iLine(iLine~=0));
%     degree(i)=sum(iLine(iLine~=0));
%     aimOrder=1:size(distMatrix,2);
%     aimOrder=aimOrder(iLine~=0);
%     sumCC=0;
%     for u=1:numel(aimOrder)
%         for j=1:numel(aimOrder)
%             if u~=j && distMatrix(aimOrder(u),aimOrder(j))~=0
%                 sumCC=sumCC+(iLine(aimOrder(u))+iLine(aimOrder(j)))/2;
%             end
%         end
%     end
%     cc(i)=sumCC/(linkNum(i)-1)/degree(i);
% end
% for i=1:size(distMatrix,2)
%     iLine=distMatrix(i,:);
%     linkMatrix=distMatrix(iLine,iLine);
%     clusterCoefficient(i)=(numel(linkMatrix(linkMatrix==1))-degree(i))/(degree(i)*(degree(i)-1));
% end
% maxDegree=max(linkNum);
% for i=1:maxDegree;
%     aimDegree=linkNum==i;
%     meanCC(i)=mean(cc(aimDegree));
% end
% figure;plot(linkNum,cc,'LineStyle','none','Marker','o','Color',[1,0,0])
%  hold on;plot(degree,cc,'LineStyle','none','Marker','o','Color',[0,0,1])

% branch mix
maxColorNum=500;
colorChoose=zeros(maxColorNum,3);
colorChoose(:,1)=0.5;
colorChoose(1:300,:)=colormap(jet(300));
degree=[];
branchIndex=bacteriaFrameInfo{aimFrame}.bacteriaInfo(:,5);
for i=1:size(distMatrix,2)
    iLine=distMatrix(i,:);
    degree(i,1)=numel(iLine(iLine==1));
    attachBranch=branchIndex(iLine);
    bacBranch=branchIndex(i);
    attachBranch=unique(attachBranch);
    attachBranch(attachBranch==bacBranch)=[];
    attachNum(i,1)=numel(attachBranch);
end
% attachNum=attachNum./degree;
attachNum(degree==0)=[];
degree(degree==0)=[];
attachNum=attachNum./degree;
mixDegree=cat(2,attachNum,degree);

% uNiqueMatrix=zeros(15,60);
% for i=1:size(mixDegree,1)
%     uNiqueMatrix(mixDegree(i,1)+1,mixDegree(i,2)+1)=uNiqueMatrix(mixDegree(i,1)+1,mixDegree(i,2)+1)+1;
% end
% data=[];
% for i=1:size(uNiqueMatrix,1)
%     for j=1:size(uNiqueMatrix,2)
%         data=[data;j-1,i-1,uNiqueMatrix(i,j)];
%     end
% end
mixResult=twoDimGaussianForMix(mixDegree);
result.mixResult=mixResult;
h1=createSurfacefigure(mixResult.y',mixResult.x',mixResult.Image/max(max(mixResult.Image))*maxColorNum,colorChoose);

% Rg distribution
% figure;
% branchIndex=bacteriaFrameInfo{end}.bacteriaInfo(:,5);
% uniBranch=unique(branchIndex);
% orderNum=1:size(branchIndex,1);
% [bioTree,branchList,~,~]=myBiograph_new2(bioTree);
% n=0;
% for i=1:max(uniBranch)
%     aimOrder=orderNum(branchIndex==i);
%     if numel(aimOrder)>=2
%         n=n+1;
%         aimBranch=branchList(i,:);
%         allRoot=bioTree{aimBranch(1)}.node{aimBranch(2)}.allRoot;
% %                 degreeRadius(n,2)=allRoot(1,1);
%         degreeRadius(n,2)=numel(aimOrder);
%         icentroid=bacteriaFrameInfo{end}.centroidInfo(aimOrder,:);
%         properCentroid=[mean(icentroid(:,1)),mean(icentroid(:,2))];
%         degreeRadius(n,1)=(mean(pdist2(icentroid,properCentroid).^2)).^0.5;
%     end
% end
% bacteriaLength=40;
% colorChoose=zeros(400,3);
% colorChoose(:,1)=0.5;
% colorChoose(1:300,:)=colormap(jet(300));
% degreeRadius(:,1)=degreeRadius(:,1)/bacteriaLength;
% degreeRad=twoDimGaussianForRg(degreeRadius);
% result.degreeRad=degreeRad;
% h2=createSurfaceRgfigure(degreeRad.y',degreeRad.x',degreeRad.Image/max(max(degreeRad.Image))*400,colorChoose);

% branch-branch mix
% amp=20000;
% colorChoose=zeros(1000,3);
% colorChoose(:,1)=0.5;
% colorChoose(1:300,:)=colormap(jet(300));
% degree=[];
% branchIndex=bacteriaFrameInfo{end}.bacteriaInfo(:,5);
% orderNum=1:numel(branchIndex);
% n=0;
% for i=1:max(branchIndex)
%     aimOrder=orderNum(branchIndex==i);
%     if ~isempty(aimOrder)
%         n=n+1;
%         branchBacNumToMixNum(n,2)=numel(aimOrder);
%         linkNum=[];
%         for j=1:numel(aimOrder)
%             iLink=distMatrix(aimOrder(j),:);
%             linkNum=[linkNum;branchIndex(iLink)];
%         end
%         linkNum=unique(linkNum);
%         linkNum(linkNum==i)=[];
%         branchBacNumToMixNum(n,1)=numel(linkNum);
%     end
% end
% mixResult=twoDimGaussianForMix(mixDegree);
% h=createSurfacefigure(mixResult.y',mixResult.x',mixResult.Image*amp,colorChoose);
end
%% for mix
function result=twoDimGaussianForMix(mixDegree)
xMin=0;
yMin=0;
xMax=1;
% xMax=1;
yMax=60;
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
end
function h=createSurfacefigure(xdata,ydata,zdata,colorChoose)
%CREATEFIGURE(ZDATA1,YDATA1,XDATA1,CDATA1)
%  ZDATA1:  surface zdata
%  YDATA1:  surface ydata
%  XDATA1:  surface xdata
%  CDATA1:  surface cdata

%  Auto-generated by MATLAB on 16-Apr-2013 19:57:26

% Create figure
figure1 = figure;
% Create axes
axes1 = axes('Parent',figure1,...
    'Position',[0.3546875 0.29847701266569 0.277395833333333 0.483895130968695],...
    'FontSize',20,...
    'FontName','Times New Roman');
xlim(axes1,[0 60]);
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
xlabel('number of bacteria nearby','FontSize',20,...
    'FontName','Times New Roman');

% Create ylabel
ylabel('different branch num','FontSize',20);
h=gcf;
end

%% for Rg
function result=twoDimGaussianForRg(mixDegree)
xMin=0;
yMin=0;
% xMax=30;
% yMax=21000;
% sigma1=1;
% interval1=0.1;
% sigma2=700;
% interval2=70;
xMax=30;
yMax=60;
sigma1=1;
interval1=0.1;
sigma2=3;
interval2=0.3;
rowPos=xMin:interval1:xMax;
rowPos=rowPos';
columnPos=yMin:interval2:yMax;
resultImage=zeros(numel(rowPos),numel(columnPos));
for i=1:size(mixDegree,1)
%     disp(i)
    rowGauss=1/(2*pi)^0.5/sigma1*exp(-(rowPos-mixDegree(i,1)).^2/2/sigma1^2);
    columnGauss=1/(2*pi)^0.5/sigma2*exp(-(columnPos-mixDegree(i,2)).^2/2/sigma2^2);
    resultImage=resultImage+rowGauss*columnGauss*mixDegree(i,2);
end
resultImage=resultImage/sum(mixDegree(:,2));
result.Image=resultImage;
result.numPoint=size(mixDegree,1);
result.x=rowPos';
result.y=columnPos;
end
function h=createSurfaceRgfigure(xdata,ydata,zdata,colorChoose)
%CREATEFIGURE(ZDATA1,YDATA1,XDATA1,CDATA1)
%  ZDATA1:  surface zdata
%  YDATA1:  surface ydata
%  XDATA1:  surface xdata
%  CDATA1:  surface cdata

%  Auto-generated by MATLAB on 16-Apr-2013 19:57:26

% Create figure
figure1 = figure;
% Create axes
axes1 = axes('Parent',figure1);
% xlim(axes1,[0 12000]);
% ylim(axes1,[0 30]);
xlim(axes1,[0 60]);
ylim(axes1,[0 30]);
hold(axes1,'all');

% Create surface
surface('Parent',axes1,'ZData',zdata,'YData',ydata,'XData',xdata,'CData',zdata,...
    'LineStyle','none');
colorChoose=colorChoose(1:ceil(max(max(zdata))),:);
colormap(colorChoose)
% Create title
title(strcat('mixNumSurface'));

% Create xlabel
xlabel('degreeNum','FontSize',20);

% Create ylabel
ylabel('radius of gyration(/40pixel)','FontSize',20);
h=gcf;
end