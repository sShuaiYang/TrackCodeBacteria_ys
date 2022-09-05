function  gRAndDivTimePlot(gRDiv,scale,dirSave)
% gRAndDivTimeAnalysis获得gRDiv
%然后进行plot
% gRDiv
%第1列 field
%第2列 代数
%第3列 出生菌长
%第4列 分裂时菌长
%第5列 出生时面积
%第6列 分裂时面积
%第7列 生长率
%第8列 division time

genrLogical=false(max(gRDiv(:,2))-1,size(gRDiv,1));
for iGenr=2:max(gRDiv(:,2))
    genrLogical(iGenr-1,:)=gRDiv(:,2)==iGenr;    
end
cmap=colormap(lines(10));

figure,
xData=gRDiv(:,3)*scale;
yData=gRDiv(:,4)*scale;
for iGenr=2:max(gRDiv(:,2))
    scatter(xData(genrLogical(iGenr-1,:)),yData(genrLogical(iGenr-1,:)),50,cmap(iGenr,:),'filled','MarkerFaceAlpha',0.15);
    hold on
end
z = bindata(xData, yData,15,[0.8,6]);
hold on
plot(z(1,:),z(2,:),'MarkerSize',8,'Marker','square','LineWidth',2,'Color',[1 0 0]);
xlabel('LI (birth Length, um)');
ylabel('LF (divison Length, um)');
set(gca,'FontSize',12)
xlim(gca,[0 6]);
ylim(gca,[0 12]);
hold on, line([0,6],[0,12],'Color','black')
fname='LI vs LF';
saveas(gcf,[dirSave,'\',fname,'.fig']);
saveas(gcf,[dirSave,'\',fname,'.tif']);
hold off
close all


figure,
xData=gRDiv(:,3)*scale;
yData=(gRDiv(:,4)-gRDiv(:,3))*scale;
for iGenr=2:max(gRDiv(:,2))    
    scatter(xData(genrLogical(iGenr-1,:)),yData(genrLogical(iGenr-1,:)),50,cmap(iGenr,:),'filled','MarkerFaceAlpha',0.15);
    hold on
end
z = bindata(xData, yData,15,[0.8,6]);
hold on
plot(z(1,:),z(2,:),'MarkerSize',8,'Marker','square','LineWidth',2,'Color',[1 0 0]);
xlabel('LI (birth Length, um)');
ylabel('LF-LI (adding Length, um)');
set(gca,'FontSize',12)
xlim(gca,[0 6]);
ylim(gca,[0 6]);
hold on, line([0,6],[0,6],'Color','black');
fname='LI vs LF-LI';
saveas(gcf,[dirSave,'\',fname,'.fig']);
saveas(gcf,[dirSave,'\',fname,'.tif']);
hold off
close all


figure,
xData=gRDiv(:,5)*(scale^2);
yData=gRDiv(:,6)*(scale^2);
for iGenr=2:max(gRDiv(:,2))    
    scatter(xData(genrLogical(iGenr-1,:)),yData(genrLogical(iGenr-1,:)),50,cmap(iGenr,:),'filled','MarkerFaceAlpha',0.15);
    hold on
end
z = bindata(xData, yData,20);
hold on
plot(z(1,:),z(2,:),'MarkerSize',8,'Marker','square','LineWidth',2,'Color',[1 0 0]);
xlabel('LIa (birth Area, um^2)');
ylabel('LFa (division Area, um^2)');
set(gca,'FontSize',12)
xlim(gca,[0 3]);
ylim(gca,[0 6]);
hold on, line([0,3],[0,6],'Color','black');
fname='LIa vs LFa';
saveas(gcf,[dirSave,'\',fname,'.fig']);
saveas(gcf,[dirSave,'\',fname,'.tif']);
hold off
close all

figure,
xData=gRDiv(:,5)*(scale^2);
yData=(gRDiv(:,6)-gRDiv(:,5))*(scale^2);
for iGenr=2:max(gRDiv(:,2))    
    scatter(xData(genrLogical(iGenr-1,:)),yData(genrLogical(iGenr-1,:)),50,cmap(iGenr,:),'filled','MarkerFaceAlpha',0.15);
    hold on
end
z = bindata(xData, yData,20);
hold on
plot(z(1,:),z(2,:),'MarkerSize',8,'Marker','square','LineWidth',2,'Color',[1 0 0]);
xlabel('LIa (birth Area, um^2)');
ylabel('LFa-LIa (adding Area, um^2)');
set(gca,'FontSize',12)
xlim(gca,[0 3]);
ylim(gca,[0 3]);
hold on, line([0,3],[0,3],'Color','black');
fname='LIa vs LFa-LIa';
saveas(gcf,[dirSave,'\',fname,'.fig']);
saveas(gcf,[dirSave,'\',fname,'.tif']);
hold off
close all

figure,
xData=gRDiv(:,8);
yData=log(gRDiv(:,6)./gRDiv(:,5))./gRDiv(:,7);
filtLogical=yData>0&yData<150;
for iGenr=2:max(gRDiv(:,2))    
    scatter(xData(genrLogical(iGenr-1,:)&filtLogical'),yData(genrLogical(iGenr-1,:)&filtLogical'),50,cmap(iGenr,:),'filled','MarkerFaceAlpha',0.15);
    hold on
end
z = bindata(xData(filtLogical'), yData(filtLogical'),15,[0,150]);
hold on
plot(z(1,:),z(2,:),'MarkerSize',8,'Marker','square','LineWidth',2,'Color',[1 0 0]);
xlabel('Div time (min)');
ylabel('log(LFa/LIa)/gR (min)');
set(gca,'FontSize',12)
xlim(gca,[0 200]);
ylim(gca,[0 200]);
hold on, line([0,200],[0,200],'Color','black');
fname='Div vs gR';
saveas(gcf,[dirSave,'\',fname,'.fig']);
saveas(gcf,[dirSave,'\',fname,'.tif']);
hold off
close all

end

%% bin data
function z = bindata(x, y,nbin, xRange)

x=x(:);
y=y(:);

if nargin<3 || isempty(nbin)
    nbin = 50;
end
if nargin<4 || isempty(xRange)
    xRange=[min(x),max(x)];
end
% nbin=25;
% x1=linspace(min(x),max(x),nbin);
x1=linspace(xRange(1),xRange(2),nbin);
% x1=linspace(0.002,0.05,nbin);

z=zeros(2,size(x1,2)-1);
for i=1:size(x1,2)-1
    [row,col] = find(x>=x1(i)&x<x1(i+1));
    z(1,i)=(x1(i)+x1(i+1))/2;
    if isempty(row)
        continue
    end
    temp=[];
    for j=1:numel(row)
        temp(j)=y(row(j),col(j));
    end
    z(2,i)=mean(temp);
end
z=z(:,z(2,:)~=0);
end



