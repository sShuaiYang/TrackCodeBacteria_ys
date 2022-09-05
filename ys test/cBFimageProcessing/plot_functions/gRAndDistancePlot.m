function  gRAndDistancePlot(gRPos_class,dirSave)
%函数gRAndDistanceAnalysis获得gRPos_class后计算
%gR与细菌直接距离关系作图
% plot function
for i=1:numel(gRPos_class)
    if ~isfield (gRPos_class{i},'plotData')
        continue
    end
    if size(gRPos_class{i}.plotData,1)<=1000 %数据点数目不能太少 50个
        continue
    else % 如果数目够多，进行nbin plot作图
        if size(gRPos_class{i}.plotData,1)>10000
            nbin=500;
        else
            nbin=100;
        end
    end
    xData=gRPos_class{i}.plotData(:,4);
    yData=gRPos_class{i}.plotData(:,3)/mean(gRPos_class{i}.plotData(:,3));%normalized by <g(r1)*g(r2)>
    z=bindata(xData, yData, nbin);
    figure,
    plot(z(1,:),z(2,:),'MarkerSize',12,'Marker','.','LineWidth',1,'LineStyle','none','Color',[1 0 0])
    %     set(gca,'FontSize',12,'YScale','log');
    set(gca,'FontSize',12,'YScale','linear');
    xlabel('r (um)');
    ylabel('<g(r1).g(r2)>');
    title(['Generation',num2str(i)]);
    fname = ['gR',num2str(i),' vs r'];
    saveas(gcf,[dirSave,'\',fname,'.fig']);
    saveas(gcf,[dirSave,'\',fname,'.tif']);
end
close all

for i=1:numel(gRPos_class)
    if ~isfield (gRPos_class{i},'plotData')
        continue
    end
    distCut=50;%distance less than 50um
    distLogical=gRPos_class{i}.plotData(:,4)<distCut;
    if sum(distLogical)>50 %在小于50um的前提下判断数据点数目
        
        xData = gRPos_class{i}.plotData(distLogical,4);
        yData = gRPos_class{i}.plotData(distLogical,3)/mean(gRPos_class{i}.plotData(distLogical,3));%normalized by <g(r1)*g(r2)>
        
        nbin = 200;%distCut/nbin = 0.25； 相当于0.25um的距离间隔
        z = bindata(xData, yData, nbin);
%         figure, %对数坐标作图
%         plot(z(1,:),z(2,:),'MarkerSize',12,'Marker','.','LineWidth',1,'LineStyle','none','Color',[1 0 0])
%         set(gca,'FontSize',12,'YScale','log');
%         %         set(gca,'FontSize',12,'YScale','linear');
%         xlabel('r (um)');
%         ylabel('<g(r1).g(r2)>');
%         title(['Generation',num2str(i)]);
%         fname=['gR',num2str(i),' vs r_',num2str(distCut),'um','_bin',num2str(nbin)];
%         saveas(gcf,[dirSave,'\',fname,'.fig']);
%         saveas(gcf,[dirSave,'\',fname,'.tif']);
        figure,
        plot(z(1,:),z(2,:),'MarkerSize',12,'Marker','.','LineWidth',1,'LineStyle','--','Color',[1 0 0])
        xlabel('r (um)');
        ylabel('<g(r1).g(r2)>');
        title(['Generation',num2str(i)]);
        fname1=['gR',num2str(i),' vs r_',num2str(distCut),'um','_bin',num2str(nbin),'_linear'];
        set(gca,'FontSize',12,'YScale','linear');
%         set(gca,'FontSize',12,'YScale','log');
        saveas(gcf,[dirSave,'\',fname1,'.fig']);
        saveas(gcf,[dirSave,'\',fname1,'.tif']);
        
        nbin = 100;%distCut/nbin = 0.5； 相当于0.5um的距离间隔
        z = bindata(xData, yData, nbin);
        figure,
        plot(z(1,:),z(2,:),'MarkerSize',12,'Marker','.','LineWidth',1,'LineStyle','--','Color',[1 0 0])
        xlabel('r (um)');
        ylabel('<g(r1).g(r2)>');
        title(['Generation',num2str(i)]);
        fname1=['gR',num2str(i),' vs r_',num2str(distCut),'um','_bin',num2str(nbin),'_linear'];
        set(gca,'FontSize',12,'YScale','linear');
        saveas(gcf,[dirSave,'\',fname1,'.fig']);
        saveas(gcf,[dirSave,'\',fname1,'.tif']);
        
        nbin = 50;%distCut/nbin = 1； 相当于1um的距离间隔
        z = bindata(xData, yData, nbin);
        figure,
        plot(z(1,:),z(2,:),'MarkerSize',12,'Marker','.','LineWidth',1,'LineStyle','--','Color',[1 0 0])
        xlabel('r (um)');
        ylabel('<g(r1).g(r2)>');
        title(['Generation',num2str(i)]);
        fname1=['gR',num2str(i),' vs r_',num2str(distCut),'um','_bin',num2str(nbin),'_linear'];
        set(gca,'FontSize',12,'YScale','linear');
        saveas(gcf,[dirSave,'\',fname1,'.fig']);
        saveas(gcf,[dirSave,'\',fname1,'.tif']);
        
        nbin = 25;%distCut/nbin = 2； 相当于2um的距离间隔
        z = bindata(xData, yData, nbin);
        figure,
        plot(z(1,:),z(2,:),'MarkerSize',12,'Marker','.','LineWidth',1,'LineStyle','--','Color',[1 0 0])
        xlabel('r (um)');
        ylabel('<g(r1).g(r2)>');
        title(['Generation',num2str(i)]);
        fname1=['gR',num2str(i),' vs r_',num2str(distCut),'um','_bin',num2str(nbin),'_linear'];
        set(gca,'FontSize',12,'YScale','linear');
        saveas(gcf,[dirSave,'\',fname1,'.fig']);
        saveas(gcf,[dirSave,'\',fname1,'.tif']);
        
    end
    
end
close all
end
%% bindata 作图

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
x1 = linspace(xRange(1),xRange(2),nbin);
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

