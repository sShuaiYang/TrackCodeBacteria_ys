function  deltgRAndDistancePlot(gRPos_class,dirSave)
%函数gRAndDistanceAnalysis获得gRPos_class后计算
%gR与细菌直接距离关系作图
% plot function

for i=1:numel(gRPos_class)
    if ~isfield (gRPos_class{i},'plotData')
        continue
    end
    distCut=50;%distance less than 50um
    distLogical=gRPos_class{i}.plotData(:,4)<distCut;
    if sum(distLogical)>50 %在小于50um的前提下判断数据点数目
        meangR = mean(gRPos_class{i}.gRPos(:,3));
        
        gR1_norm = gRPos_class{i}.plotData(distLogical,1) - meangR;
        gR2_norm = gRPos_class{i}.plotData(distLogical,2) - meangR;
        yData = gR1_norm .* gR2_norm;
        xData = gRPos_class{i}.plotData(distLogical,4);

        
        nbin = 100;%distCut/nbin = 0.5； 相当于0.5um的距离间隔
        z = bindata(xData, yData, nbin);
        figure,
        plot(z(1,:),z(2,:),'MarkerSize',12,'Marker','.','LineWidth',1,'LineStyle','--','Color',[1 0 0])
        xlabel('r (um)');
%         ylabel(strcat('<','\Delta','g(r1).','\Delta','g(r2)>'));
        ylabel('<\Deltag(r1).\Deltag(r2)>');
        title(['Generation',num2str(i)]);
        fname1 = ['deltgR',num2str(i),' vs r'];
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

