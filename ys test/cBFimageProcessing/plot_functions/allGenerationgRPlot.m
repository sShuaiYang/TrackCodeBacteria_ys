function [growthDataAll,gRPos_class,gRDiv]=allGenerationgRPlot(growthDataAll,dirSave)
% dirAll='\\192.168.1.12\e\2019-12-23 PAO1_IP32_100x ys';
% dirSave=strcat(dirAll,'\growthAnalysis');
% mkdir(dirSave);


[growthDataAll] = allgR_generationPlot(growthDataAll,dirSave);
allgR_kinPlot(growthDataAll,dirSave);
[growthDataAll,gRPos_class] = gRAndDistanceAnalysis(growthDataAll,dirSave);
[growthDataAll,gRDiv] = gRAndDivTimeAnalysis(growthDataAll,dirSave);
growthDataAll{1}.gRPos_class=gRPos_class;
growthDataAll{1}.gRDiv=gRDiv;
end
%% plot funciton-所有视野细菌生长率按代数作图
function [growthDataAll] = allgR_generationPlot(growthDataAll,dirSave)
avergegR_genr = zeros (numel(growthDataAll{1}.gRGetAll),4);%每一代细菌的平均生长率统计
figure,
for i=1:numel(growthDataAll{1}.gRGetAll)
    templogic=growthDataAll{1}.gRGetAll{i}(2,:)>0;
    scatter(growthDataAll{1}.gRGetAll{i}(1,templogic),growthDataAll{1}.gRGetAll{i}(2,templogic),50,'filled','MarkerFaceAlpha',0.3)
    hold on 
    avergegR_genr(i,1) = i;%代数
    avergegR_genr(i,2) = mean(growthDataAll{1}.gRGetAll{i}(2,templogic));%每一代平均gR
    avergegR_genr(i,3) = std(growthDataAll{1}.gRGetAll{i}(2,templogic));%每一代gR的方差
    avergegR_genr(i,4) = avergegR_genr(i,3)/avergegR_genr(i,2);%每一代gR的方差/平均值
end
hold on
plot (avergegR_genr(:,1),avergegR_genr(:,2),'MarkerSize',8,'Marker','square','LineWidth',2,'Color',[1 0 0])
ylabel('gR (min-1)');
xlabel('Generation');
xlim(gca,[0.9 i]);
ylim(gca,[0.001 0.1]);
set(gca,'FontSize',12,'XTick',(1:i),'YScale','log');
saveas(gcf,[dirSave,'\','gR_generation','.fig']);
saveas(gcf,[dirSave,'\','gR_generation','.tif']);
hold off
close all
growthDataAll{1}.avergegR_genr=avergegR_genr;
end
%% plot funciton-子代细菌与母代细菌生长率关系
function allgR_kinPlot(growthDataAll,dirSave)
figure,
adjAll=[];%adjacent gR
for i=1:numel(growthDataAll{1}.gR_kinAll)
    templogic=growthDataAll{1}.gR_kinAll{i}(:,1)>0&growthDataAll{1}.gR_kinAll{i}(:,2)>0;
    scatter(growthDataAll{1}.gR_kinAll{i}(templogic,1),growthDataAll{1}.gR_kinAll{i}(templogic,2),50,'filled','MarkerFaceAlpha',0.5)
    hold on 
    adjAll=[adjAll;growthDataAll{1}.gR_kinAll{i}(templogic,:)];    
end
z = bindata(adjAll(:,1), adjAll(:,2),25,[0.002,0.05]);
hold on
plot(z(1,:),z(2,:),'MarkerSize',8,'Marker','square','LineWidth',2,'Color',[1 0 0])
ylabel('g(n) min-1');
xlabel('g(n+1)  min-1');
title(['correff=',num2str(corr(adjAll(:,1),adjAll(:,2)),'%.03f')]);
set(gca,'FontSize',12)
xlim(gca,[0.001 0.05]);
ylim(gca,[0.001 0.05]);
saveas(gcf,[dirSave,'\','gR_kin','.fig']);
saveas(gcf,[dirSave,'\','gR_kin','.tif']);
hold off
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

