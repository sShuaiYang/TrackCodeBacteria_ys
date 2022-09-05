ms=10; %marker size
alpha=0.08;% points 透明度

k = 1;
for i=1:11
    
    rg(:,k) = dataAll{i}(1:3000,2)./dataAll{i}(1:3000,1);
    k = k+1;
end

for i = 1:numel(dataAll)
    a = dataAll{i};
    TF = a(:,1)>200;
    a = a(TF,:);
    dataAll{i} = a;
end

for i = 1:numel(dataAll)
    a = dataAll{i};
    TF = a(:,2)<20;
    a = a(TF,:);
    dataAll{i} = a;
end

% k = 1;
% for i=1:11
%     if i~= 3
%         rg(:,k) = dataAll{i}(1:3000,2)./dataAll{i}(1:3000,1);
%         k = k+1;
%     end
% end
% indexDL = [    1     1     0     1     0     1     0     0     0     1  ;
%      1     0     1     0     1     0     1     0     1     0     ];
% indexDL = [1,2,4,6,8;1,3,5,7,9];
% indexDL = [1,3,5,7,9;2,4,6,8,10];
indexDL = [1,2,4,6,8,10;1,3,5,7,9,11;];
% color1=[0,0.45,0.74];%蓝色
% color2=[0.85,0.33,0.1];%橙色

% color2=[0.93,0.69,0.13];
% color2=[126,157,164]/255;
% color2=[44,123,182]/255;
% color1=[18,13,10]/255;
% color1=[20,76,112]/255;
% color2=[119,127,54]/255;

color1=[89,87,87]/255;
color2=[72,124,177]/255;
%light data
figure, 
subplot(2,1,1),
scatter(0.75+(1.25-0.75)*rand(1,1500),rg(1:2:3000,indexDL(1,1)),ms,color1,'filled','MarkerFaceAlpha',alpha);
hold on,
scatter(1.75+(2.25-1.75)*rand(1,1500),rg(1:2:3000,indexDL(1,2)),ms,color1,'filled','MarkerFaceAlpha',alpha);
hold on,
scatter(2.75+(3.25-2.75)*rand(1,1500),rg(1:2:3000,indexDL(1,3)),ms,color1,'filled','MarkerFaceAlpha',alpha);
hold on,
scatter(3.75+(4.25-3.75)*rand(1,1500),rg(1:2:3000,indexDL(1,4)),ms,color1,'filled','MarkerFaceAlpha',alpha);
hold on,
scatter(4.75+(5.25-4.75)*rand(1,1500),rg(1:2:3000,indexDL(1,5)),ms,color1,'filled','MarkerFaceAlpha',alpha);
hold on,
scatter(5.75+(6.25-5.75)*rand(1,1500),rg(1:2:3000,indexDL(1,6)),ms,color1,'filled','MarkerFaceAlpha',alpha);
hold on
% scatter(6.75+(7.25-6.75)*rand(1,1500),rg(1:2:3000,12),ms,[0,0.45,0.74],'filled','MarkerFaceAlpha',alpha);
% hold on 
boxplot(rg(1:2:3000,indexDL(1,:)),'Colors',color1,'symbol','r','Whisker',1);
set(gca,'YScale','log')
% ylim([0.05 1])
ylim([0.005 0.1])
box(gca,'off');
%dark data
alpha=0.08;% points 透明度
subplot(2,1,2),
% hold on
scatter(0.75+(1.25-0.75)*rand(1,1500),rg(2:2:3000,indexDL(2,1)),ms,color2,'filled','MarkerFaceAlpha',alpha);
hold on,
scatter(1.75+(2.25-1.75)*rand(1,1500),rg(2:2:3000,indexDL(2,2)),ms,color2,'filled','MarkerFaceAlpha',alpha);
hold on,
scatter(2.75+(3.25-2.75)*rand(1,1500),rg(2:2:3000,indexDL(2,3)),ms,color2,'filled','MarkerFaceAlpha',alpha);
hold on,
scatter(3.75+(4.25-3.75)*rand(1,1500),rg(2:2:3000,indexDL(2,4)),ms,color2,'filled','MarkerFaceAlpha',alpha);
hold on,
scatter(4.75+(5.25-4.75)*rand(1,1500),rg(2:2:3000,indexDL(2,5)),ms,color2,'filled','MarkerFaceAlpha',alpha);
hold on,
scatter(5.75+(6.25-5.75)*rand(1,1500),rg(2:2:3000,indexDL(2,6)),ms,color2,'filled','MarkerFaceAlpha',alpha);
% hold on
% scatter(6.75+(7.25-6.75)*rand(1,1500),rg(2:2:3000,13),ms,[0.85,0.33,0.1],'filled','MarkerFaceAlpha',alpha);
hold on
boxplot(rg(2:2:3000,indexDL(2,:)),'Colors',color2,'symbol','r','Whisker',1);
% ylim([0 0.6])
set(gca,'YScale','log')
% ylim([0.05 1])
ylim([0.005 0.1])
box(gca,'off');