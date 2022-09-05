% %boxplot script
% indexDL = [    1     1     0     1     0     1     0     0     0     1     0   1;
%      1     0     1     0     1     0     1     0     1     0     1  0];
% 
% indexDL=logical(indexDL);
% figure, scatter(0.75+(1.25-0.75)*rand(1,1500),rg(1:3000,1),50,[0,0.45,0.74],'filled','MarkerFaceAlpha',0.01);
% 
% figure,subplot(2,1,1),boxplot(rg(1:2:3000,indexDL(1,:)),'Colors',[0,0.45,0.74],'symbol','r','Whisker',1);
% ylim([0 0.35])
% subplot(2,1,2),boxplot(rg(2:2:3000,indexDL(2,:)),'Colors',[0.85,0.33,0.1],'symbol','r','Whisker',1)
% ylim([0 0.35])
% 
ms=10; %marker size
alpha=0.08;% points 透明度

figure, 
subplot(2,1,1),
scatter(0.75+(1.25-0.75)*rand(1,1500),rg(1:2:3000,1),ms,[0,0.45,0.74],'filled','MarkerFaceAlpha',alpha);
hold on,
scatter(1.75+(2.25-1.75)*rand(1,1500),rg(1:2:3000,2),ms,[0,0.45,0.74],'filled','MarkerFaceAlpha',alpha);
hold on,
scatter(2.75+(3.25-2.75)*rand(1,1500),rg(1:2:3000,4),ms,[0,0.45,0.74],'filled','MarkerFaceAlpha',alpha);
hold on,
scatter(3.75+(4.25-3.75)*rand(1,1500),rg(1:2:3000,6),ms,[0,0.45,0.74],'filled','MarkerFaceAlpha',alpha);
hold on,
scatter(4.75+(5.25-4.75)*rand(1,1500),rg(1:2:3000,8),ms,[0,0.45,0.74],'filled','MarkerFaceAlpha',alpha);
hold on,
scatter(5.75+(6.25-5.75)*rand(1,1500),rg(1:2:3000,10),ms,[0,0.45,0.74],'filled','MarkerFaceAlpha',alpha);
hold on
% scatter(6.75+(7.25-6.75)*rand(1,1500),rg(1:2:3000,12),ms,[0,0.45,0.74],'filled','MarkerFaceAlpha',alpha);
% hold on 
boxplot(rg(1:2:3000,indexDL(1,:)),'Colors',[0,0.45,0.74],'symbol','r','Whisker',1);
ylim([0 0.35])

% subplot(2,1,2),
hold on
scatter(0.75+(1.25-0.75)*rand(1,1500),rg(2:2:3000,1),ms,[0.85,0.33,0.1],'filled','MarkerFaceAlpha',alpha);
hold on,
scatter(1.75+(2.25-1.75)*rand(1,1500),rg(2:2:3000,3),ms,[0.85,0.33,0.1],'filled','MarkerFaceAlpha',alpha);
hold on,
scatter(2.75+(3.25-2.75)*rand(1,1500),rg(2:2:3000,5),ms,[0.85,0.33,0.1],'filled','MarkerFaceAlpha',alpha);
hold on,
scatter(3.75+(4.25-3.75)*rand(1,1500),rg(2:2:3000,7),ms,[0.85,0.33,0.1],'filled','MarkerFaceAlpha',alpha);
hold on,
scatter(4.75+(5.25-4.75)*rand(1,1500),rg(2:2:3000,9),ms,[0.85,0.33,0.1],'filled','MarkerFaceAlpha',alpha);
hold on,
scatter(5.75+(6.25-5.75)*rand(1,1500),rg(2:2:3000,11),ms,[0.85,0.33,0.1],'filled','MarkerFaceAlpha',alpha);
% hold on
% scatter(6.75+(7.25-6.75)*rand(1,1500),rg(2:2:3000,13),ms,[0.85,0.33,0.1],'filled','MarkerFaceAlpha',alpha);
hold on
boxplot(rg(2:2:3000,indexDL(2,:)),'Colors',[0.85,0.33,0.1],'symbol','r','Whisker',1);
ylim([0 0.35])

%% 一条直线点图
% figure,
% subplot(2,1,1),
% scatter(1*ones(1,1500),rg(1:2:3000,1),ms,[0,0.45,0.74],'filled','MarkerFaceAlpha',alpha);
% hold on,
% scatter(2*ones(1,1500),rg(1:2:3000,2),ms,[0,0.45,0.74],'filled','MarkerFaceAlpha',alpha);
% hold on,
% scatter(3*ones(1,1500),rg(1:2:3000,4),ms,[0,0.45,0.74],'filled','MarkerFaceAlpha',alpha);
% hold on,
% scatter(4*ones(1,1500),rg(1:2:3000,6),ms,[0,0.45,0.74],'filled','MarkerFaceAlpha',alpha);
% hold on,
% scatter(5*ones(1,1500),rg(1:2:3000,8),ms,[0,0.45,0.74],'filled','MarkerFaceAlpha',alpha);
% hold on,
% scatter(6*ones(1,1500),rg(1:2:3000,10),ms,[0,0.45,0.74],'filled','MarkerFaceAlpha',alpha);
% hold on
% boxplot(rg(1:2:3000,indexDL(1,:)),'Colors',[0,0.45,0.74],'symbol','r','Whisker',1);
% ylim([0 0.35])
% 
% % subplot(2,1,2),
% hold on
% scatter(1*ones(1,1500),rg(2:2:3000,1),ms,[0.85,0.33,0.1],'filled','MarkerFaceAlpha',alpha);
% hold on,
% scatter(2*ones(1,1500),rg(2:2:3000,3),ms,[0.85,0.33,0.1],'filled','MarkerFaceAlpha',alpha);
% hold on,
% scatter(3*ones(1,1500),rg(2:2:3000,5),ms,[0.85,0.33,0.1],'filled','MarkerFaceAlpha',alpha);
% hold on,
% scatter(4*ones(1,1500),rg(2:2:3000,7),ms,[0.85,0.33,0.1],'filled','MarkerFaceAlpha',alpha);
% hold on,
% scatter(5*ones(1,1500),rg(2:2:3000,9),ms,[0.85,0.33,0.1],'filled','MarkerFaceAlpha',alpha);
% hold on,
% scatter(6*ones(1,1500),rg(2:2:3000,11),ms,[0.85,0.33,0.1],'filled','MarkerFaceAlpha',alpha);
% hold on
% boxplot(rg(2:2:3000,indexDL(2,:)),'Colors',[0.85,0.33,0.1],'symbol','r','Whisker',1);
% ylim([0 0.35])
% % scatter(0.75+(1.25-0.75)*rand(1,1500),rg(1:2:3000,1),50,[0,0.45,0.74],'filled','MarkerFaceAlpha',0.01);
% % scatter(0.75+(1.25-0.75)*rand(1,1500),rg(1:2:3000,1),50,[0,0.45,0.74],'filled','MarkerFaceAlpha',0.01);