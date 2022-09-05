ms=10; %marker size
alpha=0.08;% points Í¸Ã÷¶È

figure, 
subplot(2,1,1),

scatter(0.75+(1.25-0.75)*rand(1,1500),rgLight(2:2:3000,1),ms,[0.85,0.33,0.1],'filled','MarkerFaceAlpha',alpha);
hold on,
scatter(1.75+(2.25-1.75)*rand(1,1500),rgLight(2:2:3000,2),ms,[0.85,0.33,0.1],'filled','MarkerFaceAlpha',alpha);
hold on,
scatter(2.75+(3.25-2.75)*rand(1,1500),rgLight(2:2:3000,3),ms,[0.85,0.33,0.1],'filled','MarkerFaceAlpha',alpha);
hold on,
scatter(3.75+(4.25-3.75)*rand(1,1500),rgLight(2:2:3000,4),ms,[0.85,0.33,0.1],'filled','MarkerFaceAlpha',alpha);
hold on,
scatter(4.75+(5.25-4.75)*rand(1,1500),rgLight(2:2:3000,5),ms,[0.85,0.33,0.1],'filled','MarkerFaceAlpha',alpha);
hold on,
scatter(5.75+(6.25-5.75)*rand(1,1500),rgLight(2:2:3000,6),ms,[0.85,0.33,0.1],'filled','MarkerFaceAlpha',alpha);
% hold on
% scatter(6.75+(7.25-6.75)*rand(1,1500),rgLight(2:2:3000,13),ms,[0.85,0.33,0.1],'filled','MarkerFaceAlpha',alpha);
hold on
boxplot(rgLight(2:2:3000,:),'Colors',[0.85,0.33,0.1],'symbol','r','Whisker',1);
ylim([0 0.35])