% ms=25;
% figure,
% scatter(normrnd(400,100,1,numel(dataAll{1}(:,2))),dataAll{1}(:,2)./dataAll{1}(:,1),ms,[0,0.45,0.74],'filled','MarkerFaceAlpha',0.03);
% hold on
% scatter(normrnd(1200,100,1,numel(dataAll{2}(:,2))),dataAll{2}(:,2)./dataAll{2}(:,1),ms,[0,0.45,0.74],'filled','MarkerFaceAlpha',0.03);
% hold on
% scatter(normrnd(2000,100,1,numel(dataAll{4}(:,2))),dataAll{4}(:,2)./dataAll{4}(:,1),ms,[0,0.45,0.74],'filled','MarkerFaceAlpha',0.03);
% hold on
% scatter(normrnd(2800,100,1,numel(dataAll{6}(:,2))),dataAll{6}(:,2)./dataAll{6}(:,1),ms,[0,0.45,0.74],'filled','MarkerFaceAlpha',0.03);
% hold on
% scatter(normrnd(3600,100,1,numel(dataAll{8}(:,2))),dataAll{8}(:,2)./dataAll{8}(:,1),ms,[0,0.45,0.74],'filled','MarkerFaceAlpha',0.03);
% hold on
% scatter(normrnd(1200,100,1,numel(dataAll{3}(:,2))),dataAll{3}(:,2)./dataAll{3}(:,1),ms,[0.85,0.33,0.1],'filled','MarkerFaceAlpha',0.03);
% hold on
% scatter(normrnd(2000,100,1,numel(dataAll{5}(:,2))),dataAll{5}(:,2)./dataAll{5}(:,1),ms,[0.85,0.33,0.1],'filled','MarkerFaceAlpha',0.03);
% hold on
% scatter(normrnd(2800,100,1,numel(dataAll{7}(:,2))),dataAll{7}(:,2)./dataAll{7}(:,1),ms,[0.85,0.33,0.1],'filled','MarkerFaceAlpha',0.03);
% hold on
% scatter(normrnd(3600,100,1,numel(dataAll{9}(:,2))),dataAll{9}(:,2)./dataAll{9}(:,1),ms,[0.85,0.33,0.1],'filled','MarkerFaceAlpha',0.03);

% ms=35;
% points=3000;%取多少个点作图
% subplot(2,1,2)
% scatter(normrnd(400,100,1,points),dataAll{1}(1:2:points*2,2)./dataAll{1}(1:2:points*2,1),ms,[0,0.45,0.74],'filled','MarkerFaceAlpha',0.03);
% hold on
% scatter(normrnd(1200,100,1,points),dataAll{2}(1:points,2)./dataAll{2}(1:points,1),ms,[0,0.45,0.74],'filled','MarkerFaceAlpha',0.03);
% hold on
% scatter(normrnd(2000,100,1,points),dataAll{4}(1:points,2)./dataAll{4}(1:points,1),ms,[0,0.45,0.74],'filled','MarkerFaceAlpha',0.03);
% hold on
% scatter(normrnd(2800,100,1,points),dataAll{6}(1:points,2)./dataAll{6}(1:points,1),ms,[0,0.45,0.74],'filled','MarkerFaceAlpha',0.03);
% hold on
% scatter(normrnd(3600,100,1,points),dataAll{8}(1:points,2)./dataAll{8}(1:points,1),ms,[0,0.45,0.74],'filled','MarkerFaceAlpha',0.03);
% hold on
% scatter(normrnd(400,100,1,points),dataAll{1}(2:2:points*2,2)./dataAll{1}(2:2:points*2,1),ms,[0.85,0.33,0.1],'filled','MarkerFaceAlpha',0.03);
% hold on
% scatter(normrnd(1200,100,1,points),dataAll{3}(1:points,2)./dataAll{3}(1:points,1),ms,[0.85,0.33,0.1],'filled','MarkerFaceAlpha',0.03);
% hold on
% scatter(normrnd(2000,100,1,points),dataAll{5}(1:points,2)./dataAll{5}(1:points,1),ms,[0.85,0.33,0.1],'filled','MarkerFaceAlpha',0.03);
% hold on
% scatter(normrnd(2800,100,1,points),dataAll{7}(1:points,2)./dataAll{7}(1:points,1),ms,[0.85,0.33,0.1],'filled','MarkerFaceAlpha',0.03);
% hold on
% scatter(normrnd(3600,100,1,points),dataAll{9}(1:points,2)./dataAll{9}(1:points,1),ms,[0.85,0.33,0.1],'filled','MarkerFaceAlpha',0.03);

hold(subplot(2,1,1),'on'),plot([800,800],[0,0.4],'Color',[0.15,0.15,0.15],'LineStyle','--','LineWidth' ,0.2)
plot([1600,1600],[0,0.4],'Color',[0.15,0.15,0.15],'LineStyle','--','LineWidth' ,0.2)
plot([2400,2400],[0,0.4],'Color',[0.15,0.15,0.15],'LineStyle','--','LineWidth' ,0.2)
plot([3200,3200],[0,0.4],'Color',[0.15,0.15,0.15],'LineStyle','--','LineWidth' ,0.2)

hold(subplot(2,1,2),'on'),plot([800,800],[0.01,1],'Color',[0.15,0.15,0.15],'LineStyle','--','LineWidth' ,0.2)
plot([1600,1600],[0.01,1],'Color',[0.15,0.15,0.15],'LineStyle','--','LineWidth' ,0.2)
plot([2400,2400],[0.01,1],'Color',[0.15,0.15,0.15],'LineStyle','--','LineWidth' ,0.2)
plot([3200,3200],[0.01,1],'Color',[0.15,0.15,0.15],'LineStyle','--','LineWidth' ,0.2)
