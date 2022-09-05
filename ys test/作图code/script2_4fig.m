figure,
subplot(2,5,1),[scdata] = datadensityplot(dataAll{1}(:,1),dataAll{1}(:,2),1);
xlim(subplot(2,5,1),[0 2000]);
ylim(subplot(2,5,1),[0 300])
% set(subplot(2,5,1),'XMinorTick','on','XScale','log','YMinorTick','on','YScale','log');

subplot(2,5,2),[scdata] = datadensityplot(dataAll{2}(:,1),dataAll{2}(:,2),1);
xlim(subplot(2,5,2),[0 2000]);
ylim(subplot(2,5,2),[0 300])
% set(subplot(2,5,2),'XMinorTick','on','XScale','log','YMinorTick','on','YScale','log');

subplot(2,5,3),[scdata] = datadensityplot(dataAll{4}(:,1),dataAll{4}(:,2),1);
xlim(subplot(2,5,3),[0 2000]);
ylim(subplot(2,5,3),[0 300])
% set(subplot(2,5,3),'XMinorTick','on','XScale','log','YMinorTick','on','YScale','log');

subplot(2,5,4),[scdata] = datadensityplot(dataAll{6}(:,1),dataAll{6}(:,2),1);
xlim(subplot(2,5,4),[0 2000]);
ylim(subplot(2,5,4),[0 300])
% set(subplot(2,5,4),'XMinorTick','on','XScale','log','YMinorTick','on','YScale','log');

subplot(2,5,5),[scdata] = datadensityplot(dataAll{8}(:,1),dataAll{8}(:,2),1);
xlim(subplot(2,5,5),[0 2000]);
ylim(subplot(2,5,5),[0 300])
% set(subplot(2,5,5),'XMinorTick','on','XScale','log','YMinorTick','on','YScale','log');

% hold (subplot(2,5,2),'on');[scdata] = datadensityplot(dataAll{3}(:,1),dataAll{3}(:,2),2);
% hold (subplot(2,5,3),'on');[scdata] = datadensityplot(dataAll{5}(:,1),dataAll{5}(:,2),2);
% hold (subplot(2,5,4),'on');[scdata] = datadensityplot(dataAll{7}(:,1),dataAll{7}(:,2),2);
% hold (subplot(2,5,5),'on');[scdata] = datadensityplot(dataAll{9}(:,1),dataAll{9}(:,2),2);