% points=3000;
% figure,
% subplot(2,5,1),out = scatplot_modifiedys(dataAll{1}(1:points,1),dataAll{1}(1:points,2),'re');
% xlim(subplot(2,5,1),[20 3000]);
% ylim(subplot(2,5,1),[1 600])
% set(subplot(2,5,1),'XMinorTick','on','XScale','log','YMinorTick','on','YScale','log');
% 
% subplot(2,5,2),out = scatplot_modifiedys(dataAll{2}(1:points,1),dataAll{2}(1:points,2),'re');
% xlim(subplot(2,5,2),[20 3000]);
% ylim(subplot(2,5,2),[1 600])
% set(subplot(2,5,2),'XMinorTick','on','XScale','log','YMinorTick','on','YScale','log');
% 
% subplot(2,5,3),out = scatplot_modifiedys(dataAll{4}(1:points,1),dataAll{4}(1:points,2),'re');
% xlim(subplot(2,5,3),[20 3000]);
% ylim(subplot(2,5,3),[1 600])
% set(subplot(2,5,3),'XMinorTick','on','XScale','log','YMinorTick','on','YScale','log');
% 
% subplot(2,5,4),out = scatplot_modifiedys(dataAll{6}(1:points,1),dataAll{6}(1:points,2),'re');
% xlim(subplot(2,5,4),[20 3000]);
% ylim(subplot(2,5,4),[1 600])
% set(subplot(2,5,4),'XMinorTick','on','XScale','log','YMinorTick','on','YScale','log');
% 
% subplot(2,5,5),out = scatplot_modifiedys(dataAll{8}(1:points,1),dataAll{8}(1:points,2),'re');
% xlim(subplot(2,5,5),[20 3000]);
% ylim(subplot(2,5,5),[1 600])
% set(subplot(2,5,5),'XMinorTick','on','XScale','log','YMinorTick','on','YScale','log');
% % 
hold (subplot(2,5,2),'on');out = scatplot_modifiedys(dataAll{3}(1:points,1),dataAll{3}(1:points,2),'re');
hold (subplot(2,5,3),'on');out = scatplot_modifiedys(dataAll{5}(1:points,1),dataAll{5}(1:points,2),'re');
hold (subplot(2,5,4),'on');out = scatplot_modifiedys(dataAll{7}(1:points,1),dataAll{7}(1:points,2),'re');
hold (subplot(2,5,5),'on');out = scatplot_modifiedys(dataAll{9}(1:points,1),dataAll{9}(1:points,2),'re');