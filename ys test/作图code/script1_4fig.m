% figure,out = scatplot_modifiedys(dataAll{1}(:,1),dataAll{1}(:,2));
% hold on,out = scatplot_modifiedys(dataAll{2}(:,1)+2000,dataAll{2}(:,2));
% hold on,out = scatplot_modifiedys(dataAll{4}(:,1)+4000,dataAll{2}(:,2));
% hold on,out = scatplot_modifiedys(dataAll{6}(:,1)+6000,dataAll{6}(:,2));
% hold on,out = scatplot_modifiedys(dataAll{8}(:,1)+8000,dataAll{8}(:,2));

hold on,out = scatplot_modifiedys(dataAll{3}(:,1)+2000,dataAll{3}(:,2));
hold on,out = scatplot_modifiedys(dataAll{5}(:,1)+4000,dataAll{5}(:,2));
hold on,out = scatplot_modifiedys(dataAll{7}(:,1)+6000,dataAll{7}(:,2));
hold on,out = scatplot_modifiedys(dataAll{9}(:,1)+8000,dataAll{5}(:,2));
ylim(gca,[1 1000])
set(gca,'YMinorTick','on','YScale','log');