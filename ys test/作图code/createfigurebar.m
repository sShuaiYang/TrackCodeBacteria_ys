function createfigurebar(xvector1, yvector1, X1, Y1)
%CREATEFIGURE(xvector1, yvector1, X1, Y1)
%  XVECTOR1:  bar xvector
%  YVECTOR1:  bar yvector
%  X1:  plot x 数据的向量
%  Y1:  plot y 数据的向量

%  由 MATLAB 于 02-Jun-2023 16:53:21 自动生成

% 创建 figure
figure1 = figure;

% 创建 axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');
colororder([0 0.447 0.741]);

% 激活坐标区的 left 侧
yyaxis(axes1,'left');
% 创建 bar
bar(xvector1,yvector1,...
    'FaceColor',[0.552941176470588 0.827450980392157 0.780392156862745],...
    'EdgeColor','none',...
    'BarWidth',0.5);

% 创建 ylabel
ylabel('销售额（亿美元）');

% 设置其余坐标区属性
set(axes1,'YColor',[0 0.447 0.741]);
% 激活坐标区的 right 侧
yyaxis(axes1,'right');
% 创建 plot
plot(X1,Y1,...
    'MarkerFaceColor',[0.850980392156863 0.325490196078431 0.0980392156862745],...
    'Marker','o',...
    'LineWidth',1.5);

% 创建 ylabel
ylabel('同比增长率');

% 设置其余坐标区属性
set(axes1,'YColor',[0.85 0.325 0.098],'YTick',...
    [-0.3 -0.2 -0.1 0 0.1 0.2 0.3],'YTickLabel',...
    {'-30%','-20%','-10%','0','10%','20%','30%'});
% 创建 title
title('2013-2023年全球前道涂胶显影设备销售额及增长','Editing','on','FontSize',14);

hold(axes1,'off');
% 设置其余坐标区属性
set(axes1,'FontName','微软雅黑','FontSize',12,'LineWidth',1,'XTick',...
    [1 2 3 4 5 6 7 8 9 10 11],'XTickLabel',...
    {'2013','2014','2015','2016','2017','2018','2019','2020E','2021E','2022E','2023E'},...
    'YGrid','on');
