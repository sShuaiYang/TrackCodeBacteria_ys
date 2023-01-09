function barcreatefigure20221013(ymatrix1,ymatrix2,labelNames)
%CREATEFIGURE(ymatrix1)
%  YMATRIX1:  bar 矩阵数据

%  由 MATLAB 于 16-Oct-2022 19:00:03 自动生成
color1 = [0 0.45 0.75];
% color1 = [0.47,0.67,0.19];
color2 = [0.85 0.33 0.10];

% 创建 figure
figure1 = figure('Position',[565,556,1288,710],'Units','pixels');

tiledlayout(2,1,'TileSpacing','tight');
% 创建 axes
ax = nexttile;
hold(ax,'on');

% 使用 bar 的矩阵输入创建多行
bar1 = bar(ymatrix1,'EdgeColor','none','Parent',ax);
barFaceColor = color1;

set(bar1(1),'DisplayName','CT0100','FaceAlpha',1.0,...
    'FaceColor',barFaceColor);
set(bar1(2),'DisplayName','CT0300','FaceAlpha',0.9,...
    'FaceColor',barFaceColor);
set(bar1(3),'DisplayName','CT0500','FaceAlpha',0.8,...
    'FaceColor',barFaceColor);
set(bar1(4),'DisplayName','CT0900','FaceAlpha',0.7,...
    'FaceColor',barFaceColor);
set(bar1(5),'DisplayName','CT1800','FaceAlpha',0.6,...
    'FaceColor',barFaceColor);
set(bar1(6),'DisplayName','CT2400','FaceAlpha',0.4,...
    'FaceColor',barFaceColor);
set(bar1(7),'DisplayName','Dark',...
    'FaceColor',[0.50 0.50 0.50]);
set(bar1(8),'DisplayName','Light',...
    'FaceColor',[0.64 0.078 0.18]);

% 创建 ylabel
ylabel('sfmSRatio (a.u.) ','FontName','微软雅黑');
% xTickNames = {'P87-P42','P91-P42','P91-P72Q','P92-P71','P92-P72','P92-P74'};
% xTickNames = {'P87-P71Q','P87-P71','P87-P72Q','P87-P72','P91-P71','P91-P74Q'};
xTickNames = labelNames;
hold(ax,'off');
% 设置其余坐标区属性
set(ax,'Tickdir','out');
set(ax,'FontSize',12,'XTick',[1:numel(xTickNames)],'XTickLabel',...
    xTickNames,'TickLabelInterpreter','none');
% 创建 legend
legend1 = legend(ax,'show');
% set(legend1,...
%     'Position',[0.81 0.46 0.077 0.32]);

hold(ax,'off');

%% CymSRatio 作图
% 创建 axes
ax = nexttile;
hold(ax,'on');

% 使用 bar 的矩阵输入创建多行
bar1 = bar(ymatrix2,'EdgeColor','none','Parent',ax);
barFaceColor = color2;

set(bar1(1),'DisplayName','CT0100','FaceAlpha',1.0,...
    'FaceColor',barFaceColor);
set(bar1(2),'DisplayName','CT0300','FaceAlpha',0.9,...
    'FaceColor',barFaceColor);
set(bar1(3),'DisplayName','CT0500','FaceAlpha',0.8,...
    'FaceColor',barFaceColor);
set(bar1(4),'DisplayName','CT0900','FaceAlpha',0.7,...
    'FaceColor',barFaceColor);
set(bar1(5),'DisplayName','CT1800','FaceAlpha',0.6,...
    'FaceColor',barFaceColor);
set(bar1(6),'DisplayName','CT2400','FaceAlpha',0.4,...
    'FaceColor',barFaceColor);
set(bar1(7),'DisplayName','Dark',...
    'FaceColor',[0.50 0.50 0.50]);
set(bar1(8),'DisplayName','Light',...
    'FaceColor',[0.64 0.078 0.18]);

% 创建 ylabel
ylabel('CymSRatio (a.u.) ','FontName','微软雅黑');
% xTickNames = {'P87-P42','P91-P42','P91-P72Q','P92-P71','P92-P72','P92-P74'};
% xTickNames = {'P87-P71Q','P87-P71','P87-P72Q','P87-P72','P91-P71','P91-P74Q'};
xTickNames = labelNames;
hold(ax,'off');
% 设置其余坐标区属性
set(ax,'Tickdir','out');
set(ax,'FontSize',12,'XTick',[1:numel(xTickNames)],'XTickLabel',...
    xTickNames,'TickLabelInterpreter','none');
% 创建 legend
legend1 = legend(ax,'show');
% set(legend1,...
%     'Position',[0.81 0.46 0.077 0.32]);

hold(ax,'off');
