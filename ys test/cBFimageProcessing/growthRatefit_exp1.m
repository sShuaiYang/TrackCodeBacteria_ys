%% growth rate fit
function [gR,fitresult,gof] = growthRatefit_exp1(t, dataPickI)
% Shuai Yang 2021/01/10
% 用 y = a*exp(b*x) 进行拟合生长率 u = b
% 时间的初始需要为0
% dL/L = ut;L = L0*e^(ut);
% ln(L) = ln(L0) + ut;
% (L2-L1)/L1 = e^(ut2-ut1)-1;
% e^u(t2-t1) = L2/L1;u = ln(L2/L1)/(t2-t1);
% doubling time DT = log(2)/u;
% dataPickI 为细菌长度或者细菌面积随时间的变化
% timeSeries 为时间序列 hr /min /s

% row number

t = double(t(:));
t = t-t(1);
dataPickI = double(dataPickI(:));

% % 对细菌的生长信息进行smooth
dataPickI = smooth(dataPickI,'lowess');
dataPickI = smooth(dataPickI,'rlowess');

% 系统函数 对数据进行整理
[xData, yData] = prepareCurveData(t, dataPickI );

% Set up fittype and options.
ft = fittype( 'exp1');
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Robust = 'Bisquare';

warning('on');
% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
gR = fitresult.b;

% % Plot fit with data.
% figure( 'Name', 'fit_Result' );
% h = plot( fitresult, xData, yData );
% legend( h, 'b vs. a', 'fit_Result', 'Location', 'NorthEast', 'Interpreter', 'none' );
% set(gca,'YScale','log');

end