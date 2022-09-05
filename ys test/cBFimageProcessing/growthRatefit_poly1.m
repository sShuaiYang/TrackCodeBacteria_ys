%% growth rate fit
function [gR,fitperiodata,gof] = growthRatefit_poly1( t, dataPickI )
% Shuai Yang 2021/01/10
% ln(L) = ln(L0) + ut;
% 用 y = p1*x+p2 进行拟合生长率 u = p1
% 初始t不一定需要为0 影响p1值
% dL/L = ut;L = L0*e^(ut);
% ln(L) = ln(L0) + ut;
% (L2-L1)/L1 = e^(ut2-ut1)-1;
% e^u(t2-t1) = L2/L1;u = ln(L2/L1)/(t2-t1);
% doubling time DT = log(2)/u;
% dataPickI 为细菌长度或者细菌面积随时间的变化
% timeSeries 为时间序列 hr /min /s
% 也可参考 growthRatefit_exp1
% row number
t = t(:);
dataPickI = dataPickI(:);
% 取对数后 smooth
dataPickI = log(dataPickI);
dataPickI = smooth(dataPickI,'lowess');
dataPickI = smooth(dataPickI,'rlowess');

% Set up fittype and options.
ft = fittype( 'poly1' );
opts = fitoptions( ft );
opts.Lower = [-Inf -Inf];
opts.Robust = 'Bisquare';
opts.Upper = [Inf Inf];
warning('off');
% Fit model to data.
[fitperiodata,gof] = fit(t,dataPickI,ft,opts);
gR = fitperiodata.p1;
end