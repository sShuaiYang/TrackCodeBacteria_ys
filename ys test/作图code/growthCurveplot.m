function [sta_data] = growthCurveplot(ODdata)
%OD散点平滑连接起来
%   ODdata 默认按行排列，不同行数据求误差
x=0:2:10;
y_ori=ODdata()';%转换成列向量
y=mean(y_ori(:,1:end),2);
err=std(y_ori(:,1:end),0,2);
errorbar(x,y,err,'LineStyle','none','LineWidth',1.5);

sta_data(:,1)=y;
sta_data(:,2)=err;

% 利用 spline或 pchip 生成插值，进行平滑连接
hold on
xq1 = x(1):.01:x(end);
p = pchip(x,y,xq1);
s = spline(x,y,xq1);
plot(xq1,p,'-','LineWidth',1.5)
% plot(x,y,'o',xq1,p,'-',xq1,s,'-.')
% plot(x,y,'o',xq1,p,'-',xq1,s,'-.')
end

