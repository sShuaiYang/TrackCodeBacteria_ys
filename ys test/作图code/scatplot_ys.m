function out = scatplot_ys(x,y,method,radius)
% x,y 取对数，更关注于在对数坐标下的变化幅度,后期作图变换回原xy，x=10.^x;y=10.^y;

% x=log10(x);
% y=log10(y);
% 用法 out = scatplot_ys(x,y,'re');
% 默认参数 method 不用输入 默认改为‘circles’，radius 也不用输入 默认已定义
N=100;%网格的尺寸
n=5; %二维卷积滤波的矩阵
po=3; %plot类型
ms=10;% markersize
%新增density 的method 'rectangle'，半径根据当前数据点的5%幅度作为检测范围

%主程序主要根据函数scatplot得到
%%
%以下为使用说明
% Scatter plot with color indicating data density
%
% USAGE:
%   out = scatplot(x,y,method,radius,N,n,po,ms)
%   out = scatplot(x,y,dd)
%
% DESCRIPTION:
%   Draws a scatter plot with a colorscale
%   representing the data density computed
%   using three methods
%
% INPUT VARIABLES:
%   x,y - are the data points
%   method - is the method used to calculate data densities:
%       'circles' - uses circles with a determined area
%               centered at each data point
%       'squares' - uses squares with a determined area
%               centered at each data point
%       'voronoi' - uses voronoi cells to determin data densities
%               default method is 'voronoi'
%   radius - is the radius used for the circles or squares
%       used to calculate the data densities if
%       (Note: only used in methods 'circles' and 'squares'
%           default radius is sqrt((range(x)/30)^2 + (range(y)/30)^2)
%   N - is the size of the square mesh (N x N) used to
%       filter and calculate contours
%       default is 100
%   n - is the number of coeficients used in the 2-D
%       running mean filter
%       default is 5
%       (Note: if n is length(2), n(2) is tjhe number of
%       of times the filter is applied)
%   po - plot options:
%       0 - No plot
%       1 - plots only colored data points (filtered)
%       2 - plots colored data points and contours (filtered)
%       3 - plots only colored data points (unfiltered)
%       4 - plots colored data points and contours (unfiltered)
%           default is 1
%   ms - uses this marker size for filled circles
%       default is 4
%
% OUTPUT VARIABLE:
%   out - structure array that contains the following fields:
%       dd - unfiltered data densities at (x,y)
%       ddf - filtered data densities at (x,y)
%       radius - area used in 'circles' and 'squares'
%               methods to calculate densities
%       xi - x coordenates for zi matrix
%       yi - y coordenates for zi matrix
%       zi - unfiltered data densities at (xi,yi)
%       zif - filtered data densities at (xi,yi)
%       [c,h] = contour matrix C as described in
%           CONTOURC and a handle H to a contourgroup object
%       hs = scatter points handles
%
%Copy-Left, Alejandro Sanchez-Barba, 2005
%%
if nargin==0
    out=scatplotdemo;
    return
end
% if nargin<3 | isempty(method)
%     method = 'vo';
% end
%将默认的半径方式改为圆形'circles'
if nargin<3 || isempty(method)
    method = 'ci';
end
if isnumeric(method)
    gsp(x,y,method,2)
    return
else
    method = method(1:2);
end
if nargin<4 || isempty(radius)
    radius = sqrt((range(x)/400)^2 + (range(y)/400)^2);
    %     radius = sqrt((range(x)/50)^2 + (range(y)/50)^2);
end
%%
%以下修改 用默认参数
% if nargin<4 | isempty(n)
%     n = 5; %number of filter coefficients
% end
% if nargin<5 | isempty(radius)
%     radius = sqrt((range(x)/30)^2 + (range(y)/30)^2);
% end
% if nargin<6 | isempty(po)
%     po = 1; %plot option
% end
% if nargin<7 | isempty(ms)
%     ms = 4; %markersize
% end
% if nargin<8 | isempty(N)
%     N = 100; %length of grid
% end
%%
%Correct data if necessary
x = x(:);
y = y(:);
%Asuming x and y match
idat = isfinite(x);
x = x(idat);
y = y(idat);
% figure;
holdstate = ishold;
if holdstate==0
    cla
end
hold on
%--------- Caclulate data density ---------
dd = datadensity(x,y,method,radius);
%------------- Gridding -------------------
xi = repmat(linspace(min(x),max(x),N),N,1);
yi = repmat(linspace(min(y),max(y),N)',1,N);
% [xi,yi]=meshgrid(linspace(min(x),max(x),N),linspace(min(y),max(y),N));
zi = griddata(x,y,dd,xi,yi);
%----- Bidimensional running mean filter -----
zi(isnan(zi)) = 0;
coef = ones(n(1),1)/n(1);
zif = conv2(coef,coef,zi,'same');
if length(n)>1
    for k=1:n(2)
        zif = conv2(coef,coef,zif,'same');
    end
end
%-------- New Filtered data densities --------
ddf = griddata(xi,yi,zif,x,y);
%----------- Plotting --------------------
switch po
    case {1,2}
        if po==2
            [c,h] = contour(xi,yi,zif);
            out.c = c;
            out.h = h;
        end %if
        hs = gsp(x,y,ddf,ms);
        out.hs = hs;
        colorbar
    case {3,4}
        if po>3
            [c,h] = contour(xi,yi,zi);
            out.c = c;
        end %if
        hs = gsp(x,y,dd,ms);
        out.hs = hs;
        colorbar
end %switch
%------Relocate variables and place NaN's ----------
dd(idat) = dd;
dd(~idat) = NaN;
ddf(idat) = ddf;
ddf(~idat) = NaN;
%--------- Collect variables ----------------
out.dd = dd;
out.ddf = ddf;
out.radius = radius;
out.xi = xi;
out.yi = yi;
out.zi = zi;
out.zif = zif;
if ~holdstate
    hold off
end
return
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function out=scatplotdemo
po = 3;
method = 'circles';
N = [];
n = [];
ms = 5;
x = randn(2000,1);
y = randn(2000,1);
radius = sqrt((range(x)/50)^2 + (range(y)/50)^2);
out=scatplot_ys(x,y,method,radius);

return
%~~~~~~~~~~ Data Density ~~~~~~~~~~~~~~
function dd = datadensity(x,y,method,r)
%Computes the data density (points/area) of scattered points
%Striped Down version
%
% USAGE:
%   dd = datadensity(x,y,method,radius)
%
% INPUT:
%   (x,y) -  coordinates of points
%   method - either 'squares','circles', or 'voronoi'
%       default = 'voronoi'
%   radius - Equal to the circle radius or half the square width
Ld = length(x);
dd = zeros(Ld,1);
switch method %Calculate Data Density
    % ys add rectangle method 
    case 're'  %---- Using rectangle ----
        r=0.1;%百分比作为检测幅度，e.g.5% or 10%
        for k=1:Ld
            dd(k) = sum( x>((1-r)*x(k)) & x<((1+r)*x(k)) & y>((1-r)*y(k)) & y<((1+r)*y(k)) );
            %每个数据点搜寻的范围都在变化，可以视情况用面积进行归一
            area=(r*2*x(k))*(r*2*y(k));
            dd(k)=dd(k)/area;
        end %for
    case 'sq'  %---- Using squares ----
        for k=1:Ld
            dd(k) = sum( x>(x(k)-r) & x<(x(k)+r) & y>(y(k)-r) & y<(y(k)+r) );
        end %for
        area = (2*r)^2;
        dd = dd/area;
    case 'ci'
        for k=1:Ld
            dd(k) = sum( sqrt((x-x(k)).^2 + (y-y(k)).^2) < r );
        end
        area = pi*r^2;
        dd = dd/area;
    case 'vo'  %----- Using voronoi cells ------
        [v,c] = voronoin([x,y]);
        for k=1:length(c)
            %If at least one of the indices is 1,
            %then it is an open region, its area
            %is infinity and the data density is 0
            if all(c{k}>1)
                a = polyarea(v(c{k},1),v(c{k},2));
                dd(k) = 1/a;
            end %if
        end %for
end %switch
return
%~~~~~~~~~~ Graf Scatter Plot ~~~~~~~~~~~
%%
%以下为原始作图函数
% function varargout = gsp(x,y,c,ms)
% %Graphs scattered poits
% map = colormap;
% ind = fix((c-min(c))/(max(c)-min(c))*(size(map,1)-1))+1;
% h = [];
% %much more efficient than matlab's scatter plot
% for k=1:size(map,1)
%     if any(ind==k)
%         h(end+1) = line('Xdata',x(ind==k),'Ydata',y(ind==k), ...
%             'LineStyle','none','Color',map(k,:), ...
%             'Marker','.','MarkerSize',ms);
%     end
% end
% if nargout==1
%     varargout{1} = h;
% end
% return
%%
%稍作修改
% function varargout = gsp(x,y,c,ms)
% %Graphs scattered poits
% % 如果xy初始化时取了对数，此时需要变换回原来的数值
% % x=10.^x;
% % y=10.^y;
% colorindex=2;
% if colorindex==1
%     %蓝色
%     color1=[0,0.45,0.74];
%     color2=[0.95,0.98,1];
% end
% if colorindex==2
%     %橙黄色
%     color1=[0.85,0.33,0.1];
%     color2=[1,0.92,0.89];
% end
% I=64;
% %color分成线性等分成64份
% mymap(:,1)=linspace(color2(1),color1(1),I);
% mymap(:,2)=linspace(color2(2),color1(2),I);
% mymap(:,3)=linspace(color2(3),color1(3),I);
% 
% % %color分成对数等分成64份
% % mymap(:,1)=logspace(log10(color2(1)),log10(color1(1)),I);
% % mymap(:,2)=logspace(log10(color2(2)),log10(color1(2)),I);
% % mymap(:,3)=logspace(log10(color2(3)),log10(color1(3)),I);
% 
% map = colormap(mymap);
% ind = fix((c-min(c))/(max(c)-min(c))*(size(map,1)-1))+1;
% h = [];
% %much more efficient than matlab's scatter plot
% for k=1:size(map,1)
%     if any(ind==k)
%         h(end+1) = line('Xdata',x(ind==k),'Ydata',y(ind==k), ...
%             'LineStyle','none','Color',map(k,:), ...
%             'Marker','o','MarkerSize',ms,...
%             'MarkerFaceColor',map(k,:));
%     end
% end
% if nargout==1
%     varargout{1} = h;
% end
% return
%%
% 修改用matlab scatter 函数作图
function varargout = gsp(x,y,c,ms)
%Graphs scattered poits
% delete NaN data of c.用‘re’计算时，如果用x或y等于0时，面积会等于0，约化后
%会出现分母等于0
x=x(~isnan(c));
y=y(~isnan(c));
c=c(~isnan(c));
%define the color
colorindex=2;
if colorindex==1
    %蓝色
    color1=[0,0.45,0.74];
    color2=[0.95,0.98,1];
end
if colorindex==2
    %橙黄色
    color1=[0.85,0.33,0.1];
    color2=[1,0.92,0.89];
end
I=64;
% color分成线性等分成64份
mymap(:,1)=linspace(color2(1),color1(1),I);
mymap(:,2)=linspace(color2(2),color1(2),I);
mymap(:,3)=linspace(color2(3),color1(3),I);

% %color分成对数等分成64份
% mymap(:,1)=logspace(log10(color2(1)),log10(color1(1)),I);
% mymap(:,2)=logspace(log10(color2(2)),log10(color1(2)),I);
% mymap(:,3)=logspace(log10(color2(3)),log10(color1(3)),I);

map = colormap(mymap);
map = colormap(jet(256));
ind = fix((c-min(c))/(max(c)-min(c))*(size(map,1)-1))+1;
h = [];
h=scatter(x,y,ms,map(ind,:),'filled','MarkerFaceAlpha',0.5);
if nargout==1
    varargout{1} = h;
end
return
%%
