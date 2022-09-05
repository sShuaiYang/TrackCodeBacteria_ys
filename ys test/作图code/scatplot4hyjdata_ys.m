function out = scatplot4hyjdata_ys(x,y)

% 此程序只用于hyj数据分析，根据程序scatplot_ys修改而来
N=100;%网格的尺寸
n=5; %二维卷积滤波的矩阵
po=3; %plot类型
ms=8;% markersize
method='ci';
radius = sqrt((5000/800)^2 + (10000/800)^2);

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
        r=0.05;%百分比作为检测幅度，e.g.5% or 10%
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
%         area = pi*r^2;%直接统计数据for hyj data
%         dd = dd/area;%直接统计数据for hyj data
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
% % 以下为原始作图函数
function varargout = gsp(x,y,c,ms)
%Graphs scattered poits
map = colormap(jet(256));

% ind = fix((c-min(c))/(max(c)-min(c))*(size(map,1)-1))+1;
nummax=max(c);
ind = fix((c-min(c))/(nummax-min(c))*(size(map,1)-1))+1;
h = [];
%much more efficient than matlab's scatter plot
for k=1:size(map,1)
    if any(ind==k)
        h(end+1) = line('Xdata',x(ind==k),'Ydata',y(ind==k), ...
            'LineStyle','none','Color',map(k,:), ...
            'Marker','.','MarkerSize',ms);
    end
end
if nargout==1
    varargout{1} = h;
end
return
