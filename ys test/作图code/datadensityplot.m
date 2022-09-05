function [scdata] = datadensityplot(x,y,colorindex)
%在矩阵方格里统计个数，计算数目的密度,两列数据 第一列x值，第二列y值；
% X must be a matrix with two columns.
% x,y 取对数，更关注于在对数坐标下的变化幅度,后期作图变换回原xy，x=10.^x;y=10.^y;
% x=log10(x);
% y=log10(y);

X(:,1)=x;
X(:,2)=y;

% colorindex=1;% 1代表蓝色 2 代表橙黄色

ms=25;%marker size

n=100; %默认取100*100格子计算，可以自定义

[N,C]=hist3(X,[n,n]);%二维平面格子里找到数据点的数目
N(N==0)=NaN;
k=1;
scdata=[];
for ix=1:n
    for iy=1:n
        if ~isnan(N(ix,iy))
            scdata(k,1)=C{1}(ix);
            scdata(k,2)=C{2}(iy);
            scdata(k,3)=N(ix,iy);
            k=k+1;
        end
    end
end

% [xi,yi]=meshgrid(C{1},C{2}); 
% zi = griddata(x,y,N,xi,yi);

scx=scdata(:,1);
scy=scdata(:,2);
c=scdata(:,3);


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
%color分成线性等分成64份
mymap(:,1)=linspace(color2(1),color1(1),I);
mymap(:,2)=linspace(color2(2),color1(2),I);
mymap(:,3)=linspace(color2(3),color1(3),I);
map = colormap(mymap);
ind = fix((c-min(c))/(max(c)-min(c))*(size(map,1)-1))+1;
h = [];
h=scatter(scx,scy,ms,map(ind,:),'filled');


end

