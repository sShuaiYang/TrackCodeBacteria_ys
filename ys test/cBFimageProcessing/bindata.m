%% bindata 作图

function z = bindata(x, y,nbin, xRange)
% Shuai Yang
% 2022/01/11
% 可以用disdiscretize 函数更方便 可以自定义边界和bin数目
% [Y,E] = discretize(sf,N);
%  for i = 1:N;TF = Y==i;b(i)= mean(mS(TF));c(i) = mean(E(i:i+1));end
% Shuai Yang 2022/5/11

x=x(:);
y=y(:);

if nargin<3 || isempty(nbin)
    nbin = 50;
end
if nargin<4 || isempty(xRange)
    xRange=[min(x),max(x)];
end
% nbin=25;
% x1=linspace(min(x),max(x),nbin);
x1 = linspace(xRange(1),xRange(2),nbin);
% x1=linspace(0.002,0.05,nbin);

z=zeros(2,size(x1,2)-1);
for i=1:size(x1,2)-1
    [row,col] = find(x>=x1(i)&x<x1(i+1));
    z(1,i)=(x1(i)+x1(i+1))/2;
    if isempty(row)
        continue
    end
    temp=[];
    for j=1:numel(row)
        temp(j)=y(row(j),col(j));
    end
    z(2,i)=mean(temp);
end
z=z(:,z(2,:)~=0);
end

% figure，
% plot(z(1,:),z(2,:),'MarkerSize',12,'Marker','.','LineWidth',1,'LineStyle','--','Color',[1 0 0])