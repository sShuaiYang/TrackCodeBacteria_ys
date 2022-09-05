function cellLengthAndIntensityPlot(scale,dataCollect,dirSave)
% scale = 0.065;
L = dataCollect.MajorAxisLength*scale;
R = dataCollect.MinorAxisLength* scale/2;
V = (4/3)*pi*(R.^3) + pi*(R.^2).*(L-2*R);% volume 
TF = L>=1&L<=10;%长度筛选 [1,10]um
L = L(TF);
V = V(TF);

intData = dataCollect.intmScarletI;
intData = intData(TF);
fname = 'cellSizeVsInt';
% f1 = figure('WindowState','maximized','Name',fname);
f1 = figure('Name',fname);
subplot1 = subplot(2,2,1,'Parent',f1);
z = bindata(L, intData,46,[1,10]);
plot(z(1,:),z(2,:),'MarkerSize',12,'Marker','.','LineWidth',1,'LineStyle','--','Color',[1 0 0]);
% errorbar(z(1,:),z(2,:),z(3,:));
set(subplot1,'FontSize',12);
xlim([2,8]);
% ylim([20,50])
axis square
box(subplot1,'off');
xlabel('cell Length (um)');
ylabel('RFP intensity (a.u.)');

subplot2 = subplot(2,2,2,'Parent',f1);
z = bindata(V, intData,100,[0,15]);
plot(z(1,:),z(2,:),'MarkerSize',12,'Marker','.','LineWidth',1,'LineStyle','--','Color',[1 0 0]);
% errorbar(z(1,:),z(2,:),z(3,:));
set(subplot2,'FontSize',12,'XScale','log');
xlim([0,15]);
% ylim([20,50])
axis square
box(subplot2,'off');
xlabel('cell Volume (um^3)');
ylabel('RFP intensity (a.u.)');

subplot3 = subplot(2,2,3,'Parent',f1);
histogram2(L,intData,'BinMethod','scott','DisplayStyle','tile');
set(subplot3,'FontSize',12);
ylim([0,500]);
axis square
box(subplot3,'off');
xlabel('cell Length (um)');
ylabel('RFP intensity (a.u.)');

subplot4 = subplot(2,2,4,'Parent',f1);
histogram2(V,intData,'BinMethod','scott','DisplayStyle','tile');
set(subplot4,'FontSize',12);
ylim([0,500]);
axis square
box(subplot4,'off');
xlabel('cell Volume (um^3)');
ylabel('cell Num (a.u.)');

saveas(f1,[dirSave,'\',fname,'.fig']);
saveas(f1,[dirSave,'\',fname,'.tif'])
close all
end

%% bindata 作图
function z = bindata(x, y, nbin, xRange)

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
x1=linspace(xRange(1),xRange(2),nbin);
% x1=linspace(0.002,0.05,nbin);

z=zeros(3,size(x1,2)-1);
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
    z(2,i) = mean(temp);
    z(3,i) = std(temp);
end
z=z(:,z(2,:)~=0);
end