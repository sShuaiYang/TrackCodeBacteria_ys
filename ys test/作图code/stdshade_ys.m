function stdshade_ys(amatrix,alpha,acolor,F,smth)
%适用于数据点比较密时，将误差棒(std或者sem)作成阴影的,行代表实验重复次数
% - amatrix define the data， coloumn define datapoint from different
% timepoint, row define repeated times.
% - acolor defines the used color (default is red)
% - F assignes the used x axis (default is steps of 1).start from 0.
% - alpha defines transparency of the shading (default is no shading and black mean line)
% - smth defines the smoothing factor (default is no smooth)

if nargin==0
    stdshadedemo
    return
end

if exist('acolor','var')==0 || isempty(acolor)
    acolor='k';
end

if exist('F','var')==0 || isempty(F)
    F=0:size(amatrix,1)-1;
else
    F=(0:size(amatrix,1)-1)*F;
end
F=F';

if exist('smth','var')
    if isempty(smth)
        smth=1;
    end
else
    smth=1;
end

amean=mean(amatrix,2);
amean=smooth(amean,smth);
astd=std(amatrix,0,2); % to get std shading
% astd=std(amatrix,0,2)/sqrt(size(amatrix,2)); % to get sem shading

if exist('alpha','var')==0 || isempty(alpha)
    alpha=0.2;
    fill([F; flipud(F)],[amean+astd; flipud(amean-astd)], acolor, 'FaceAlpha', alpha,'linestyle','none');
    acolor='r';
else
    fill([F; flipud(F)],[amean+astd; flipud(amean-astd)], acolor, 'FaceAlpha', alpha,'linestyle','none');
end



if ishold==0
    check=true;
else
    check=false;
end

hold on;
plot(F,amean,'color',acolor,'linewidth',1.5); %% change color or linewidth to adjust mean line

if check
    hold off;
end

end

function stdshadedemo
x = 0:10:500;
y_true =  30*sind(x) + x/10;
sigma = 3;
y_measured(:,1) = y_true + sigma*rand(size(x)).*x/100;
y_measured(:,2) = y_true -sigma*rand(size(x)).*x/100;
y_measured(:,3) = y_true +0.5*sigma*rand(size(x)).*x/100;
y_measured(:,4) = y_true -0.5*sigma*rand(size(x)).*x/100;
stdshade_ys(y_measured,[],[],10,10);
end
