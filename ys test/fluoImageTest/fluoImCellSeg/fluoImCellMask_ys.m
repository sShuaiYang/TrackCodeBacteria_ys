function [fluoImMask] = fluoImCellMask_ys(fluoImage)
%% 利用edge图像进行边缘检测
% watershed transform进行segementaion
% Shuai Yang 2020.10.11
% fluoImSeg_edge_Test.mlx;segCellTest_Watershedtransform.mlx;
% segCellTest_splitKinkyCells.mlx;segCellTest_splitFatCells.mlx

%%
% get mask
% &&& not suitable for  weak fluorescence &&&%
% &&& not suitable for  strong fluo loci in cells &&& %
fluoImMask = false(size(fluoImage));
I0 = fluoImage;
[~,BG] = substractBackGround(I0);%获得背景数值 backGround signal
% I2 = imgaussfilt(I0,2);%100x maybe sigma = 2; 60x sigma = [0.5,1];
I2 = imgaussfilt(I0,1);

if max(double(I2(:))) > (BG + 80) %weak fluo cause err

    x = double(I2);
    s = sort(x(:));
    % small = s(25);
    % case for imdrift correction pixelvalue = 0
    small = BG;
    %     big= s(end-25);
    signal = s-(small+40);% 40 could change
    signal(signal<=0) = [];
    if isempty(signal)
        return
    end
    signal = signal+(small+40);
    %     big = signal(round(length(signal)*0.93));% 0.93 could change
    big = ImBigSignalGet(signal);
    rescaled = (x - small)/(big - small);
    rescaled(rescaled<0) = 0;
    I3 = uint16(10000*rescaled);


    % use a threshold to find where the cells are
    I3 = I3 - 2000;% default 100x 2000
    bw = edge(I3,'log',0);

    bw1 = imclearborder (bw);
    bw1 = imfill( bw1,'holes');
    bw1 = bwareaopen( bw1, 80, 4);
    fluoImMask = bw1;

end

end
%% sub function
function [big_signal] = ImBigSignalGet(signal)
% 荧光图像mask获得用到的子函数
% 用于荧光图像中大光强信号的获得Shuai Yang 2021.10.15

signal = sort(double(signal(:)),'descend' );
pixelSpace = 1:100:numel(signal);%每100个pixel计算平均值
% case for <100 pixel
if numel(pixelSpace) < 2
    big_signal = mean(signal(:));
    return
end

M = zeros(numel(pixelSpace)-1,1);

for i = 1:numel(pixelSpace)-1
    M(i) = mean(signal(pixelSpace(i):pixelSpace(i+1)));
end

if numel(M) < 5
    big_signal = M(1);
    return
end

n = 5;
CV = std(M(1:n))/mean(M(1:n));
while CV > 0.03 % equal CV > 0.03 && n >= 1
    n = n-1;
    CV = std(M(1:n))/mean(M(1:n));
end

if n == 1
    big_signal = mean(M(1:5));
else
    big_signal = mean(M(1:n));
end
end