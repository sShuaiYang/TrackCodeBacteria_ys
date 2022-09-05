function [fluoImMask] = fluoImCellMask_edge_w60x(fluoImage)
%% 利用edge图像进行边缘检测
% 用edge获得mask
% 原先code处理60x water objective 数据有问题 改进
% water60xFluoImMask_ys
% Shuai Yang 2021/12/28

fluoImMask = false(size(fluoImage));
I0 = fluoImage;
[~,BG] = substractBackGround(I0);%获得背景数值 backGround signal
small = BG;
big = ImBigSignalGet(I0);

sigma =  1;% gaussfilt sigma; could change

% edge Filter 看情况使用 细菌内部不均匀100x图像可能会出现多余的edge线条
% edgeFilter = (ones(5,5)).*-1;edgeFilter(3,3) = 24;%created edgefilter
I1 = imgaussfilt(I0-BG, sigma);
% I1 = imfilter(I1,edgeFilter);
% I1 = imgaussfilt(I1, sigma);


if big > small %weak fluo cause err

    thre = prctile((small:big),10); % 百分位数
    bw = edge(I1-thre,'log',0);

    bw1 = imclearborder (bw);
    bw1 = imfill( bw1,'holes');
    bw1 = bwareaopen( bw1, 40, 4);
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