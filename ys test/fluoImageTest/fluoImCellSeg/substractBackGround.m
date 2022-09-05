function [I,BG] = substractBackGround(I)
%用于减除图像的背景 Shuai Yang 2020.09.09
sz = size(I);
ia = I == 0;% 找出 pixelValue=0 的位置赋值非0的最小值
if sum(ia,'all') ~= 0
    I(ia) = min(I(~ia));
end
I1 = imgaussfilt(I,2);%gaussblur 
I1 = I1(round(sz(1)/2-sz(1)/4):round(sz(1)/2+sz(1)/4), ...
    round(sz(2)/2-sz(2)/4):round(sz(2)/2+sz(2)/4));
I1 = double(I1);
I1 = sort(I1(:));
pixelSpace = 1:100:numel(I1);%每100个pixel计算平均值
M = zeros(numel(pixelSpace)-1,1);
for i = 1:numel(pixelSpace)-1
    M(i) = mean(I1(pixelSpace(i):pixelSpace(i+1)));
end

n = 5;
CV = std(M(1:n))/mean(M(1:n));
while CV > 0.03 % equal CV > 0.03 && n >= 1
    n = n-1;
    CV = std(M(1:n))/mean(M(1:n));
end
BG = mean(M(1:n));
I = I - BG;% 用原始的图像减背景，blur
end