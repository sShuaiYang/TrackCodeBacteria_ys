function [fluoImMask] = fluoImCellMask_Otsu_w60x(fluoImage)
%% 利用Otsu图像进行mask获得
% 用edge获得mask
% 原先code处理60x water objective 数据有问题 改进
% water60xFluoImMask_ys
% Shuai Yang 2021/12/28

I0 = fluoImage;
[~,BG] = substractBackGround(I0);%获得背景数值 backGround signal

sigma =  1;% gaussfilt sigma; could change
edgeFilter = (ones(5,5)).*-1;edgeFilter(3,3) = 24;%created edgefilter

I1 = imgaussfilt(I0-BG, sigma);
I1 = imfilter(I1,edgeFilter);
I1 = imgaussfilt(I1, sigma);

level = graythresh(I1); % OtSu threshold get
bw = imbinarize(I1,level/2.6);% 2.6 could change
bw1 = imclearborder (bw);
bw1 = imfill( bw1,'holes');
bw1 = bwareaopen( bw1, 40, 4);
fluoImMask = bw1;

end