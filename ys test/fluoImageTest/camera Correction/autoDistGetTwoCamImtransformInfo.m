function [transformInfo] = autoDistGetTwoCamImtransformInfo(Im_cam1,Im_cam2,dirFile)
% 通过质心距离来自动获取图像校正点 selectedFixedPoints;selectedMovingPoints
% Im_cam1 相机1拍摄；Im_cam2 相机2拍摄
% Shuai Yang 2022.06.24
% camera 1 拍摄的为fixed channels，camera2 拍摄的为moving channels
% clc
% close all
%y轴镜像对称，2;x轴镜像对称1；
Im_cam2 = flip(Im_cam2,1);%图像镜像对称 Micromanger 2; MATLAB 1

disp('Auto distance get camera transformInfo by two Cam images')
allChannels={'BF1','PVD','Venus','sfGFP','mScarletI','CyOFP','TDsmURFP'};
channels_cam1 = {'BF1','PVD','sfGFP'};% camera 1拍摄的通道 短波长
channels_cam2 = {'mScarletI','CyOFP','Venus','TDsmURFP'};% camera 2拍摄的通道


%初始化
for iChannel = 1:numel(allChannels)
    transformInfo.(allChannels{iChannel}) = NaN;
end

%'PhC','Fluo','Binary'
[fixedImageMask,ImType_cam1] = imageTypeGetAndSeg(Im_cam1); 
disp(['Im_Cam1 is ',ImType_cam1, 'Image'])
[movingImageMask,ImType_cam2] = imageTypeGetAndSeg(Im_cam2);
disp(['Im_Cam2 is ',ImType_cam2, 'Image'])

% 距离阈值可以根据图像偏移大小进行自定义
FpMpdist_thre = 30;% 固定点和移动点之间最大偏移量 unit pixel distance threshold
% neighPtDist_thre = 32;% 确保周围多少范围内没有其它pts unit pixel
neighPtDist_thre = FpMpdist_thre*1.6;

[selectedMovingPoints,selectedFixedPoints] = ...
    autoDistGetCalibrationPoints_ys(fixedImageMask,movingImageMask,FpMpdist_thre,neighPtDist_thre);
tform = getTransformInfoAndCheckResult(selectedMovingPoints,...
    selectedFixedPoints,Im_cam2,Im_cam1);

for iChannel = 1:numel(allChannels)
    if ismember(allChannels{iChannel},channels_cam1)
        transformInfo.(allChannels{iChannel}) = 0;
    end
    if ismember(allChannels{iChannel},channels_cam2)
        transformInfo.(allChannels{iChannel}) = tform;
    end
end
save(strcat(dirFile,'\transformInfo.mat'),'transformInfo');
end

%% 
function [selectedMovingPoints,selectedFixedPoints] = ...
autoDistGetCalibrationPoints_ys(fixedImageMask,movingImageMask,FpMpdist_thre,neighPtDist_thre)
% fixedImage,movingImage;binay images
% Shuai Yang 2022/6/24

% % 距离阈值可以根据图像偏移大小进行自定义
% FpMpdist_thre = 50;% 固定点和移动点之间的质心距离 unit pixel distance threshold
% neighPtDist_thre = 80;% 确保周围多少范围内没有其它pts unit pixel

fp = [];
mp = [];
% movingImageMask = flip(movingImageMask,1);
imsz = size(fixedImageMask);
linkPtsImage = false(imsz);

stats1 = regionprops(fixedImageMask,'Centroid');
stats2 = regionprops(movingImageMask,'Centroid');
centroid_Fixed = cat(1,stats1.Centroid);
centroid_Moving = cat(1,stats2.Centroid);

if numel(stats1) < 12 || numel(stats2) < 12
    disp('No enough points')
    return
end

for iFpts = 1:numel(stats1)
    oneFpt  =  centroid_Fixed(iFpts,:);
    D1 = pdist2(oneFpt,centroid_Fixed);
    D2 = pdist2(oneFpt,centroid_Moving);
    B1 = mink(D1,2);% 找到距离最小的前两位 最小是0 自己
    if B1(2) < neighPtDist_thre % 确保周围距离阈值内没有其他菌
        continue
    end
    [M,I] = min(D2);

    if M > FpMpdist_thre % fix 与 mov pt 距离不能太远
        continue
    end

    if sum(D2 < neighPtDist_thre) > 2 % 距离fp的小于阈值的mp数目只有一个
        continue
    end

    candiMpt = centroid_Moving(I,:);
    D3 = pdist2(candiMpt,centroid_Moving);
    B2 = mink(D3,2);
    if B2(2) < neighPtDist_thre % 确保备选的movPt周围距离阈值内没有其他点
        continue
    end

    fp = cat(1, fp ,oneFpt);
    mp = cat(1, mp ,candiMpt);

    oneFpt = round(oneFpt);
    candiMpt = round(candiMpt);
    xSpace = abs(oneFpt(1)-candiMpt(1))+1;
    ySpace = abs(oneFpt(2)-candiMpt(2))+1;
    npts = xSpace;% number of points
    if xSpace < ySpace
        npts = ySpace;
    end

    x = linspace(oneFpt(1),candiMpt(1),npts);
    y = linspace(oneFpt(2),candiMpt(2),npts);
    x = round(x);y = round(y);
    linkPtsInd = sub2ind(imsz,y,x);
    linkPtsImage(linkPtsInd) = true;
end

% selectedMovingPoints mp;selectedFixedPoints fp；
selectedFixedPoints = fp;
selectedMovingPoints = mp;

tform = fitgeotrans(selectedMovingPoints,selectedFixedPoints,'projective');
Rfixed = imref2d(size(movingImageMask));%fixed image size
cmmIregistered = imwarp(movingImageMask,tform,'OutputView',Rfixed);

C1 = cat(3,uint8(fixedImageMask*255),uint8(movingImageMask*255),uint8(linkPtsImage*255));
figure,imshow(C1);
% C1 = cat(3,fixedImageMask,movingImageMask,linkPtsImage);
% C2 = cat(3,fixedImageMask,cmmIregistered,fixedImageMask*0);
% figure,
% imshowpair(C1,C2,'montage')
% title('befor correction (left)- afer correction (right)')

end

