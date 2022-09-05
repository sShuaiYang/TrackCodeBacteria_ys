function [transformInfo] = getTwoCamImtransformInfo(Im_cam1,Im_cam2,dirFile)
%根据两个相机拍摄的明场(一般为相差PhC)获得相机的transform 信息
% Im_cam1 相机1拍摄；Im_cam2 相机2拍摄
% Shuai Yang 2020.05.27
% camera 1 拍摄的为fixed channel，camera2 拍摄的为moving channel

warning on
isDualCam = true;
% true or false 判断是否双相机拍摄 双相机图像需要旋转 单相机不需要
% 双相机肯定需要校正 单相机不一定
if isDualCam
    %y轴镜像对称，2;x轴镜像对称1；
    Im_cam2 = flip(Im_cam2,1);%图像镜像对称%matlab 垂直
    % Im_cam2 = flip(Im_cam2,2); %用micromanger拍摄时水平镜像对称
end

disp('Get camera transformInfo by two Cam images')
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


[selectedMovingPoints,selectedFixedPoints]=...
    getCalibrationPoints_ys(fixedImageMask,movingImageMask);
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
function [selectedMovingPoints,selectedFixedPoints] = getCalibrationPoints_ys(fixedImage,movingImage)
% fixedImage movingImage 为binary图像 二值图像
% Shuai Yang 2022/6/15
imsz = size(fixedImage);
stats1 = regionprops(fixedImage,'Centroid');
stats2 = regionprops(movingImage,'Centroid');
centroid_Fixed = cat(1,stats1.Centroid);
centroid_Moving = cat(1,stats2.Centroid);

centerMarkFixedImage = false(imsz);
centerMarkMovingImage = false(imsz);

ind_centerFixed = sub2ind(imsz,round(centroid_Fixed(:,2)),round(centroid_Fixed(:,1)));
ind_centerMoving = sub2ind(imsz,round(centroid_Moving(:,2)),round(centroid_Moving(:,1)));

centerMarkFixedImage(ind_centerFixed) = true;  
centerMarkMovingImage(ind_centerMoving)= true;  

markedFixedImage = imoverlay(uint8(fixedImage*255),centerMarkFixedImage,[0,1,0]);
markedMovingImage = imoverlay(uint8(movingImage*255),centerMarkMovingImage,[0,1,0]);

[mp,fp] = cpselect(markedMovingImage,markedFixedImage,'Wait',true);

% selectedMovingPoints mp;selectedFixedPoints fp；
selectedMovingPoints = mp;
selectedFixedPoints = fp;

[fp_c,mp_c] = cptsCorrectionBycellMask(fp,mp,fixedImage,movingImage, ...
    stats1,stats2);

selectedMovingPoints = mp_c;
selectedFixedPoints = fp_c;

end
%% calibration points correction
function [fp_c,mp_c] = cptsCorrectionBycellMask(fp,mp,fixedImage,movingImage, ...
    stats1,stats2)
% 不再要求精确选择细菌的质心了 只需要保证选在细菌内部就可以
% Shuai Yang 2022/6/15

L_Fixed = bwlabel(fixedImage);
L_Moving = bwlabel(movingImage);

centroid_Fixed = cat(1,stats1.Centroid);
centroid_Moving = cat(1,stats2.Centroid);

fp_c = zeros(size(fp)); % correct
mp_c = zeros(size(mp));
for iSpts = 1:size(fp,1) % selected points

    % 根据最近的质心位置 对选取的点进行校正
    D1 = pdist2(fp(iSpts,:),centroid_Fixed);
    k1 = D1 == min(D1);
    fp_c(iSpts,:) = stats1(k1).Centroid;

    D2 = pdist2(mp(iSpts,:),centroid_Moving);
    k2 = D2 == min(D2);
    mp_c(iSpts,:) = stats2(k2).Centroid;

    % 判断校正的点和最初选取的点是不是一个细菌 根据bwlabel的数值
    bw1 = false(size(fixedImage));
    bw2 = false(size(movingImage));
    fp_check = cat(1,fp(iSpts,:),fp_c(iSpts,:));
    mp_check = cat(1,mp(iSpts,:),mp_c(iSpts,:));

    ind1 = sub2ind(size(bw1),round(fp_check(:,2)),round(fp_check(:,1)));
    bw1(ind1) = true;
    ind2 = sub2ind(size(bw2),round(mp_check(:,2)),round(mp_check(:,1)));
    bw2(ind2) = true;
    n1 = L_Fixed(bw1);
    n2 = L_Moving(bw2);

    % unique(n1) 唯一且等于 find(k1);unique(n2) 唯一且等于 find(k2);
    % 保证校正和最初选择的点是唯一一个菌 
    if sum(find(k1)-n1) ~= 0|| sum(find(k2)-n2) ~= 0
        msg = 'Some of selectedpoints are not located inside the cellmask';
        warning(msg)
    end
end
end


