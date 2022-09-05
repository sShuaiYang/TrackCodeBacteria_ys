function [transformInfo] = imtransformInfoGetByTwoOriImages(Im_cam1,Im_cam2,dirFile)
% 根据获得原始图像转化为8bit图像，直接展示进行选点校正 不在标记质心
% Shuai Yang 2021.10.10
%根据两个相机拍摄的图像获得相机的transform 信息
% BF_cam1 相机1拍摄；BF_cam2 相机2拍摄
% camera 1 拍摄的为fixed channel，camera2 拍摄的为moving channel
% Im_cam1 = import_tiff_stack('');
% Im_cam2 = import_tiff_stack('');
%y轴镜像对称，2;x轴镜像对称1;
Im_cam2 = flip(Im_cam2,1);%图像镜像对称%用micromanger拍摄时水平镜像对称%matlab 垂直
disp('Get camera transformInfo by two Original images')
allChannels = {'BF1','PVD','Venus','sfGFP','mScarletI','CyOFP','TDsmURFP'};
channels_cam1 = {'BF1','PVD','sfGFP'};% camera 1拍摄的通道 短波长
channels_cam2 = {'mScarletI','CyOFP','Venus','TDsmURFP'};% camera 2拍摄的通道
%初始化
for iChannel = 1:numel(allChannels)
    transformInfo.(allChannels{iChannel}) = NaN;
end

fixedImage = Im_cam1;
movingImage = Im_cam2;
% 转化为8bit 图像
Im_cam1 = rescale(double(Im_cam1))*255;
Im_cam1 = uint8(Im_cam1);
Im_cam2 = rescale(double(Im_cam2))*255;
Im_cam2 = uint8(Im_cam2);

[mp,fp] = cpselect(Im_cam2,Im_cam1,'Wait',true);
% selectedMovingPoints mp;selectedFixedPoints fp；
selectedMovingPoints = mp;
selectedFixedPoints = fp;

tform = getTransformInfoAndCheckResult(selectedMovingPoints,...
    selectedFixedPoints,movingImage,fixedImage);

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
function transformInfo = getTransformInfoAndCheckResult(selectedMovingPoints,selectedFixedPoints,movingImage,fixedImage)

tform = fitgeotrans(selectedMovingPoints,selectedFixedPoints,'projective');
Rfixed = imref2d(size(movingImage));%fixed image size
mIregistered = imwarp(movingImage,tform,'OutputView',Rfixed);%moving image registered

transformInfo = tform;

I0 = fixedImage - 100;
I0 = uint8(double(I0-min(I0(:)))/double(max(I0(:))-min(I0(:)))*255);
fixedImage = I0;

I0 = mIregistered - 100;
I0 = uint8(double(I0-min(I0(:)))/double(max(I0(:))-min(I0(:)))*255);
mIregistered = I0;
% C = imfuse(fixedImage,mIregistered,'falsecolor','Scaling','joint','ColorChannels',[2 1 0]);
C = cat(3,fixedImage,mIregistered,fixedImage*0);
figure, imshow(C,[])
end
