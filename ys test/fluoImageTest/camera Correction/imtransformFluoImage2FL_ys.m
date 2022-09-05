function transformInfo = imtransformFluoImage2FL_ys(dirFile,fixedChannel)
%此函数对所有channel进行校正，
%从目前IP33的数据来看，同一相机拍摄的图像重合很好，因此只需要相机校正，并不需要所有通道校正
%详见imtransformFluoImage2FL_basedOnCam 2020.05.13 ys
%transformInfo校正信息获得
% 用于不同通道的荧光图像对齐 需要选择一个fixchannel
%其他则为movingchannel 
%其他channel以CyOFP的图像为准进行对齐 2020.04.01 ys
disp('Align fluorescent imgaes that are from different channel')
allChannels={'BF1','PVD','Venus','sfGFP','mScarletI','CyOFP','TDsmURFP'};
channels_cam1 = {'BF1','PVD','sfGFP'};% camera 1拍摄的通道 短波长
channels_cam2 = {'mScarletI','CyOFP','Venus','TDsmURFP'};% camera 2拍摄的通道
%加入了PhC相差图像可以来做mask，通道暂时为命名为BF1，camera 1 拍摄2020.05.13 ys

if ismember(fixedChannel,allChannels)
    disp([fixedChannel, 32, 'used as the fixed channel'])
else
    disp ( 'input error of fixedChannel')
    return
end

if ismember(fixedChannel,channels_cam1)
    fixedChanLib = channels_cam1;%固定的channels library
else
    fixedChanLib = channels_cam2;
end
for iChannel = 1:numel(allChannels)
    transformInfo.(allChannels{iChannel}) = NaN;
end
transformInfo.(fixedChannel) = 0;
% eval(['transformInfo.',fixedChannel,'= 0;']);%将fixedChannel的值设置为0；没有拍摄的通道为NaN

myfield = 1;
myframe = 1;
dirField = strcat(dirFile,'\field',num2str(myfield,'%.4d'));
channelList = dir(dirField);

fixedImage=import_tiff_stack(strcat(dirField,'\',fixedChannel,'\image',fixedChannel,num2str(myframe,'%.5d'),'.tif'));
if strcmp(fixedChannel,'BF1')
    [fixedImageMask] = phaseContrastImageProcessing_ys(fixedImage);
else
    [fixedImageMask] = fluoImageProcessing_ystest(fixedImage); %以fixedImage为基准 计算其他视野的相机偏移量
    % 10 为intensity threshold
end


for iChannel = 1 : length(channelList)
    if(isequal(channelList(iChannel).name,'.')||... % 去除系统自带的两个隐文件夹
            isequal(channelList(iChannel).name,'..')||...
            ~channelList(iChannel).isdir)                % 去除遍历中不是文件夹的
        continue;
    end
    
    if ismember(channelList(iChannel).name,allChannels ) && ~strcmp(channelList(iChannel).name,fixedChannel)
        imageChannel=channelList(iChannel).name;        
        switch imageChannel
            case 'BF1'
                imageBF1=import_tiff_stack(strcat(dirField,'\BF1\imageBF1',num2str(myframe,'%.5d'),'.tif'));
                if ~ismember(imageChannel,fixedChanLib) 
                     imageBF1 = flip(imageBF1,2);%两个相机拍摄图像y轴镜像对称
                end
                maskimageBF1 = phaseContrastImageProcessing_ys(imageBF1);
                [selectedMovingPoints, selectedFixedPoints] = ...
                    getCalibrationPoints_ys(fixedImageMask,maskimageBF1);
                tform=getTransformInfoandCheckResult(selectedMovingPoints,...
                    selectedFixedPoints,maskimageBF1,fixedImageMask);
                transformInfo.BF1 = tform;
                
            case 'mScarletI'
                imagemScarletI = import_tiff_stack(strcat(dirField,'\mScarletI\imagemScarlet',num2str(myframe,'%.5d'),'.tif'));
                if ~ismember(imageChannel,fixedChanLib)
                    imagemScarletI = flip(imagemScarletI,2);%两个相机拍摄图像y轴镜像对称
                end
                maskImagemScarletI = fluoImageProcessing_ystest(imagemScarletI);
                [selectedMovingPoints, selectedFixedPoints] = ...
                    getCalibrationPoints_ys(fixedImageMask,maskImagemScarletI);
                tform=getTransformInfoandCheckResult(selectedMovingPoints,...
                    selectedFixedPoints,maskImagemScarletI,fixedImageMask);   
                transformInfo.mScarletI = tform;
                
            case 'sfGFP'
                imagesfGFP = import_tiff_stack(strcat(dirField,'\sfGFP\imagesfGFP',num2str(myframe,'%.5d'),'.tif'));
                if ~ismember(imageChannel,fixedChanLib)
                    imagesfGFP = flip(imagesfGFP,2);%两个相机拍摄图像y轴镜像对称
                end
                maskImagesfGFP = fluoImageProcessing_ystest(imagesfGFP);
                [selectedMovingPoints, selectedFixedPoints]=...
                    getCalibrationPoints_ys(fixedImageMask,maskImagesfGFP);
                tform = getTransformInfoandCheckResult(selectedMovingPoints,...
                    selectedFixedPoints,maskImagesfGFP,fixedImageMask);
                transformInfo.sfGFP = tform;
                
            case 'CyOFP'
                imageCyOFP = import_tiff_stack(strcat(dirField,'\CyOFP\imageCyOFP',num2str(myframe,'%.5d'),'.tif'));
                if ~ismember(imageChannel,fixedChanLib)
                    imageCyOFP = flip(imageCyOFP,2);%两个相机拍摄图像y轴镜像对称
                end
                maskImageCyOFP = fluoImageProcessing_ystest(imageCyOFP);
                [selectedMovingPoints, selectedFixedPoints]=...
                    getCalibrationPoints_ys(fixedImageMask,maskImageCyOFP);
                tform = getTransformInfoandCheckResult(selectedMovingPoints,...
                    selectedFixedPoints,maskImageCyOFP,fixedImageMask);
                transformInfo.CyOFP = tform;
                
            case 'Venus'
                imageVenus = import_tiff_stack(strcat(dirField,'\Venus\imageVenus',num2str(myframe,'%.5d'),'.tif'));
                if ~ismember(imageChannel,fixedChanLib)
                    imageVenus = flip(imageVenus,2);%两个相机拍摄图像y轴镜像对称
                end
                maskImageVenus = fluoImageProcessing_ystest(imageVenus); 
                [selectedMovingPoints,selectedFixedPoints]=...
                    getCalibrationPoints_ys(fixedImageMask,maskImageVenus);
                tform = getTransformInfoandCheckResult(selectedMovingPoints,...
                    selectedFixedPoints,maskImageVenus,fixedImageMask);
                transformInfo.Venus = tform;
                
            case 'PVD'
                imagePVD = import_tiff_stack(strcat(dirField,'\Venus\imageVenus',num2str(myframe,'%.5d'),'.tif'));
                if ~ismember(imageChannel,fixedChanLib)
                    imagePVD = flip(imagePVD,2);%两个相机拍摄图像y轴镜像对称
                end
                maskImagePVD = fluoImageProcessing_ystest(imagePVD);
                [selectedMovingPoints,selectedFixedPoints]=...
                    getCalibrationPoints_ys(fixedImageMask,maskImagePVD);
                tform = getTransformInfoandCheckResult(selectedMovingPoints,...
                    selectedFixedPoints,maskImagePVD,fixedImageMask);
                transformInfo.PVD = tform;
                
            case'TDsmURFP'
                imageTDsmURFP = import_tiff_stack(strcat(dirField,'\TDsmURFP\imageTDsmURFP',num2str(myframe,'%.5d'),'.tif'));
                if ~ismember(imageChannel,fixedChanLib)
                    imageTDsmURFP = flip(imageTDsmURFP,2);%两个相机拍摄图像y轴镜像对称
                end
                maskImageTDsmURFP = fluoImageProcessing_ystest(imageTDsmURFP);
                [selectedMovingPoints,selectedFixedPoints]=...
                    getCalibrationPoints_ys(fixedImageMask,maskImageTDsmURFP);
                tform = getTransformInfoandCheckResult(selectedMovingPoints,...
                    selectedFixedPoints,maskImageTDsmURFP,fixedImageMask);
                transformInfo.TDsmURFP = tform;
        end
    end
end
save(strcat(dirFile,'\transformInfo.mat'),'transformInfo');
end

%%
function transformInfo=getTransformInfoandCheckResult(selectedMovingPoints,selectedFixedPoints,movingImageMask,fixedImageMask)

tform = fitgeotrans(selectedMovingPoints,selectedFixedPoints,'projective');
Rfixed = imref2d(size(movingImageMask));%fixed image size
mIregistered = imwarp(movingImageMask,tform,'OutputView',Rfixed);%moving image registered
registeredCheckImage=segrationImage(fixedImageMask,mIregistered);
figure,imshow(registeredCheckImage);
transformInfo= tform;

end
%%
function [selectedMovingPoints,selectedFixedPoints]=getCalibrationPoints_ys(fixedImage,movingImage)

stats1 = regionprops(fixedImage,'Centroid');
stats2 = regionprops(movingImage,'Centroid');

centerMarkFixedImage = false(size(fixedImage));
centerMarkMovingImage = false(size(movingImage));

for i= 1:numel(stats1)
    centerMarkFixedImage(round(stats1(i).Centroid(2)),round(stats1(i).Centroid(1))) = true;    
end
for i= 1:numel(stats2)
    centerMarkMovingImage(round(stats2(i).Centroid(2)),round(stats2(i).Centroid(1))) = true;    
end

markedFixedImage = segrationImage(fixedImage,centerMarkFixedImage);
markedMovingImage = segrationImage(movingImage,centerMarkMovingImage);

[mp,fp] = cpselect(markedMovingImage,markedFixedImage,'Wait',true);
% selectedMovingPoints mp;selectedFixedPoints fp；
selectedMovingPoints = mp;
selectedFixedPoints = fp;

end


%%
function  afterProcessingImages=segrationImage(beforeProcessingImages,maskImages) % this function can overlay your mask and orignal images 
maskColor=[0,1,0];
imageType='uint8'; %here you can change your image type
afterProcessingImages=zeros(size(beforeProcessingImages,1),size(beforeProcessingImages,2),3,size(beforeProcessingImages,3),imageType);
parfor iframe=1:size(beforeProcessingImages,3)
%     afterProcessingImages(:,:,iframe)=immultiply((~bwperim(maskImages(:,:,iframe),4)),beforeProcessingImages(:,:,iframe));
%      afterProcessingImages(:,:,iframe)=immultiply(~maskImages(:,:,iframe),beforeProcessingImages(:,:,iframe));
afterProcessingImages(:,:,:,iframe)=imoverlay(beforeProcessingImages(:,:,iframe),maskImages(:,:,iframe),maskColor);
end
end
%%
%xy shift correction
function [imageStack,bestPosition]= fluoImageCorrection(imageStack)
% 针对荧光图像的时间序列的xy矫正
% 输入变量为二值图像，核心矫正函数已经优化过
imageSize=size(imageStack);
if min(imageSize(1:2))>=1500
    imageStackNew=imageStack(200:1500,200:1500,:);
else
    imageStackNew=imageStack;
end
imageStackNew1(:,:,2:size(imageStackNew,3)+1)=imageStackNew;
parfor i=2:size(imageStack,3)
    bestPosition(i,:)=caculateCrossCorrelationForImage(imageStackNew1(:,:,i),imageStackNew(:,:,i),15);
end
imageStack=imageCorrectionWithBestPosition(imageStack,bestPosition);
% gfpImage=imageCorrectionWithBestPosition(gfpImage,bestPosition);
% rfpImage=imageCorrectionWithBestPosition(rfpImage,bestPosition);
end
function bestPosition=caculateCrossCorrelationForImage(image1,image2,step)
% for calculater the cross correlation of two image
% step is the searching range
% backGround should be calculater or set by user
[x,y]=meshgrid((-step:step)',(-step:step)');
x=x(:);
y=y(:);
correlationMatrix=zeros(size(x,1),1);
parfor i=1:numel(x)
    se=translate(strel(1),[x(i),y(i)]);
    image2New=imdilate(image2,se);  % 巧用imdilate实现平移
    sumImage=image1 & image2New;    % 利用逻辑矩阵的乘法相当于&
    correlationMatrix(i)=sum(sum(sumImage)); % 直接对逻辑矩阵求和速度比较快
end
bestPosition=[x(correlationMatrix==max(correlationMatrix)),y(correlationMatrix==max(correlationMatrix))];
bestPosition=bestPosition(1,:);
end
function image=imageCorrectionWithBestPosition(image,bestPosition)
% 已知漂移量后进行的较正
bestPositionAccumulation=zeros(size(bestPosition));
for i=2:size(bestPosition,1)
    bestPositionAccumulation(i,:)=bestPosition(i,:)+bestPositionAccumulation(i-1,:);
end
for i=2:size(bestPosition,1)
    se=translate(strel(1),bestPositionAccumulation(i,:));
    image(:,:,i)=imdilate(image(:,:,i),se);
end
end