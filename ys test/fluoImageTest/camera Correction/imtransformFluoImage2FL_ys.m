function transformInfo = imtransformFluoImage2FL_ys(dirFile,fixedChannel)
%�˺���������channel����У����
%��ĿǰIP33������������ͬһ��������ͼ���غϺܺã����ֻ��Ҫ���У����������Ҫ����ͨ��У��
%���imtransformFluoImage2FL_basedOnCam 2020.05.13 ys
%transformInfoУ����Ϣ���
% ���ڲ�ͬͨ����ӫ��ͼ����� ��Ҫѡ��һ��fixchannel
%������Ϊmovingchannel 
%����channel��CyOFP��ͼ��Ϊ׼���ж��� 2020.04.01 ys
disp('Align fluorescent imgaes that are from different channel')
allChannels={'BF1','PVD','Venus','sfGFP','mScarletI','CyOFP','TDsmURFP'};
channels_cam1 = {'BF1','PVD','sfGFP'};% camera 1�����ͨ�� �̲���
channels_cam2 = {'mScarletI','CyOFP','Venus','TDsmURFP'};% camera 2�����ͨ��
%������PhC���ͼ���������mask��ͨ����ʱΪ����ΪBF1��camera 1 ����2020.05.13 ys

if ismember(fixedChannel,allChannels)
    disp([fixedChannel, 32, 'used as the fixed channel'])
else
    disp ( 'input error of fixedChannel')
    return
end

if ismember(fixedChannel,channels_cam1)
    fixedChanLib = channels_cam1;%�̶���channels library
else
    fixedChanLib = channels_cam2;
end
for iChannel = 1:numel(allChannels)
    transformInfo.(allChannels{iChannel}) = NaN;
end
transformInfo.(fixedChannel) = 0;
% eval(['transformInfo.',fixedChannel,'= 0;']);%��fixedChannel��ֵ����Ϊ0��û�������ͨ��ΪNaN

myfield = 1;
myframe = 1;
dirField = strcat(dirFile,'\field',num2str(myfield,'%.4d'));
channelList = dir(dirField);

fixedImage=import_tiff_stack(strcat(dirField,'\',fixedChannel,'\image',fixedChannel,num2str(myframe,'%.5d'),'.tif'));
if strcmp(fixedChannel,'BF1')
    [fixedImageMask] = phaseContrastImageProcessing_ys(fixedImage);
else
    [fixedImageMask] = fluoImageProcessing_ystest(fixedImage); %��fixedImageΪ��׼ ����������Ұ�����ƫ����
    % 10 Ϊintensity threshold
end


for iChannel = 1 : length(channelList)
    if(isequal(channelList(iChannel).name,'.')||... % ȥ��ϵͳ�Դ����������ļ���
            isequal(channelList(iChannel).name,'..')||...
            ~channelList(iChannel).isdir)                % ȥ�������в����ļ��е�
        continue;
    end
    
    if ismember(channelList(iChannel).name,allChannels ) && ~strcmp(channelList(iChannel).name,fixedChannel)
        imageChannel=channelList(iChannel).name;        
        switch imageChannel
            case 'BF1'
                imageBF1=import_tiff_stack(strcat(dirField,'\BF1\imageBF1',num2str(myframe,'%.5d'),'.tif'));
                if ~ismember(imageChannel,fixedChanLib) 
                     imageBF1 = flip(imageBF1,2);%�����������ͼ��y�᾵��Գ�
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
                    imagemScarletI = flip(imagemScarletI,2);%�����������ͼ��y�᾵��Գ�
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
                    imagesfGFP = flip(imagesfGFP,2);%�����������ͼ��y�᾵��Գ�
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
                    imageCyOFP = flip(imageCyOFP,2);%�����������ͼ��y�᾵��Գ�
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
                    imageVenus = flip(imageVenus,2);%�����������ͼ��y�᾵��Գ�
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
                    imagePVD = flip(imagePVD,2);%�����������ͼ��y�᾵��Գ�
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
                    imageTDsmURFP = flip(imageTDsmURFP,2);%�����������ͼ��y�᾵��Գ�
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
% selectedMovingPoints mp;selectedFixedPoints fp��
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
% ���ӫ��ͼ���ʱ�����е�xy����
% �������Ϊ��ֵͼ�񣬺��Ľ��������Ѿ��Ż���
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
    image2New=imdilate(image2,se);  % ����imdilateʵ��ƽ��
    sumImage=image1 & image2New;    % �����߼�����ĳ˷��൱��&
    correlationMatrix(i)=sum(sum(sumImage)); % ֱ�Ӷ��߼���������ٶȱȽϿ�
end
bestPosition=[x(correlationMatrix==max(correlationMatrix)),y(correlationMatrix==max(correlationMatrix))];
bestPosition=bestPosition(1,:);
end
function image=imageCorrectionWithBestPosition(image,bestPosition)
% ��֪Ư��������еĽ���
bestPositionAccumulation=zeros(size(bestPosition));
for i=2:size(bestPosition,1)
    bestPositionAccumulation(i,:)=bestPosition(i,:)+bestPositionAccumulation(i-1,:);
end
for i=2:size(bestPosition,1)
    se=translate(strel(1),bestPositionAccumulation(i,:));
    image(:,:,i)=imdilate(image(:,:,i),se);
end
end