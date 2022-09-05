function timeLapsexyShiftCorrectionForMultiFluoChannels(dirFile,fixedChannel)
% time lapse 拍摄 荧光图像会出现xy的漂移
%  通过计算fixedChannel的偏移量 进行校正
%认为其他通道的xy偏移量与fixedChannel一致
% Shuai Yang 2020.04.03
disp('xy shift correction for time lapse images of different fluo channels')
allChannels = {'BF1','sfGFP','mScarletI','Venus','PVD','CyOFP','TDsmURFP'};
% fluoChannels = {'sfGFP','mScarletI','Venus','PVD','CyOFP','TDsmURFP'};

if ismember(fixedChannel,allChannels)
    disp([fixedChannel, 32, 'used as the fixed channel'])
else
    disp ( 'input error of fixedChannel')
    return
end

fieldList = dir(dirFile);
for iField = 1:length(fieldList)-2
    if ~strcmp(fieldList(iField+2).name(1:5),'field')
        continue
    end
    disp (fieldList(iField+2).name);
    dirField = strcat(dirFile,'\',fieldList(iField+2).name);
    dirFixedChannel = strcat(dirField,'\',fixedChannel);
    fixedImList = dir(dirFixedChannel);
    trackingNum = length(fixedImList)-3;
    bestPositionSeries = timeLapseBestPositionGet(fixedChannel,dirFixedChannel);
    load([dirFixedChannel,'\frameInfo.mat'])
    frameInfo_tracking = frameInfo;

    for iChannel = 1:numel(allChannels)
%         if strcmp(fluoChannels{iChannel},fixedChannel)
%             continue
%         end % fixed channels 也放到这一步一起校正
        dirImage=[dirField,'\',allChannels{iChannel}];
        try
            load([dirImage,'\frameInfo.mat'])
        catch ME
            continue
        end
        
        frameInfo_fluo = frameInfo;
        disp (['Start',32,allChannels{iChannel},32,'correciton']);
        imageList = dir(dirImage);
        imageNum = length(imageList)-3;
        fluoImStack = zeros(2048,2048,imageNum,'uint16');
        newPositionSeries = newBestPositionGetForUnsynFluo(trackingNum,bestPositionSeries,frameInfo_tracking,frameInfo_fluo);
        newPositionSeries = newPositionSeries(1:imageNum,:);
        parfor iImage= 1:imageNum
            fluoImStack (:,:,iImage) = import_tiff_stack(strcat(dirImage,'\',imageList(iImage+3).name));
        end
        fluoImStack = imageCorrectionWithBestPosition(fluoImStack,newPositionSeries);
        
        dirNewImage = strcat(dirField,'\',allChannels{iChannel},'_xycorrected');
        mkdir(dirNewImage)
        
        for  iImage = 1:imageNum
            imageFluo = fluoImStack(:,:,iImage);
            %                         imwrite(imageFluo,strcat(dirImage,'\',imageList(iImage+3).name));
            imwrite(imageFluo,strcat(dirNewImage,'\',imageList(iImage+3).name));
        end
        %                 load([dirImage,'\frameInfo.mat'])
        %                 save([dirNewImage,'\','frameInfo.mat'],'frameInfo');
        
        
    end
end

end

%%
function bestPositionSeries = timeLapseBestPositionGet(fixedChannel,dirFixedChannel)
% 根据fixedchannel来计算bestPositionSeries
%然后先对fixedchannel的图像进行xy shift校正
% Shuai Yang
% disp(['BestPositionSeries get and',32,fixedChannel,32, 'correction']);
disp(['BestPositionSeries get of ',32,fixedChannel,32, 'Channel'])
fixedImList = dir(dirFixedChannel);
imageNum = length(fixedImList)-3;

imageType='uint16';
imStack = zeros(2048,2048,imageNum,imageType);

parfor iImage= 1 : imageNum
    imStack (:,:,iImage)=import_tiff_stack(strcat(dirFixedChannel,'\',fixedImList(iImage+3).name));
end

if strcmp(fixedChannel,'BF1')
    [imStackMask] = phaseContrastImageProcessing_ys(imStack);
else
    [imStackMask] = fluoImageProcessing_ystest(imStack); %以fixedImage为基准 计算其他视野的相机偏移量
end

[~,bestPositionSeries] = fluoImageCorrection(imStackMask);
% imStack = imageCorrectionWithBestPosition(imStack,bestPositionSeries);
% for  iImage = 1:imageNum
%     fixedImage = imStack(:,:,iImage);
%     imwrite(fixedImage,strcat(dirFixedChannel,'\',fixedImList(iImage+3).name));
% end

clear imageStack imageStackMask fixedImage

end
function newPositionSeries = newBestPositionGetForUnsynFluo(trackingNum,bestPositionSeries,frameInfo_tracking,frameInfo_fluo)
if size(frameInfo_tracking,1) == size(frameInfo_fluo,1)
    newPositionSeries = bestPositionSeries;
    return;
end
% 寻找tracking的时间窗口内，所拍摄的荧光图像的数目
fluo2TrackingIdx = zeros (1, size(frameInfo_fluo,1));
for iFluo = 1 :size (frameInfo_fluo,1)
    time_lag = abs(etime (frameInfo_fluo(iFluo,1:6),frameInfo_tracking(:,1:6))/60);
    %min 与frameInfo_tracking拍摄的时间间隔绝对值
    trckingIdx = find(time_lag==min(time_lag));
    fluo2TrackingIdx(iFluo) = trckingIdx(end);
    %use trckingIdx(end),for the case of trckingIdx 不唯一 ;case:1,1,2,2...
end
fluo2TrackingIdx = fluo2TrackingIdx(fluo2TrackingIdx <= trackingNum);
[~,checkFrame] = find(diff(fluo2TrackingIdx)==0);
fluo2TrackingIdx(checkFrame+1) = 0;% case for:两张荧光同时对应一张tracking
%赋值为0，数组索引就会出错，便于debug
newPositionSeries = bestPositionSeries(fluo2TrackingIdx,:);
end
%% xyshift correction
function [imageStack,bestPosition] = fluoImageCorrection(imageStack)
% 针对荧光图像的时间序列的xy矫正
% 输入变量为二值图像，核心矫正函数已经优化过
% imageSize=size(imageStack);
% if min(imageSize(1:2))>=1500
%     imageStackNew=imageStack(1024-500:1024+500,1024-500:1024+500,:);
% %     imageStackNew=imageStack(600:1100,600:1100,:);
% else
%     imageStackNew=imageStack;
% end
[imageStackNew] = cropImageForxyShiftCorrection(imageStack,400);% 2020,01.03ys
imageStackNew1(:,:,2:size(imageStackNew,3)+1)=imageStackNew;
parfor i=2:size(imageStack,3)
    bestPosition(i,:)=caculateCrossCorrelationForImage(imageStackNew1(:,:,i),imageStackNew(:,:,i),15);%偏移最大[-15,15];ys2020.01.03
end
% imageStack=imageCorrectionWithBestPosition(imageStack,bestPosition);
% bestPositionAccumulation=zeros(size(bestPosition));
% for i=2:size(bestPosition,1)
%     bestPositionAccumulation(i,:)=bestPosition(i,:)+bestPositionAccumulation(i-1,:);
% end
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
function [imageStackNew]=cropImageForxyShiftCorrection(imageStack,pixelRange)
%imageStack为二值图像 2020.01.03ys 寻找图像切割的范围
%   pixelRange=250;
imageSize = size(imageStack);
if min(imageSize(1:2))>=1000
    xCrop = [fix(imageSize(1)/2)-pixelRange,fix(imageSize(1)/2)+pixelRange];
    yCrop = [fix(imageSize(2)/2)-pixelRange,fix(imageSize(2)/2)+pixelRange];
    tempStack = imageStack(xCrop(1):xCrop(2),yCrop(1):yCrop(2),:);
    stats = regionprops(tempStack(:,:,1),'PixelList');
    if numel(stats)<7
        pixelRange = pixelRange+200;
        if pixelRange>=min(imageSize(1:2))/2
            imageStackNew = imageStack;
        else
            [imageStackNew] = cropImageForxyShiftCorrection(imageStack,pixelRange);
        end
    else
        imageStackNew = tempStack;
    end
    
else
    imageStackNew = imageStack;
end

end
