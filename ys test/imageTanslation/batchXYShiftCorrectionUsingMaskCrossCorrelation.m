function [XYtransArray] = batchXYShiftCorrectionUsingMaskCrossCorrelation(dirFile,fixedChannel)
% time lapse 拍摄 荧光图像会出现xy的漂移
%  通过计算fixedChannel的偏移量 进行校正
%认为其他通道的xy偏移量与fixedChannel一致
%对图像生成mask 利用mask的图像相关计算偏移量
% Shuai Yang 2021.07.01

disp('XY shift correction using Mask Cross Correlation')

allChannels = {'BF1','sfGFP','mScarletI','Venus','PVD','CyOFP','TDsmURFP'};
% fluoChannels = {'sfGFP','mScarletI','Venus','PVD','CyOFP','TDsmURFP'};

if ismember(fixedChannel,allChannels)
    disp([fixedChannel, 32, 'used as the fixed channel'])
else
    disp ( 'input error of fixedChannel')
    return
end

fieldList = dir([dirFile,filesep,'field*']);
for iField = 1:length(fieldList)
    disp (fieldList(iField).name);
    
    dirField = strcat(dirFile,'\',fieldList(iField).name);    
    
    dirFixedChannel = strcat(dirField,'\',fixedChannel); 
    load([dirFixedChannel,'\frameInfo.mat'])
    frameInfo_fixed = frameInfo;
    
    imList = dir([dirFixedChannel,filesep,'*.tif']);    
    priorImList = imList(1:end-1);
    afterImList = imList(2:end);
    fixedImNum = length(imList);
    
    xyShift_Array = zeros(numel(imList)-1,2);
    for iImage = 1:numel(priorImList)
        priorIm = imread([dirField,filesep,fixedChannel,filesep,priorImList(iImage).name]);
        afterIm = imread([dirField,filesep,fixedChannel,filesep,afterImList(iImage).name]);
        if strcmp(fixedChannel,'BF1')
            priorImMask = phCImProcessing_supperSegger(priorIm);
            afterImMask = phCImProcessing_supperSegger(afterIm);
        else
            priorImMask = fluoImCellSeg_ys(priorIm);
            afterImMask = fluoImCellSeg_ys(afterIm);
        end
        
        bestPosition = caculateCrossCorrelationForImage(priorImMask,afterImMask,50);%偏移最大[-50,50];ys2020.01.03      
        xyShift_Array(iImage,:) = bestPosition;
        disp([num2str(iImage),' frame']);
        
    end
    cumShift_Array = cumsum(xyShift_Array,1);
    XYtransArray = cat(1,[0,0],cumShift_Array);% 1st frame 不用校正 赋值[0,0]
    
    for iChannel = 1:numel(allChannels)
        dirImage = [dirField,'\',allChannels{iChannel}];
        % 通过try判断有哪些channel
        try
            load([dirImage,'\frameInfo.mat'])
        catch ME
            continue
        end        
        frameInfo_check = frameInfo;
        channelName = allChannels{iChannel};
        disp (['Start',32,channelName,32,'XY shift correciton']);
        
        imageList = dir([dirImage,filesep,'*.tif']);
       
        XYtransArrayNew = newXYtransArrayGetForUnsynFluo(fixedImNum,XYtransArray,frameInfo_fixed,frameInfo_check);
        
        dirImNew = [dirField,filesep,allChannels{iChannel},'_New2'];
        if ~isfolder(dirImNew)
            mkdir(dirImNew)
        end
        
        parfor iImage = 1:numel(imageList)
            imIn = imread([dirField,filesep,channelName,filesep,imageList(iImage).name]);
            imOut = imtranslate(imIn, XYtransArrayNew(iImage,:),'FillValues', 0);
            imwrite(imOut,[dirImNew,filesep,imageList(iImage).name]);
        end        
       
    end
end

end
%%
function newXYtransArray = newXYtransArrayGetForUnsynFluo(fixedImNum,XYtransArray,frameInfo_fixed,frameInfo_check)

if size(frameInfo_fixed,1) == size(frameInfo_check,1)
    newXYtransArray = XYtransArray;
    return;
end

% 寻找tracking的时间窗口内，所拍摄的其他通道图像的数目
%一般Fixedchannel的图像数据点最多

fluo2TrackingIdx = zeros(1, size(frameInfo_check,1));
for iIm = 1 :size (frameInfo_check,1)
    time_lag = abs(etime (frameInfo_check(iIm,1:6),frameInfo_fixed(:,1:6))/60);
    %min 与frameInfo_tracking拍摄的时间间隔绝对值
    trackingIdx = find(time_lag==min(time_lag));
    fluo2TrackingIdx(iIm) = trackingIdx(end);
    %use trckingIdx(end),for the case of trckingIdx 不唯一 ;case:1,1,2,2...
end
fluo2TrackingIdx = fluo2TrackingIdx(fluo2TrackingIdx <= fixedImNum);
[~,checkFrame] = find(diff(fluo2TrackingIdx)==0);
fluo2TrackingIdx(checkFrame+1) = 0;% case for:两张荧光同时对应一张tracking
%赋值为0，数组索引就会出错，便于debug
newXYtransArray = XYtransArray(fluo2TrackingIdx,:);
end


function bestPosition = caculateCrossCorrelationForImage(image1,image2,step)
% for calculater the cross correlation of two image
% step is the searching range

% crop image , for fasr calculaiton
imageStack = cat(3,image1,image2);
[imageStackNew] = cropImageForxyShiftCorrection(imageStack,400);
image1 = imageStackNew(:,:,1);
image2 = imageStackNew(:,:,2);

[x,y] = meshgrid((-step:step)',(-step:step)');
x=x(:);
y=y(:);
correlationMatrix = zeros(size(x,1),1);
parfor i = 1:numel(x)
    se = translate(strel(1),[x(i),y(i)]);
    image2New=imdilate(image2,se);  % 巧用imdilate实现平移
    sumImage=image1 & image2New;    % 利用逻辑矩阵的乘法相当于&
    correlationMatrix(i)=sum(sum(sumImage)); % 直接对逻辑矩阵求和速度比较快
end
bestPosition=[x(correlationMatrix==max(correlationMatrix)),y(correlationMatrix==max(correlationMatrix))];
bestPosition=bestPosition(1,:);
end

function [imageStackNew] = cropImageForxyShiftCorrection(imageStack,pixelRange)
%imageStack为二值图像 2020.01.03ys 寻找图像切割的范围
%   pixelRange = 400;
imageSize = size(imageStack);
if min(imageSize(1:2)) >= 1000
    xCrop = [fix(imageSize(1)/2)-pixelRange,fix(imageSize(1)/2)+pixelRange];
    yCrop = [fix(imageSize(2)/2)-pixelRange,fix(imageSize(2)/2)+pixelRange];
    tempStack = imageStack(xCrop(1):xCrop(2),yCrop(1):yCrop(2),:);
    stats = regionprops(tempStack(:,:,1),'PixelList');
    if numel(stats) < 7
        pixelRange = pixelRange+200;
        if pixelRange >= min(imageSize(1:2))/2
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