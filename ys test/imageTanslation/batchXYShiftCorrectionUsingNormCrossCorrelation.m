function [XYtransArray] = batchXYShiftCorrectionUsingNormCrossCorrelation(dirFile,fixedChannel)
% time lapse 拍摄 荧光图像会出现xy的漂移 通过计算fixedChannel的偏移量 进行校正
%认为其他通道的xy偏移量与fixedChannel一致
% Shuai Yang 2021.06.30 采用 NormalizedCrossCorrelation进行偏移量计算
%相比FT Phase correlation 计算速度慢些
%Reference
% ImageRegistrationUsingNormalizedCrossCorrelation.mlx
%对于视野中只有细菌数目较少的图像，适合用mask来做相关
tic;
disp('XY shift correction using Normalized Cross Correlation')
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
    parfor iImage = 1:numel(priorImList)
        priorIm = imread([dirField,filesep,fixedChannel,filesep,priorImList(iImage).name]);
        afterIm = imread([dirField,filesep,fixedChannel,filesep,afterImList(iImage).name]);
        
        %         % 细菌太少的话 计算出mask图像 再做相关 效果好
        %         priorIm = fluoImCellMask_ys(priorIm);
        %         afterIm = fluoImCellMask_ys(afterIm);
        
        c = normxcorr2(priorIm,afterIm);
        
        %         [ypeak,xpeak] = find(c == max(c(:)));
        %         yShift = size(priorIm,1) - ypeak;
        %         xShift = size(priorIm,2) - xpeak;
        %         xyShift_Array(iImage,:) = [xShift,yShift];
        [~, imax] = max(abs(c(:)));
        [ypeak, xpeak] = ind2sub(size(c),imax(1));
        corr_offset = [(size(priorIm,2)-xpeak)
            (size(priorIm,1)-ypeak)];
        x_offset = corr_offset(1);
        y_offset = corr_offset(2);
        xyShift_Array(iImage,:) = [x_offset,y_offset];
       
        
        
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
        
        dirImNew = [dirField,filesep,allChannels{iChannel},'_NCC'];
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
toc;
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