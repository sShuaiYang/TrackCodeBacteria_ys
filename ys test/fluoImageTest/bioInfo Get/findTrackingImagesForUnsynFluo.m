function [fluo2TrackingIdx] = findTrackingImagesForUnsynFluo(dirField)
% UnsynFluo unsynchronized fluorescent shotting
%case 1 针对PhC图像(Tracking)拍摄的时间间隔和荧光图像时间间隔不同步
%case 2：生成Tracking文件时，截取了前面一定的时间区域，并非全部

fluoChannels = {'sfGFP','mScarletI','Venus','PVD','CyOFP','TDsmURFP'};
dirTracking = strcat(dirField,'\Tracking');
trackingImList = dir([dirTracking,filesep,'imageTracking*']);
trackingNum = length(trackingImList);%tracking文件的实际数目
load([dirTracking,'\frameInfo.mat']);
frameInfo_tracking = frameInfo;
channelList = dir(dirField);
for iChannel = 1:length(channelList)
    if(ismember(channelList(iChannel).name,{'.','..'})||... % 去除系统自带的两个隐文件夹
            ~channelList(iChannel).isdir)                % 去除遍历中不是文件夹的
        continue;
    end
    
    if ismember(channelList(iChannel).name,fluoChannels )
        dirFluoImage = strcat(dirField,'\',channelList(iChannel).name);
        load([dirFluoImage,'\frameInfo.mat']);
        frameInfo_fluo = frameInfo;
        break;
    end
end

%如果frameInfo的长度相等 直接返回
if length (frameInfo_fluo) == length (frameInfo_tracking)
    fluo2TrackingIdx = 1:trackingNum;
    return
end

% 寻找tracking的时间窗口内，所拍摄的荧光图像的数目
fluo2TrackingIdx = zeros (1,size (frameInfo_fluo,1));
for iFluo = 1 :size (frameInfo_fluo,1)
    time_lag = abs(etime (frameInfo_fluo(iFluo,1:6),frameInfo_tracking(:,1:6))/60);
    %min 与frameInfo_tracking拍摄的时间间隔绝对值
    trckingIdx = find(time_lag == min(time_lag));
    fluo2TrackingIdx(iFluo) = trckingIdx(end);
    %use trckingIdx(end),for the case of trckingIdx 不唯一 ;case:1,1,2,2... 
end
fluo2TrackingIdx = fluo2TrackingIdx(fluo2TrackingIdx <= trackingNum);

end
