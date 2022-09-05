function [XYtransArray] = batchXYShiftCorrectionUsingFTPhaseCorrelation(dirFile,fixedChannel)
% time lapse 拍摄 荧光图像会出现xy的漂移 通过计算fixedChannel的偏移量 进行校正
%认为其他通道的xy偏移量与fixedChannel一致
% Shuai Yang 2021.06.30 采用 Fourier transform 空间的FFT upsampled cross
% correlation进行计算偏移量
%Reference
% 参考 bactch multiField channel register test1.mlx 参考文献 Manuel
% Guizar-Sicairos, Samuel T. Thurman, and James R. Fienup,"Efficient
% subpixel image registration algorithms," Opt. Lett. 33, 156-158 (2008).
% 核心函数dftregistration
% 不适合 视野细菌数目只有几个的图像
tic;

disp('XY shift correction using Fourier Phase Correlation ')
allChannels = {'BF1','sfGFP','mScarletI','Venus','PVD','CyOFP','TDsmURFP'};
% fluoChannels = {'sfGFP','mScarletI','Venus','PVD','CyOFP','TDsmURFP'};

if ismember(fixedChannel,allChannels)
    disp([fixedChannel, 32, 'used as the fixed channel'])
else
    disp ( 'input error of fixedChannel')
    return
end

fieldList = dir([dirFile,filesep,'field*']);
usfac = 1;

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
    
    outPut_Array = zeros(numel(imList)-1,4);
    parfor iImage = 1:numel(priorImList)
        priorIm = imread([dirField,filesep,fixedChannel,filesep,priorImList(iImage).name]);
        afterIm = imread([dirField,filesep,fixedChannel,filesep,afterImList(iImage).name]);
        
        %         % 细菌太少的话 计算出mask图像 再做相关 效果好
        %         priorIm = fluoImCellMask_ys(priorIm);
        %         afterIm = fluoImCellMask_ys(afterIm);
        
        %         priorIm = im2double(priorIm);
        %         afterIm = im2double(afterIm);
        % fourier transform images
        fft_priorIm = fft2(priorIm);
        fft_afterIm = fft2(afterIm);
        
        [output, ~] = dftregistration(fft_priorIm,fft_afterIm,usfac);
        outPut_Array(iImage,:) = output;
    end
    xyshiftPixels = outPut_Array(:,end-1:end);%第3和4列是xy方向的平移pixel大小
    cumXYshift = cumsum(xyshiftPixels,1);%Cumulative 累积的xy偏移量
    cumShift_Array(:,1:2) = outPut_Array(:,1:2);%累计的Pixel shifts between images
    cumShift_Array(:,3:4) = cumXYshift;
    XYtransArray = cumXYshift(:,[2,1]);
    XYtransArray = cat(1,[0,0],XYtransArray);% 1st frame 不用校正 赋值[0,0]

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

%         dirImNew = [dirField,filesep,allChannels{iChannel},'_FT'];
%         if ~isfolder(dirImNew)
%             mkdir(dirImNew)
%         end

        parfor iImage = 1:numel(imageList)
            imIn = imread([dirImage,filesep,imageList(iImage).name]);
            imOut = imtranslate(imIn, XYtransArrayNew(iImage,:),'FillValues', 0);
            imwrite(imOut,[dirImage,filesep,imageList(iImage).name]);
            %   imwrite(imOut,[dirImNew,filesep,imageList(iImage).name]);
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