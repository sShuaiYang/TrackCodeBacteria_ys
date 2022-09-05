function [fluo2TrackingIdx] = multiFieldFluoImageBioInfoGetWithoutTracking(dirFile)
%   Shuai Yang 2020.04.04 今天加班写code
%用于多视野拍摄荧光图像用于统计荧光强度 没有经过bioTree进行细菌追踪
%可montage 也可多视野time-lapse拍摄
% allBioInfo cell结构 数目代表视野数目
%每一个cell里面为结构体 记录了此视野里不同frame的细菌信息

% scale=0.0650;% scale pixel to um

disp('bioInfo get for multifield (time-lapse) fluo image without tracking')
fluoChannels = {'sfGFP','mScarletI','Venus','PVD','CyOFP','TDsmURFP'};
%PhC的BF图像一般用作mask,不再计算光强等信息
fieldList = dir([dirFile,filesep,'field*']);
% fieldList = fieldListClean (fieldList);%只保留field的文件和两个系统文件
% allBioInfo = cell (1, length(fieldList));
for iField = 1:length(fieldList)

    dirField = strcat(dirFile,'\',fieldList(iField).name);
    if isfolder([dirField,'\Tracking']) %判断是否存在Tracking 文件夹

        % initialize fluo intensity for different channles
        % int represent intensity;
        %fieldBioInfo represnts single field bioInfo
        fieldBioInfo = struct('fieldIdx',[],'frameIdx',[],'intsfGFP',[], ...
            'intmScarletI',[],'intVenus',[],'intPVD',[],'intCyOFP',[], ...
            'intTDsmURFP',[],'BG_sfGFP',[],'BG_mScarletI',[],'BG_Venus',[], ...
            'BG_PVD',[],'BG_CyOFP',[],'BG_TDsmURFP',[],'Centroid',[], ...
            'MajorAxisLength',[],'MinorAxisLength',[],'FilledArea',[], ...
            'PixelIdxList',{},'PixelList',{});

        dirTracking = strcat(dirField,'\Tracking');
        % 2022/6/8 后存储为tif图像 采用imread 读取
        trackingList = dir([dirTracking,filesep,'imageTracking*.tif']);
        [fluo2TrackingIdx] = findTrackingImagesForUnsynFluo(dirField);
        if  size(fluo2TrackingIdx,2) <= 500 %数据少可以用parfor，否则用for 以防内存不足
            parfor iImage = 1:size(fluo2TrackingIdx,2)
                maskImName = [dirTracking,'\',trackingList(fluo2TrackingIdx(iImage)).name];
                %     maskImage = load(maskImName);
                %     maskImage = maskImage.imageTracking;
                maskImage = imread(maskImName);
                maskImage = logical(maskImage);
                [imBioInfo] = multiChannelOneImageFluoIntensityGet (fluoChannels,maskImage,dirField,iImage);
                fieldBioInfo (iImage) = imBioInfo;
            end
        else
            for iImage = 1:size(fluo2TrackingIdx,2)
                maskImName = [dirTracking,'\',trackingList(fluo2TrackingIdx(iImage)).name];
                %     maskImage = load(maskImName);
                %     maskImage = maskImage.imageTracking;
                maskImage = imread(maskImName);
                maskImage = logical(maskImage);                
                [imBioInfo] = multiChannelOneImageFluoIntensityGet (fluoChannels,maskImage,dirField,iImage);
                fieldBioInfo (iImage) = imBioInfo;
            end
        end
        bioInfo = fieldBioInfo;
        save(strcat(dirField,'\bioInfo.mat'),'bioInfo');
    else
        disp ('Tracking folder does not exist')
        return
    end

    %     allBioInfo{iField} = fieldBioInfo;
    %     save(strcat(dirFile,'\allBioInfo.mat'),'allBioInfo');
end
end
%%
function [imBioInfo] = multiChannelOneImageFluoIntensityGet (fluoChannels,maskImage,dirField,iImage)
% imBioInfo 的结构与fieldBioInfo 结构保持一致 single image bioInfo
imBioInfo.fieldIdx = [];
imBioInfo.frameIdx = [];
imBioInfo.intsfGFP=[];imBioInfo.intmScarletI=[];
imBioInfo.intVenus=[];imBioInfo.intPVD=[];
imBioInfo.intCyOFP=[];imBioInfo.intTDsmURFP=[];
imBioInfo.BG_sfGFP=[];imBioInfo.BG_mScarletI=[];
imBioInfo.BG_Venus=[];imBioInfo.BG_PVD=[];
imBioInfo.BG_CyOFP=[];imBioInfo.BG_TDsmURFP=[];
imBioInfo.Centroid=[];imBioInfo.MajorAxisLength=[];imBioInfo.MinorAxisLength=[];
imBioInfo.FilledArea=[];imBioInfo.PixelIdxList=[];
imBioInfo.PixelList=[];

imBioInfo.fieldIdx = str2double(dirField(end-3:end));
imBioInfo.frameIdx = iImage;
%菌长等基本信息获取
stats = regionprops(maskImage,'Centroid','MajorAxisLength', ...
    'MinorAxisLength','FilledArea','PixelIdxList','PixelList');
if numel(stats) > 0
    for iCell = 1 : numel(stats)
        imBioInfo.Centroid(iCell,:) = stats(iCell).Centroid;
        imBioInfo.MajorAxisLength(iCell,:) = stats(iCell).MajorAxisLength;
        imBioInfo.MinorAxisLength(iCell,:) = stats(iCell).MinorAxisLength;
        imBioInfo.FilledArea(iCell,:) = stats(iCell).FilledArea;
        imBioInfo.PixelIdxList{iCell} = stats(iCell).PixelIdxList;
        imBioInfo.PixelList{iCell} = stats(iCell).PixelList;
    end
    % 不同fluo channel 荧光强度获取
    channelList = dir(dirField);
    for iChannel = 1:length(channelList)
        if ismember(channelList(iChannel).name,{'.','..'})||... % 去除系统自带的两个隐文件夹
                ~channelList(iChannel).isdir               % 去除遍历中不是文件夹的
            continue;
        end
        
        if ismember(channelList(iChannel).name,fluoChannels )
            fluoChannel = channelList(iChannel).name;
            dirImage = strcat(dirField,'\',channelList(iChannel).name);
            imageList = dir([dirImage,filesep,'image*.tif']);
            imFluo = import_tiff_stack( strcat(dirImage,'\',imageList(iImage).name) );
            ccFluo = regionprops(maskImage,imFluo,'MeanIntensity');
            [~,BG] = substractBackGround(imFluo);

            switch fluoChannel
                case 'sfGFP'
                    imBioInfo.intsfGFP = [ccFluo.MeanIntensity]';
                    imBioInfo.BG_sfGFP = BG;
                case 'mScarletI'
                    imBioInfo.intmScarletI = [ccFluo.MeanIntensity]';
                    imBioInfo.BG_mScarletI = BG;
                case 'Venus'
                    imBioInfo.intVenus = [ccFluo.MeanIntensity]';
                    imBioInfo.BG_Venus = BG;
                case 'PVD'
                    imBioInfo.intPVD = [ccFluo.MeanIntensity]';
                    imBioInfo.BG_PVD = BG;
                case 'CyOFP'
                    imBioInfo.intCyOFP = [ccFluo.MeanIntensity]';
                    imBioInfo.BG_CyOFP = BG;
                case 'TDsmURFP'
                    imBioInfo.intTDsmURFP = [ccFluo.MeanIntensity]';
                    imBioInfo.BG_TDsmURFP = BG;
            end

        end
    end
end
end

