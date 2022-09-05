function [fluo2TrackingIdx] = fluoColonyImageBioInfoGet(dirFile)
%% 新注释
%通过一下code修改得到 2021.03.29 Shuai Yang
%主要在计算光强的时候，针对每张荧光图像减除了背景，因为体式显微镜拍摄的图像背景差别很大
%[fluo2TrackingIdx] = multiFieldFluoImageBioInfoGetWithoutTracking(dirFile)
%加入修改 imFluo = substractBackGround(imFluo);
%% 原始注释
%    2020.04.04 今天加班写code
%用于多视野拍摄荧光图像用于统计荧光强度 没有经过bioTree进行细菌追踪
%可montage 也可多视野time-lapse拍摄
% allBioInfo cell结构 数目代表视野数目
%每一个cell里面为结构体 记录了此视野里不同frame的细菌信息
%%
% scale=0.0650;% scale pixel to um

disp('bioInfo get for multifield (time-lapse) fluo image without tracking')
fluoChannels={'sfGFP','mScarletI','Venus','PVD','CyOFP','TDsmURFP'};
%PhC的BF图像一般用作mask,不再计算光强等信息
fieldList = dir([dirFile,filesep,'field*']);
fieldList = fieldListClean (fieldList);%只保留field的文件和两个系统文件
% allBioInfo = cell (1, length(fieldList));
for iField = 1:length(fieldList)
    if strcmp(fieldList(iField).name(1:5),'field')
        dirField = strcat(dirFile,'\',fieldList(iField).name);
        if isfolder([dirField,'\Tracking']) %判断是否存在Tracking 文件夹
            
            % initialize fluo intensity for different channles
            % int represent intensity;
            %fieldBioInfo represnts single field bioInfo
            fieldBioInfo = struct('fieldIdx',[],'frameIdx',[],'intsfGFP',[ ],'intmScarletI',[ ],'intVenus',[ ],...
                'intPVD',[ ],'intCyOFP',[ ],'intTDsmURFP',[ ],'Centroid',[ ],'MajorAxisLength',[ ],...
                'MinorAxisLength',[ ],'FilledArea',[ ],'PixelIdxList',{ },'PixelList',{ },...
                'pixelsfGFP',[ ],'pixelmScarletI',[ ],'pixelVenus',[ ],'pixelPVD',[ ],...
                'pixelCyOFP',[ ],'pixelTDsmURFP',[ ]);
            
            dirTracking = strcat(dirField,'\Tracking');
            trackingList = dir(dirTracking);
            [fluo2TrackingIdx] = findTrackingImagesForUnsynFluo(dirField);
            if  size(fluo2TrackingIdx,2) <= 500 %数据少可以用parfor，否则用for 以防内存不足
                parfor iImage = 1 :  size(fluo2TrackingIdx,2)
                    maskImage = load(strcat(dirTracking,'\',trackingList(fluo2TrackingIdx(iImage)+3).name));
                    maskImage = maskImage.imageTracking;
                    [imBioInfo] = multiChannelOneImageFluoIntensityGet (fluoChannels,maskImage,dirField,iImage);
                    fieldBioInfo (iImage) = imBioInfo;
                end
            else
                for iImage = 1 :  size(fluo2TrackingIdx,2)
                    maskImage = load(strcat(dirTracking,'\',trackingList(fluo2TrackingIdx(iImage)+3).name));
                    maskImage = maskImage.imageTracking;
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
imBioInfo.intVenus=[];imBioInfo.intPVD=[];imBioInfo.intCyOFP=[];imBioInfo.intTDsmURFP=[];
imBioInfo.Centroid=[];imBioInfo.MajorAxisLength=[];imBioInfo.MinorAxisLength=[];
imBioInfo.FilledArea=[];imBioInfo.PixelIdxList=[];
imBioInfo.PixelList=[];
imBioInfo.pixelsfGFP=[];imBioInfo.pixelmScarletI=[];
imBioInfo.pixelVenus=[];imBioInfo.pixelPVD=[];
imBioInfo.pixelCyOFP=[];imBioInfo.pixelTDsmURFP=[];

imBioInfo.fieldIdx = str2double(dirField(end-3:end));
imBioInfo.frameIdx = iImage;
%菌长等基本信息获取
stats=regionprops(maskImage,'Centroid','MajorAxisLength','MinorAxisLength','FilledArea','PixelIdxList','PixelList');
if numel(stats) > 0
    for iCell = 1 : numel(stats)
        imBioInfo. Centroid(iCell,:) = stats(iCell).Centroid;
        imBioInfo. MajorAxisLength(iCell,:) = stats(iCell).MajorAxisLength;
        imBioInfo.MinorAxisLength(iCell,:) = stats(iCell).MinorAxisLength;
        imBioInfo. FilledArea(iCell,:) = stats(iCell).FilledArea;
        imBioInfo. PixelIdxList{iCell} = stats(iCell).PixelIdxList;
        imBioInfo. PixelList{iCell} = stats(iCell).PixelList;
    end
    % 不同fluo channel 荧光强度获取
    channelList = dir(dirField);
    for iChannel = 1:length(channelList)
        if(isequal(channelList(iChannel).name,'.')||... % 去除系统自带的两个隐文件夹
                isequal(channelList(iChannel).name,'..')||...
                ~channelList(iChannel).isdir)                % 去除遍历中不是文件夹的
            continue;
        end
        
        if ismember(channelList(iChannel).name,fluoChannels )
            fluoChannel = channelList(iChannel).name;
            dirImage = strcat(dirField,'\',channelList(iChannel).name);
            imageList = dir([dirImage,filesep,'image*.tif']);
            imFluo = import_tiff_stack( strcat(dirImage,'\',imageList(iImage).name) );
            imFluo = substractBackGround(imFluo);
            ccFluo = regionprops(maskImage,imFluo,'MeanIntensity','PixelValues');
            %针对cellColony image 加入了'PixelValues'， Shuai Yang 2021.04.14
            
            switch fluoChannel
                case 'sfGFP'
                    for iCell = 1 : numel(ccFluo)
                        imBioInfo. intsfGFP(iCell,:) = ccFluo(iCell).MeanIntensity;
                        imBioInfo. pixelsfGFP{iCell} = ccFluo(iCell).PixelValues;
                    end
                case 'mScarletI'
                    for iCell = 1 : numel(ccFluo)
                        imBioInfo. intmScarletI(iCell,:) = ccFluo(iCell).MeanIntensity;
                        imBioInfo. pixelmScarletI{iCell} = ccFluo(iCell).PixelValues;
                    end
                case 'Venus'
                    for iCell = 1 : numel(ccFluo)
                        imBioInfo. intVenus(iCell,:) = ccFluo(iCell).MeanIntensity;
                        imBioInfo. pixelVenus{iCell} = ccFluo(iCell).PixelValues;
                    end
                case 'PVD'
                    for iCell = 1 : numel(ccFluo)
                        imBioInfo. intPVD(iCell,:) = ccFluo(iCell).MeanIntensity;
                        imBioInfo. pixelPVD{iCell} = ccFluo(iCell).PixelValues;
                    end
                case 'CyOFP'
                    for iCell = 1 : numel(ccFluo)
                        imBioInfo. intCyOFP(iCell,:) = ccFluo(iCell).MeanIntensity;
                        imBioInfo. pixelCyOFP{iCell} = ccFluo(iCell).PixelValues;
                    end
                case 'TDsmURFP'
                    for iCell = 1 : numel(ccFluo)
                        imBioInfo. intTDsmURFP(iCell,:) = ccFluo(iCell).MeanIntensity;
                        imBioInfo. pixelTDsmURFP{iCell} = ccFluo(iCell).PixelValues;
                    end
            end
            
        end
    end
end
end
%%
function fieldList = fieldListClean (fieldList)
% 只保留field 文件和两个系统文件
templogic = false (1, numel (fieldList));
for iField = 1: numel (fieldList)
    
    if(isequal(fieldList(iField).name,'.')||... % 系统自带的两个隐文件夹
            isequal(fieldList(iField).name,'..')||...
            strcmp(fieldList(iField).name(1:5),'field'))
        templogic(iField) = true;
    end
    
end
fieldList = fieldList(templogic);
end
%%
function I = substractBackGround(I)
%用于减除图像的背景 Shuai Yang 2020.09.09
sz = size(I);
I0 = I(round(sz(1)/2-sz(1)/4):round(sz(1)/2+sz(1)/4),round(sz(2)/2-sz(2)/4):round(sz(2)/2+sz(2)/4));
I0 = double(I0);
I0 = sort(I0(:));
pixelSpace = 1:100:numel(I0);%每100个pixel计算平均值
M = zeros(numel(pixelSpace)-1,1);
for i = 1:numel(pixelSpace)-1
    M(i) = mean(I0(pixelSpace(i):pixelSpace(i+1)));
end

n = 5;
CV = std(M(1:n))/mean(M(1:n));
while CV > 0.03 % equal CV > 0.03 && n >= 1
    n = n-1;
    CV = std(M(1:n))/mean(M(1:n));
end
BG = mean(M(1:n));
I = I - BG;
end
