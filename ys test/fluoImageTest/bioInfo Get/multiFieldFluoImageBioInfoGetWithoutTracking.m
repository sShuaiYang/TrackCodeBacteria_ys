function [fluo2TrackingIdx] = multiFieldFluoImageBioInfoGetWithoutTracking(dirFile)
%   Shuai Yang 2020.04.04 ����Ӱ�дcode
%���ڶ���Ұ����ӫ��ͼ������ͳ��ӫ��ǿ�� û�о���bioTree����ϸ��׷��
%��montage Ҳ�ɶ���Ұtime-lapse����
% allBioInfo cell�ṹ ��Ŀ������Ұ��Ŀ
%ÿһ��cell����Ϊ�ṹ�� ��¼�˴���Ұ�ﲻͬframe��ϸ����Ϣ

% scale=0.0650;% scale pixel to um

disp('bioInfo get for multifield (time-lapse) fluo image without tracking')
fluoChannels = {'sfGFP','mScarletI','Venus','PVD','CyOFP','TDsmURFP'};
%PhC��BFͼ��һ������mask,���ټ����ǿ����Ϣ
fieldList = dir([dirFile,filesep,'field*']);
% fieldList = fieldListClean (fieldList);%ֻ����field���ļ�������ϵͳ�ļ�
% allBioInfo = cell (1, length(fieldList));
for iField = 1:length(fieldList)

    dirField = strcat(dirFile,'\',fieldList(iField).name);
    if isfolder([dirField,'\Tracking']) %�ж��Ƿ����Tracking �ļ���

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
        % 2022/6/8 ��洢Ϊtifͼ�� ����imread ��ȡ
        trackingList = dir([dirTracking,filesep,'imageTracking*.tif']);
        [fluo2TrackingIdx] = findTrackingImagesForUnsynFluo(dirField);
        if  size(fluo2TrackingIdx,2) <= 500 %�����ٿ�����parfor��������for �Է��ڴ治��
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
% imBioInfo �Ľṹ��fieldBioInfo �ṹ����һ�� single image bioInfo
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
%�����Ȼ�����Ϣ��ȡ
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
    % ��ͬfluo channel ӫ��ǿ�Ȼ�ȡ
    channelList = dir(dirField);
    for iChannel = 1:length(channelList)
        if ismember(channelList(iChannel).name,{'.','..'})||... % ȥ��ϵͳ�Դ����������ļ���
                ~channelList(iChannel).isdir               % ȥ�������в����ļ��е�
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

