function [bioInfo] = oneShotMontageFluoImageBioInfoGet(dirFile)
%   Shuai Yang 2020.04.04 ����Ӱ�дcode
%����montage ÿ����Ұֻ����һ�εķ���
%bioInfo struct�ṹ
%�뺯��multiFieldFluoImageBioInfoGet����������
%bioInfo Ϊ�ṹ�� ��ÿһ�д�����һ����Ұϸ������Ϣ
%Ŀ���Ƕ���Ұ����parfor����

% dirFile='F:\2020-04-03-fluocode test ys';
% scale=0.0650;% scale pixel to um

disp('get bioInfo for oneshot montage fluo image')
fluoChannels={'sfGFP','mScarletI','Venus','PVD','CyOFP','TDsmURFP'};
bioInfo = struct('fieldIdx',[],'frameIdx',[],'intsfGFP',[], ...
    'intmScarletI',[],'intVenus',[],'intPVD',[],'intCyOFP',[], ...
    'intTDsmURFP',[],'BG_sfGFP',[],'BG_mScarletI',[],'BG_Venus',[], ...
    'BG_PVD',[],'BG_CyOFP',[],'BG_TDsmURFP',[],'Centroid',[], ...
    'MajorAxisLength',[],'MinorAxisLength',[],'FilledArea',[], ...
    'PixelIdxList',{},'PixelList',{});

% BG is short for backGround signal

% fieldList = dir(dirFile);
% fieldList = fieldListClean (fieldList);%ֻ����field���ļ�������ϵͳ�ļ�
% dirFieldOne = strcat(dirFile,'\','field0001');
fieldList = dir([dirFile,filesep,'field*']);
dirFieldOne = strcat(dirFile,'\',fieldList(1).name);%case for field0001 ����Ϊɾ��
% ȡ��һ����Ұ�ж��Ƿ����Tracking �ļ���
if isfolder([dirFieldOne,'\Tracking'])
    parfor iField = 1 : length (fieldList) %
        dirField = strcat(dirFile,'\',fieldList(iField).name);
        dirMask = strcat(dirField,'\Tracking');

        %         maskList = dir([dirMask,filesep,'imageTracking*.mat']);
        %         maskImage = load(strcat(dirMask,'\',maskList(1).name));
        %         % ֻload��һ��trackingͼ��
        %         [~,name,~] = fileparts(strcat(dirMask,'\',maskList(1).name));
        %         maskImage = maskImage.imageTracking;
        
        % 2022/5/25 ��洢Ϊtifͼ�� ����imread ��ȡ
        maskList = dir([dirMask,filesep,'imageTracking*.tif']);
        maskImage = imread([dirMask,filesep,maskList(1).name]);
        % ֻimread��һ��trackingͼ��
        [~,name,~] = fileparts(strcat(dirMask,'\',maskList(1).name));
        maskImage = logical(maskImage);
        iImage = str2double(name(end-4:end));%��ȡ������ͼ������mask

        [imBioInfo] = multiChannelOneImageFluoIntensityGet(fluoChannels,maskImage,dirField,iImage);
        bioInfo (iField) = imBioInfo;
    end
    return
elseif isfolder([dirFieldOne,'\segmentation']) % Delta segmentation % Shuai Yang 2022/5/9
    parfor iField = 1 : length (fieldList) %
        dirField = strcat(dirFile,'\',fieldList(iField).name);
        dirMask = strcat(dirField,'\segmentation');
        maskList = dir([dirMask,filesep,'*.tif']);
        maskImage = imread([dirMask,'\',maskList(1).name]);
        % ֻload��һ��trackingͼ��
        [~,name,~] = fileparts(strcat(dirMask,'\',maskList(1).name));
        maskImage = logical(maskImage);
        iImage = str2double(name(end-4:end));%��ȡ������ͼ������mask
        [imBioInfo] = multiChannelOneImageFluoIntensityGet(fluoChannels,maskImage,dirField,iImage);
        bioInfo (iField) = imBioInfo;
    end
    return
else
    disp ('Tracking  or segmentation folder does not exist')
    return
end
% save(strcat(dirFile,'\bioInfo.mat'),'bioInfo','-v7.3');
end

function [imBioInfo] = multiChannelOneImageFluoIntensityGet(fluoChannels,maskImage,dirField,iImage)

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

% imBioInfo.fieldIdx = iField;
imBioInfo.fieldIdx = str2double(dirField(end-3:end)); % case for ��������field
imBioInfo.frameIdx = iImage;
%�����Ȼ�����Ϣ��ȡ
stats = regionprops(maskImage,'Centroid','MajorAxisLength','MinorAxisLength','FilledArea','PixelIdxList','PixelList');
if numel(stats) > 0
    for iCell = 1 : numel(stats)
        imBioInfo. Centroid(iCell,:) = stats(iCell).Centroid;
        imBioInfo. MajorAxisLength(iCell,:) = stats(iCell).MajorAxisLength;
        imBioInfo. MinorAxisLength(iCell,:) = stats(iCell).MinorAxisLength;
        imBioInfo. FilledArea(iCell,:) = stats(iCell).FilledArea;
        imBioInfo. PixelIdxList{iCell} = stats(iCell).PixelIdxList;
        imBioInfo. PixelList{iCell} = stats(iCell).PixelList;
    end
    % ��ͬfluo channel ӫ��ǿ�Ȼ�ȡ
    channelList = dir(dirField);
    for iChannel = 1:length(channelList)
        if(ismember(channelList(iChannel).name,{'.','..'})||... % ȥ��ϵͳ�Դ����������ļ���
                ~channelList(iChannel).isdir)                % ȥ�������в����ļ��е�
            continue;
        end

        if ismember(channelList(iChannel).name,fluoChannels )
            fluoChannel = channelList(iChannel).name;
            dirImage = strcat(dirField,'\',channelList(iChannel).name);
            imageName = ['image',fluoChannel,num2str(iImage,'%05d'),'.tif'];
            imFluo = import_tiff_stack( strcat(dirImage,'\',imageName));
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

function fieldList = fieldListClean (fieldList)
% ֻ����field �ļ�������ϵͳ�ļ�
templogic = false (1, numel (fieldList));
for iField = 1: numel (fieldList)

    if(isequal(fieldList(iField).name,'.')||... % ϵͳ�Դ����������ļ���
            isequal(fieldList(iField).name,'..')||...
            strcmp(fieldList(iField).name(1:5),'field'))
        templogic(iField) = true;
    end

end
fieldList = fieldList(templogic);
end

