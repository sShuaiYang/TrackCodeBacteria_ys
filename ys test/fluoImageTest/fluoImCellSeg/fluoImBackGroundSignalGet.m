function fluoImBackGroundSignalGet(dirFile)
% 用于获取实验中荧光图像的背景 fluo image backGround get
% Shuai Yang 2021/10/28
% 目前只适用于不同荧通道 光图像数目相等的情况
disp('get fluo image backGround signal')
fluoChannels={'sfGFP','mScarletI','Venus','PVD','CyOFP','TDsmURFP'};

fieldList = dir([dirFile,filesep,'field*']);

for iField  = 1:length(fieldList)

    disp(fieldList(iField).name);

    % BKG is short for backGround
    BKGInfo = struct('fieldIdx',[],'frameIdx',[],'BG_sfGFP',[], ...
        'BG_mScarletI',[],'BG_Venus',[],'BG_PVD',[], ...
        'BG_CyOFP',[],'BG_TDsmURFP',[]);
    dirField = strcat(dirFile,'\',fieldList(iField).name);

    channelList = dir(dirField);
    fileName = {channelList.name};
    Lia = ismember(fluoChannels,fileName);
    k  = find(Lia);
    fluoImCheck = zeros(size(k)); % 检验不同通道图像数目是否相等
    for iCh = 1:numel(k)
        channelName = fluoChannels{k(iCh)};
        dirImage = [dirField,filesep,channelName];

        imageList = dir([dirImage,filesep,'image*.tif']);
        imageNum = length(imageList);
        fluoImCheck(iCh) = imageNum;
    end
    imNum = unique(fluoImCheck);
    if numel(imNum) ~= 1 % =!1 说明数目不相等
        msg = [fieldList(iField).name,32,'Abnormal !!!'];
        disp(msg)
        warning(msg)
        continue
    end

    parfor iImage = 1:imNum
        [imBKGInfo] = multiChannelOneImBKGSignalGet (fluoChannels,dirField,iImage);
        BKGInfo(iImage) = imBKGInfo;
    end
    save(strcat(dirField,'\BKGInfo.mat'),'BKGInfo');

end
end

function [imBKGInfo] = multiChannelOneImBKGSignalGet (fluoChannels,dirField,iImage)

imBKGInfo = struct('fieldIdx',[],'frameIdx',[],'BG_sfGFP',[], ...
    'BG_mScarletI',[],'BG_Venus',[],'BG_PVD',[], ...
    'BG_CyOFP',[],'BG_TDsmURFP',[]);

imBKGInfo.fieldIdx = str2double(dirField(end-3:end));
imBKGInfo.frameIdx = iImage;

channelList = dir(dirField);
fileName = {channelList.name};
Lia = ismember(fluoChannels,fileName);
k  = find(Lia);
for iCh = 1:numel(k)
    fluoChannel = fluoChannels{k(iCh)};
    dirImage = strcat(dirField,'\',fluoChannel);
    imageList = dir([dirImage,filesep,'image*.tif']);
    imFluo = imread([dirImage,filesep,imageList(iImage).name]);
    [~,BG] = substractBackGround(imFluo);
    switch fluoChannel
        case 'sfGFP'
            imBKGInfo.BG_sfGFP = BG;
        case 'mScarletI'
            imBKGInfo.BG_mScarletI = BG;
        case 'Venus'
            imBKGInfo.BG_Venus = BG;
        case 'PVD'
            imBKGInfo.BG_PVD = BG;
        case 'CyOFP'
            imBKGInfo.BG_CyOFP = BG;
        case 'TDsmURFP'
            imBKGInfo.BG_TDsmURFP = BG;
    end

end

end