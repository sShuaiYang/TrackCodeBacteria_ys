function addingBKGSignal2BioInfo(dirFile)
% Shuai Yang 2021/10/29
% 需要的数据存储结构 dirFile{field001,field0002,......}

disp('adding backGround Signal to bioInfo')
fluoChannels={'sfGFP','mScarletI','Venus','PVD','CyOFP','TDsmURFP'};
fieldList = dir([dirFile,filesep,'field*']);
for iField= 1:(length(fieldList))
    t0 = clock;

    disp(fieldList(iField).name);
    
    dirField = strcat(dirFile,'\',fieldList(iField).name);
    load([dirField,'\bioInfo.mat']);

    parfor iInfo = 1:numel(bioInfo)
        image = bioInfo(iInfo).frameIdx;
        [imBKGInfo] = multiChannelOneImBKGSignalGet (fluoChannels,dirField,image);
        bioInfo(iInfo).BG_sfGFP = imBKGInfo.BG_sfGFP;
        bioInfo(iInfo).BG_mScarletI = imBKGInfo.BG_mScarletI;
        bioInfo(iInfo).BG_Venus = imBKGInfo.BG_Venus;
        bioInfo(iInfo).BG_PVD = imBKGInfo.BG_PVD;
        bioInfo(iInfo).BG_CyOFP = imBKGInfo.BG_CyOFP;
        bioInfo(iInfo).BG_TDsmURFP = imBKGInfo.BG_TDsmURFP;
    end

    save (strcat(dirField,'\bioInfo.mat'),'bioInfo','-v7.3');
    clear bioInfo
    t1 = clock;
    disp(['共耗时:',num2str(etime(t1,t0)),'秒']);
end
end

function [imBKGInfo] = multiChannelOneImBKGSignalGet (fluoChannels,dirField,iImage)

imBKGInfo = struct('BG_sfGFP',[], ...
    'BG_mScarletI',[],'BG_Venus',[],'BG_PVD',[], ...
    'BG_CyOFP',[],'BG_TDsmURFP',[]);

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