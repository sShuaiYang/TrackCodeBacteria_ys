function singleFluoChTF = isSingleFluoChannel(dirFile)
% 判断是否只有一个荧光通道
% Shuai Yang 2021.07.27
% allChannels={'BF1','PVD','Venus','sfGFP','mScarletI','CyOFP','TDsmURFP'};
% channels_cam1 = {'BF1','PVD','sfGFP'};% camera 1拍摄的通道 短波长
% channels_cam2 = {'mScarletI','CyOFP','Venus','TDsmURFP'};% camera 2拍摄的通道

fluoChannels={'PVD','Venus','sfGFP','mScarletI','CyOFP','TDsmURFP'};
warning('on')
%%
fieldList =  dir([dirFile,filesep,'field*']);
dirField_check = [dirFile,filesep,fieldList(1).name];

chTF_fluo = false (1,numel(fluoChannels));
for iCh_fluo = 1:numel(fluoChannels)
    chName_fluo = fluoChannels{iCh_fluo};
    if isfolder([dirField_check,filesep,chName_fluo])
        chTF_fluo(iCh_fluo) = true;
    end
end

if sum(chTF_fluo) == 1
    singleFluoChTF = true;
else
    singleFluoChTF = false;
end

end