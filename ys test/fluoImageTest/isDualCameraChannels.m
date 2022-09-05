function dualCamTF = isDualCameraChannels(dirFile)
% 判断是否是双相机拍摄
% Shuai Yang 2021.07.27
% allChannels={'BF1','PVD','Venus','sfGFP','mScarletI','CyOFP','TDsmURFP'};
channels_cam1 = {'BF1','PVD','sfGFP'};% camera 1拍摄的通道 短波长
channels_cam2 = {'mScarletI','CyOFP','Venus','TDsmURFP'};% camera 2拍摄的通道
warning('on')
%%
fieldList =  dir([dirFile,filesep,'field*']);
dirField_check = [dirFile,filesep,fieldList(1).name];

% 判断是否只是cam 1 拍摄的channel
chTF_cam1 = false (1,numel(channels_cam1));
for iCh_cam1 = 1:numel(channels_cam1)
    chName_cam1 = channels_cam1{iCh_cam1};
    if isfolder([dirField_check,filesep,chName_cam1])
        chTF_cam1(iCh_cam1) = true;
    end
end


% 判断是否有cam2 拍摄的channels
chTF_cam2 = false (1,numel(channels_cam2));
for iCh_cam2 = 1:numel(channels_cam2)
    chName_cam2 = channels_cam2{iCh_cam2};
    if isfolder([dirField_check,filesep,chName_cam2])
        chTF_cam2(iCh_cam2) = true;
    end
end

if sum(chTF_cam1) == 0 ||sum(chTF_cam2) == 0
    dualCamTF = false;
else
    dualCamTF = true;
end

end