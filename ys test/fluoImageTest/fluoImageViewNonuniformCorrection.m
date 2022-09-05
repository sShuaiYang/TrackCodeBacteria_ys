function fluoImageViewNonuniformCorrection(dirFile)
%视场不均匀校正，需要读取mip中的荧光背景信息，实验需要操作拍摄
% shuai Yang 2020.05.26
disp('FluoImage view non-uniformity correction');
load([dirFile,'\mip.mat']);
fluoChannels={'sfGFP','mScarletI','Venus','PVD','CyOFP','TDsmURFP'};
% fluoChannels 需要与backGround对应，BG is short for backGround
BG_sfGFP = mip.Calib.LumencorBlueField;
BG_mScarletI = mip.Calib.LumencorTealField;
BG_Venus = mip.Calib.LumencorVioletField;
BG_PVD = mip.Calib.LumencorCyanField;
BG_CyOFP = mip.Calib.LumencorCyanField;
BG_TDsmURFP = mip.Calib.LumencorGreenField;

fieldList = dir([dirFile,filesep,'field*']);

for iField = 1:length(fieldList)
    if ~strcmp(fieldList(iField).name(1:5),'field')
        continue
    end
    dirField = strcat(dirFile,'\',fieldList(iField).name);
    channelList = dir(dirField);
    for iChannel = 1:length(channelList)
        if(isequal(channelList(iChannel).name,'.')||... % 去除系统自带的两个隐文件夹
                isequal(channelList(iChannel).name,'..')||...
                ~channelList(iChannel).isdir)                % 去除遍历中不是文件夹的
            continue;
        end
        if ~ismember(channelList(iChannel).name,fluoChannels)
            continue
        end
        fluoChannel = channelList(iChannel).name;
        dirImage = strcat(dirField,'\',channelList(iChannel).name);
        imageList = dir([dirImage,filesep,'image*']);
        for iImage = 1:length(imageList)
            if ~strcmp(imageList(iImage).name(1:5),'image')
                continue
            end
            imageFluo = import_tiff_stack(strcat(dirImage,'\',imageList(iImage).name));
            imageFluo = double(imageFluo);
            switch fluoChannel
                case 'sfGFP'
                    imageFluo = imageFluo./BG_sfGFP;
                case 'mScarletI'
                    imageFluo = imageFluo./BG_mScarletI;
                case 'Venus'
                    imageFluo = imageFluo./BG_Venus;
                case 'PVD'
                    imageFluo = imageFluo./BG_PVD;
                case 'CyOFP'
                    imageFluo = imageFluo./BG_CyOFP;
                case 'TDsmURFP'
                    imageFluo = imageFluo./BG_TDsmURFP;
            end
            imageFluo = uint16(imageFluo);
            imwrite(imageFluo,strcat(dirImage,'\',imageList(iImage).name));
        end
        
    end
end
end

%%
