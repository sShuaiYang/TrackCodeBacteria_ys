function [bioInfo] = cellFluoIntensityCorrection(bioInfo,dirFile)
%细菌获得每个通道荧光的平均强度后 进行串色校正
%需要mip文件 shuai Yang 2020.05.19

load([dirFile,'\mip.mat']);
disp('cell fluo intensity crosstalk correction')
% Ij, (the subscripts, j = a, ..., e, refer to GFP, YFP, OFP, RFP, and IRFP
% channels, respectively)
% Weighting Coefficient of Each FP in Different Fluorescent Channels
% FP GFP YFP OFP RFP IRFP
% sfGFP 1.67 3.35 0.18 9.60 × 10−5 0
% Venus 1.72 2.66 0.46 6.60 × 10−4 9.80 × 10−6
% CyOFP1 6.10 × 10−3 4.60 × 10−2 0.38 2.6 × 10−2 2.20 × 10−4
% mScarletI 3.70 × 10−2 6.00 × 10−3 0.26 0.73 5.80 × 10−3
% smURFP 2.70 × 10−5 1.90 × 10−5 2.10 × 10−5 1.30 × 10−5 0.17
% coeffMatrix = [1.67, 3.35, 0.18, 9.60E-5, 0;...
%     1.72, 2.66, 0.46, 6.60E-4, 9.80E-6;...
%     6.10E-3, 4.60E-2, 0.38, 2.60E-2, 2.20E-4;...
%     3.70E-2, 6.00E-3, 0.26, 0.73, 5.80E-3;...
%     2.70E-5, 1.90E-5, 2.10E-5, 1.30E-5, 0.17];

%G:\2020-06-11-FPs crosstalk correction ys，再IP33上用细菌重新测得串光系数
%每个蛋白在自己的通道的系数为1，<10-3 为0
% coeffMatrix = [...
%     1,       2.83E-2, 7.22E-2, 3.47E-3, 0;...
%     0.89,    1,       0.23,    0,       0;...
%     1.14E-2, 2.44E-2, 1,       0.27,    0;...
%     5.19E-3, 0,       5.86E-2, 1,       0;...
%     0,       0,       0,       0,       1];

% F:\2022-10-13-IP33 多通道拍摄串光校正 ys\crosstalk 
% 重新对GFP OFP RFP 3个通道拍摄的串光系数进行校正 细菌荧光测定
% 细菌荧光测定 IP 33 滤光系统
% 'Wheel-A','State',’1'
% 'Wheel-B','State',’1'
%FP GFP     YFP         OFP         RFP         IRFP
% coeffMatrix = [...
%     1,      2.83E-2,    0.146,       0,          0;...
%     0.89,   1,          0.23,       0,          0;...
%     6.50E-2,2.44E-2,    1,          0.156,      0;...
%     1E-2,   0,          1E-1,       1,          0;...
%     0,      0,          0,          0,          1];

% F:\2022-10-17-IP33 2通道拍摄串光校正 ys
% F:\2022-11-05-crosstalk -ys
% 对IP33 2个通道拍摄 指GFPCyOFP 或者GFPRFP 拍摄选择的滤光片系统下的串光校正
% 细菌荧光测定 IP 33 滤光系统
% 'Wheel-A','State',’4'
% 'Wheel-B','State',’3'
coeffMatrix = [...
    1,      2.83E-2,    0.065,      0,          0;...
    0.89,   1,          0.23,       0,          0;...
    1.55E-2,2.44E-2,    1,          0.123,      0;...
    1E-2,   0,          1.1E-1,       1,          0;...
    0,      0,          0,          0,          1];

if isempty(bioInfo)
    return
end

%fluoChannels_c,allFluoName,allCorrectedFluoName需要对应 目前5个通道校正
fluoChannels_c = {'sfGFP','Venus','CyOFP','mScarletI','TDsmURFP'};%无PVD，PVD暂不校正
fluoName_c = {'intsfGFP','intVenus','intCyOFP',...
    'intmScarletI','intTDsmURFP'};
correctedFluoName = {'intsfGFP_c','intVenus_c',...
    'intCyOFP_c','intmScarletI_c','intTDsmURFP_c'};% c represents for corrected

imShotPara = zeros(numel(fluoName_c),2);%column1激发光强，column2 曝光时间
for iCh = 1: numel(fluoChannels_c)
    fluoChannel = fluoChannels_c{iCh};
    switch fluoChannel
        case 'sfGFP'
            if ischar(mip.ImagePara.GFPexcitationIntensity)% case for input of 字符 char
                imShotPara(iCh,1) = str2double(mip.ImagePara.GFPexcitationIntensity);
            else
                imShotPara(iCh,1) = mip.ImagePara.GFPexcitationIntensity;
            end
            imShotPara(iCh,2) = mip.ImagePara.GFPexposeTime;
        case 'Venus'
            if ischar(mip.ImagePara.YFPexcitationIntensity)
                imShotPara(iCh,1) = str2double(mip.ImagePara.YFPexcitationIntensity);
            else
                imShotPara(iCh,1) = mip.ImagePara.YFPexcitationIntensity;
            end
            imShotPara(iCh,2) = mip.ImagePara.YFPexposeTime;
        case 'CyOFP'
            % if ischar(mip.ImagePara.cyOFPexcitationIntensity)
            %     imShotPara(iCh,1) = str2double(mip.ImagePara.cyOFPexcitationIntensity);
            % else
            %     imShotPara(iCh,1) = mip.ImagePara.cyOFPexcitationIntensity;
            % end
            % imShotPara(iCh,2) = mip.ImagePara.cyOFPexposeTime;

            %双相机CyOFP读取的是GFP的信息 Shuai Yang 2022/9/2
            if ischar(mip.ImagePara.GFPexcitationIntensity)% case for input of 字符 char
                imShotPara(iCh,1) = str2double(mip.ImagePara.GFPexcitationIntensity);
            else
                imShotPara(iCh,1) = mip.ImagePara.GFPexcitationIntensity;
            end
            imShotPara(iCh,2) = mip.ImagePara.GFPexposeTime;
        case 'mScarletI'
            if ischar(mip.ImagePara.RFPexcitationIntensity)
                imShotPara(iCh,1) = str2double(mip.ImagePara.RFPexcitationIntensity);
            else
                imShotPara(iCh,1) = mip.ImagePara.RFPexcitationIntensity;
            end
            imShotPara(iCh,2) = mip.ImagePara.RFPexposeTime;
        case 'TDsmURFP'
            if ischar(mip.ImagePara.iRFPexcitationIntensity)
                imShotPara(iCh,1) = str2double(mip.ImagePara.iRFPexcitationIntensity);
            else
                imShotPara(iCh,1) = mip.ImagePara.iRFPexcitationIntensity;
            end
            imShotPara(iCh,2) = mip.ImagePara.iRFPexposeTime;
    end
end


[channelIdx,bioInfo] = findCorrectionChannel(bioInfo,correctedFluoName);
imShotPara = imShotPara(channelIdx,:);
channelNum = sum(channelIdx(:));
% if channelNum <=3
%     coeffMatrix = coeffMatrix1;
% end
channelTag = find(channelIdx);
correctionMatrix = zeros (channelNum,channelNum);
for iCh = 1: channelNum
    correctionMatrix(iCh,:) = coeffMatrix(channelTag(iCh),channelTag);
end

for iInfo = 1:numel(bioInfo)
    
    if  isempty(bioInfo(iInfo).Centroid)
        continue
    else
        cellNum =  size((bioInfo(iInfo).Centroid),1);
        I = zeros(cellNum,channelNum);
    end
    
    for iCh = 1: channelNum
        fluoChannel = fluoName_c{channelTag(iCh)}(4:end);
        % I(:,iCh) = bioInfo(iInfo).(fluoName_c{channelTag(iCh)})-100;% 减背景
        I(:,iCh) = bioInfo(iInfo).(fluoName_c{channelTag(iCh)})-bioInfo(iInfo).(['BG_',fluoChannel]);
        I(:,iCh) = I(:,iCh)/(imShotPara(iCh,1)*imShotPara(iCh,2));%除以拍摄参数进行归一
    end
%     I_corrected = I*inv(correctionMatrix);
    I_corrected = I/correctionMatrix;% I*inv(correctionMatrix)
    for iCh = 1: channelNum
        bioInfo(iInfo).(correctedFluoName{channelTag(iCh)}) = I_corrected(:,iCh)*(imShotPara(iCh,1)*imShotPara(iCh,2));
        %校正后乘以拍摄参数进行还原
    end

end
% save(strcat(dirFile,'\bioInfo.mat'),'bioInfo','-v7.3');
end

function [channelIdx,bioInfo] = findCorrectionChannel(bioInfo,allCorrectedFluoName)
% correction channel 初始化 赋值为空
for iCh_c = 1: numel(allCorrectedFluoName)
    bioInfo(1).(allCorrectedFluoName{iCh_c}) = [];
end

channelIdx = false (1,5);
for iInfo = 1:numel(bioInfo)
    
    if  isempty(bioInfo(iInfo).Centroid)
        continue
    end
    if ~isempty (bioInfo(iInfo).intsfGFP)
        channelIdx(:,1) = true; %GFP channel Fluo
    end
    if ~isempty (bioInfo(iInfo).intVenus)
        channelIdx(:,2) = true; %YFP channel Fluo
    end
    if ~isempty (bioInfo(iInfo).intCyOFP)
        channelIdx(:,3) = true; %OFP channel Fluo
    end
    if ~isempty (bioInfo(iInfo).intmScarletI)
        channelIdx(:,4) = true; %RFP channel Fluo
    end
    if ~isempty (bioInfo(iInfo).intTDsmURFP)
        channelIdx(:,5) = true;%iRFP channel Fluo
    end
    return
end
end