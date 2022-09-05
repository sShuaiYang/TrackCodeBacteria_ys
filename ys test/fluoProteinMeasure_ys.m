% fluo protein measurement code 
% Shuai Yang 2022/01/18
% 'sfGFP','CyOFP','mCherry','mScarletI';
% sutter optional:
% obj.mmc_SutterLB.setProperty('Wheel-A','State','4');% state 4 GFP
% obj.mmc_SutterLB.setProperty('Wheel-B','State','3');% state 3 RFP


FP_name = 'mCherry';%'CyOFP','mCherry','mScarletI'; 
h = '5000nm';% SiO2 particel height
c = '500nM';% FP concentration nM

% h = 'xxx';% SiO2 particel height
% c = 'yakeliR';% FP concentration nM

Ex_Intensity = 60 ;% excitatio intensiy 10% default
mip.ImagePara.GFPexcitationIntensity = Ex_Intensity;
mip.ImagePara.RFPexcitationIntensity = Ex_Intensity;
mip.ImagePara.cyOFPexcitationIntensity = Ex_Intensity;

Ex_Time =  100:-10:10;% 10-100mS
dirFile = 'F:\2022-01-18-FP measure-ys';
for i = 1:numel(Ex_Time)

    mip.ImagePara.GFPexposeTime = Ex_Time(i);
    mip.ImagePara.RFPexposeTime = Ex_Time(i);
    mip.ImagePara.cyOFPexposeTime = Ex_Time(i);

    dirSave = [dirFile,'\',FP_name,'\_',c,'\_',h,'_ET',num2str(Ex_Time(i),'%03d'),'_EI',num2str(Ex_Intensity,'%04d')];
    if ~isfolder(dirSave)
        mkdir(dirSave)
    end
    mip.myDir = dirSave;
    switch FP_name
        case 'sfGFP'
            protocalTestGFPRFPYS(mip,MesStruct)
        case 'CyOFP'
            protocalTestGFPCyOFPYS(mip,MesStruct)
        case 'mCherry'
            protocalTestGFPRFPYS(mip,MesStruct)
        case 'mScarletI'
            protocalTestGFPRFPYS(mip,MesStruct)
    end
    pause(2)
end