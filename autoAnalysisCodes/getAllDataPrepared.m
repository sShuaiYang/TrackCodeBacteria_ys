function [allData] = getAllDataPrepared(dirStorageFile)
%GETALLDATAPREPARED Summary of this function goes here
%   Detailed explanation goes here
clc
disp('please input the experiment amount you want to add')
experimentNum=input('____:');
for i=1:experimentNum
    disp(['experiment',num2str(i,'%02.f')])
    dirFile{i}=uigetdir();
    nameList=dir(dirFile{i});
%     field{i}=1:numel(nameList)-3;
    field{i}=1:225;
end
allData.maskResult=[];
allData.maskResultTag=[];
allData.maskResultTagValue=[];
allData.trackingResult=[];
allData.trackingResultTag=[];
allData.trackingResultTagValue=[];
allData.correlationDataCYOFP=[];
allData.correlationDataGFP=[];
allData.correlationDatamScalet=[];
allData.correlationDataRFP=[];
allData.treeAll=[];
allData.treeSize=[];
allData.correlationDataCyPet=[];
allData.correlationDataVenus=[];
allData.correlationDatamAmetrine=[];
% allData���溬�������ݣ�һ���Ǵ�mask�õ���bacInfo
% ��һ��Ϊtime(min),�ڶ���ΪmajorLength,������ΪminorLength,������ΪCyOFP,������ΪGFP,������ΪmScalet,������ΪRFP,�ڰ���ΪCyPet,�ھ���ΪVenus,��ʮ��ΪAmetrineһһ��Ӧ��û�����ݵĵط���NaN��ʾ
% ��һ��Ϊtime(min),�ڶ���ΪmajorLength,������ΪminorLength,������ΪCyOFP,������ΪGFP,������ΪmScalet,������ΪRFP��һһ��Ӧ��û�����ݵĵط���NaN��ʾ���ڰ���ΪOrientation,�ھ���ΪFilledArea,
% ��ʮ��ΪgrowthRate,��ʮһ��ΪfpProduceGFP����ʮ����ΪtreeSize,��ʮ����Ϊgeneration,��ʮ����ΪCyPet,��ʮ����ΪVenus,��ʮ����ΪAmetrine����ʮ����ΪfpProduceCyOFP,��ʮ����ΪfpProducemScalet
% ��ʮ����ΪfpProduceRFP���ڶ�ʮ��ΪfpProduceCyPet���ڶ�ʮһ��ΪfpProduceVenus,�ڶ�ʮ����ΪfpProduceAmetrine
treeAll=[];
for i=1:experimentNum
    nameList=dir(dirFile{i});
    for iField=field{i}
        try
            allData1=load([dirFile{i},'\',nameList(iField+2).name,'\allData']);
            allData.maskResult=[allData.maskResult;allData1.allData.maskResult];
            allData.maskResultTagValue=[allData.maskResultTagValue;allData1.allData.maskResultTagValue];
            allData.maskResultTag=[allData.maskResultTag;allData1.allData.maskResultTag];
            allData.trackingResult=[allData.trackingResult;allData1.allData.trackingResult];
            allData.trackingResultTagValue=[allData.trackingResultTagValue;allData1.allData.trackingResultTagValue];
            allData.trackingResultTag=[allData.trackingResultTag;allData1.allData.trackingResultTag];
            allData.correlationDataCYOFP=[allData.correlationDataCYOFP;allData1.allData.correlationDataCYOFP];
            allData.correlationDataGFP=[allData.correlationDataGFP;allData1.allData.correlationDataGFP];
            allData.correlationDatamScalet=[allData.correlationDatamScalet;allData1.allData.correlationDatamScalet];
            allData.correlationDataRFP=[allData.correlationDataRFP;allData1.allData.correlationDataRFP];
            allData.correlationDataCyPet=[allData.correlationDataCyPet;allData1.allData.correlationDataCyPet];
            allData.correlationDataVenus=[allData.correlationDataVenus;allData1.allData.correlationDataVenus];
            allData.correlationDatamAmetrine=[allData.correlationDatamAmetrine;allData1.allData.correlationDatamAmetrine];
            allData.treeAll=[allData.treeAll,allData1.allData.treeAll];
            allData.treeSize=[allData.treeSize;allData1.allData.treeSize];
        end
    end
end
mkdir([dirStorageFile,'\allData']);
cd([dirStorageFile,'\allData'])
save([dirStorageFile,'\allData\allData.mat'],'allData')
end
