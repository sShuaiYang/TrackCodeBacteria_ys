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
% allData里面含两类数据，一个是从mask得到的bacInfo
% 第一列为time(min),第二列为majorLength,第三列为minorLength,第四列为CyOFP,第五列为GFP,第六列为mScalet,第七列为RFP,第八列为CyPet,第九列为Venus,第十列为Ametrine一一对应，没有数据的地方用NaN表示
% 第一列为time(min),第二列为majorLength,第三列为minorLength,第四列为CyOFP,第五列为GFP,第六列为mScalet,第七列为RFP，一一对应，没有数据的地方用NaN表示，第八列为Orientation,第九列为FilledArea,
% 第十列为growthRate,第十一列为fpProduceGFP，第十二列为treeSize,第十三列为generation,第十四列为CyPet,第十五列为Venus,第十六列为Ametrine，第十七列为fpProduceCyOFP,第十八列为fpProducemScalet
% 第十九列为fpProduceRFP，第二十列为fpProduceCyPet，第二十一列为fpProduceVenus,第二十二列为fpProduceAmetrine
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
