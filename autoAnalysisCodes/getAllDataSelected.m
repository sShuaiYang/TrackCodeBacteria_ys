function [allData]=getAllDataSelected()
clc
disp('please choose the fold where you want to get the AllData')
dirAllDataFile=uigetdir();
load([dirAllDataFile,'\allData.mat']);
disp('please choose the storage folder');
% allData���溬�������ݣ�һ���Ǵ�mask�õ���bacInfo
% ��һ��Ϊtime(min),�ڶ���ΪmajorLength,������ΪminorLength,������ΪCyOFP,������ΪGFP,������ΪmScalet,������ΪRFP,�ڰ���ΪCyPet,�ھ���ΪVenus,��ʮ��ΪAmetrineһһ��Ӧ��û�����ݵĵط���NaN��ʾ
% ��һ��Ϊtime(min),�ڶ���ΪmajorLength,������ΪminorLength,������ΪCyOFP,������ΪGFP,������ΪmScalet,������ΪRFP��һһ��Ӧ��û�����ݵĵط���NaN��ʾ���ڰ���ΪOrientation,�ھ���ΪFilledArea,
% ��ʮ��ΪgrowthRate,��ʮһ��ΪfpProduceGFP����ʮ����ΪtreeSize,��ʮ����Ϊgeneration,��ʮ����ΪCyPet,��ʮ����ΪVenus,��ʮ����ΪAmetrine����ʮ����ΪfpProduceCyOFP,��ʮ����ΪfpProducemScalet
% ��ʮ����ΪfpProduceRFP���ڶ�ʮ��ΪfpProduceCyPet���ڶ�ʮһ��ΪfpProduceVenus,�ڶ�ʮ����ΪfpProduceAmetrine
dirStorageFile=uigetdir();
variableName=checkData(allData);
i=0;
variableNew=input('please type the variable____:');
variableNew=variableName{variableNew};
while isempty(variableNew) || ~strcmp(variableNew,'over')
    i=i+1;
    variableAll{2*i-1}=variableNew;
    variableAll{2*i}=input('please input the range____:');
    variableNew=input('please type the variable____:');
end
for i=1:numel(variableAll)/2
    switch variableAll{i*2-1}
        case 'tag'
            maskLine=strcmp(allData.maskResultTag,variableAll{i*2});
            trackingLine=strcmp(allData.trackingResultTag,variableAll{i*2});
        case 'tagValue'
            maskLine=allData.maskResultTagValue>=variableAll{i*2}(1) & allData.maskResultTagValue<=variableAll{i*2}(2);
            trackingLine=allData.trackingResultTagValue>=variableAll{i*2}(1) & allData.trackingResultTagValue<=variableAll{i*2}(2);
        otherwise
            iLine=getMaskResultLine(allData.maskResult,variableAll{i*2-1});
            if ~isempty(iLine)
                maskLine=iLine>=variableAll{i*2}(1) & iLine<=variableAll{i*2}(2);
            else
                maskLine=[];
            end
            iLine=getTrackingResultLine(allData.trackingResult,variableAll{i*2-1});
            trackingLine=iLine>=variableAll{i*2}(1) & iLine<=variableAll{i*2}(2);
            if strcmp(variableAll{i*2-1},'treeSize')
                treeLine=allData.treeSize>=variableAll{i*2}(1) & allData.treeSize<=variableAll{i*2}(2);
                allData.treeAll=allData.treeAll(treeLine);
                allData.treeSize=allData.treeSize(treeLine);
                try
                    generationLine=allData.correlationDataCYOFP(:,3)>=variableAll{i*2}(1) & allData.correlationDataCYOFP(:,3)<=variableAll{i*2}(2);
                    allData.correlationDataCYOFP=allData.correlationDataCYOFP(generationLine,:);
                end
                try
                    generationLine=allData.correlationDataGFP(:,3)>=variableAll{i*2}(1) & allData.correlationDataGFP(:,3)<=variableAll{i*2}(2);
                    allData.correlationDataGFP=allData.correlationDataGFP(generationLine,:);
                end
                try
                    generationLine=allData.correlationDatamScalet(:,3)>=variableAll{i*2}(1) & allData.correlationDatamScalet(:,3)<=variableAll{i*2}(2);
                    allData.correlationDatamScalet=allData.correlationDatamScalet(generationLine,:);
                end
                try
                    generationLine=allData.correlationDataRFP(:,3)>=variableAll{i*2}(1) & allData.correlationDataRFP(:,3)<=variableAll{i*2}(2);
                    allData.correlationDataRFP=allData.correlationDataRFP(generationLine,:);
                end
                try
                    generationLine=allData.correlationDataCyPet(:,3)>=variableAll{i*2}(1) & allData.correlationDataCyPet(:,3)<=variableAll{i*2}(2);
                    allData.correlationDataCyPet=allData.correlationDataCyPet(generationLine,:);
                end
                try
                    generationLine=allData.correlationDataVenus(:,3)>=variableAll{i*2}(1) & allData.correlationDataVenus(:,3)<=variableAll{i*2}(2);
                    allData.correlationDataVenus=allData.correlationDataVenus(generationLine,:);
                end
                try
                    generationLine=allData.correlationDatamAmetrine(:,3)>=variableAll{i*2}(1) & allData.correlationDatamAmetrine(:,3)<=variableAll{i*2}(2);
                    allData.correlationDatamAmetrine=allData.correlationDatamAmetrine(generationLine,:);
                end
            end
    end
    if ~isempty(maskLine)
        allData.maskResult=allData.maskResult(maskLine,:);
        allData.maskResultTag=allData.maskResultTag(maskLine,:);
        allData.maskResultTagValue=allData.maskResultTagValue(maskLine,:);
    end
    allData.trackingResult=allData.trackingResult(trackingLine,:);
    allData.trackingResultTag=allData.trackingResultTag(trackingLine,:);
    allData.trackingResultTagValue=allData.trackingResultTagValue(trackingLine,:);
    
    
end
mkdir(dirStorageFile);
cd(dirStorageFile)
save([dirStorageFile,'\allData.mat'],'allData')
end
function iLine=getMaskResultLine(dataAll,inputVari)
switch inputVari
    case 'time'
        iLine=dataAll(:,1);
    case 'maxLen'
        iLine=dataAll(:,2);
    case 'minLen'
        iLine=dataAll(:,3);
    case 'CyOFP'
        iLine=dataAll(:,4);
    case 'GFP'
        iLine=dataAll(:,5);
    case 'mScalet'
        iLine=dataAll(:,6);
    case 'RFP'
        iLine=dataAll(:,7);
    case 'CyPet'
        iLine=dataAll(:,8);
    case 'Venus'
        iLine=dataAll(:,9);
    case 'mAmetrine'
        iLine=dataAll(:,10);
    case 'redControl'
        iLine=dataAll(:,11);
    case 'blueControl'
        iLine=dataAll(:,12);
    case 'greenControl'
        iLine=dataAll(:,13);
    otherwise
        iLine=[];
end
end
function iLine=getTrackingResultLine(dataAll,inputVari)
if isempty(dataAll)
    iLine=[];
    return
end
switch inputVari
    case 'time'
        iLine=dataAll(:,1);
    case 'maxLen'
        iLine=dataAll(:,2);
    case 'minLen'
        iLine=dataAll(:,3);
    case 'CyOFP'
        iLine=dataAll(:,4);
    case 'GFP'
        iLine=dataAll(:,5);
    case 'mScalet'
        iLine=dataAll(:,6);
    case 'RFP'
        iLine=dataAll(:,7);
    case 'angle'
        iLine=dataAll(:,8);
    case 'area'
        iLine=dataAll(:,9);
    case 'growthRate'
        iLine=dataAll(:,10);
    case 'fpProduceGFP'
        iLine=dataAll(:,11);
    case 'treeSize'
        iLine=dataAll(:,12);
    case 'generation'
        iLine=dataAll(:,13);
    case 'CyPet'
        iLine=dataAll(:,14);
    case 'Venus'
        iLine=dataAll(:,15);
    case 'mAmetrine'
        iLine=dataAll(:,16);
    case 'fpProduceCyOFP'
        iLine=dataAll(:,17);
    case 'fpProducemScalet'
        iLine=dataAll(:,18);
    case 'fpProduceRFP'
        iLine=dataAll(:,19);
    case 'fpProduceCyPet'
        iLine=dataAll(:,20);
    case 'fpProduceVenus'
        iLine=dataAll(:,21);
    case 'fpProducemAmetrine'
        iLine=dataAll(:,22);
    case 'redControl'
        iLine=dataAll(:,23);
    case 'blueControl'
        iLine=dataAll(:,24);
    case 'greenControl'
        iLine=dataAll(:,25);
end
end
function variableName=checkData(allData)
% allData���溬�������ݣ�һ���Ǵ�mask�õ���bacInfo
% ��һ��Ϊtime(min),�ڶ���ΪmajorLength,������ΪminorLength,������ΪCyOFP,������ΪGFP,������ΪmScalet,������ΪRFP,�ڰ���ΪCyPet,�ھ���ΪVenus,��ʮ��ΪAmetrineһһ��Ӧ��û�����ݵĵط���NaN��ʾ
% ��һ��Ϊtime(min),�ڶ���ΪmajorLength,������ΪminorLength,������ΪCyOFP,������ΪGFP,������ΪmScalet,������ΪRFP��һһ��Ӧ��û�����ݵĵط���NaN��ʾ���ڰ���ΪOrientation,�ھ���ΪFilledArea,
% ��ʮ��ΪgrowthRate,��ʮһ��ΪfpProduceGFP����ʮ����ΪtreeSize,��ʮ����Ϊgeneration,��ʮ����ΪCyPet,��ʮ����ΪVenus,��ʮ����ΪAmetrine����ʮ����ΪfpProduceCyOFP,��ʮ����ΪfpProducemScalet
% ��ʮ����ΪfpProduceRFP���ڶ�ʮ��ΪfpProduceCyPet���ڶ�ʮһ��ΪfpProduceVenus,�ڶ�ʮ����ΪfpProduceAmetrine,�ڶ�ʮ����ΪredControl,�ڶ�ʮ����ΪblueControl,�ڶ�ʮ����ΪgreenControl
variableName{1}='time';
variableName{2}='maxLen';
variableName{3}='minLen';
variableName{4}='CyOFP';
variableName{5}='GFP';
variableName{6}='mScalet';
variableName{7}='RFP';
variableName{8}='CyPet';
variableName{9}='Venus';
variableName{10}='mAmetrine';
variableName{11}='angle';
variableName{12}='area';
variableName{13}='growthRate';
variableName{14}='tag';
variableName{15}='tagValue';
variableName{16}='treeSize';
variableName{17}='generation';
variableName{18}='fpProduceCyOFP';
variableName{19}='fpProduceGFP';
variableName{20}='fpProducemScalet';
variableName{21}='fpProduceRFP';
variableName{22}='fpProduceCyPet';
variableName{23}='fpProduceVenus';
variableName{24}='fpProducemAmetrine';
variableName{25}='redControl';
variableName{26}='blueControl';
variableName{27}='greenControl';
maskImage=allData.maskResult;
trackingResult=allData.trackingResult;
chooseIndex=true(27,1);
for i=4:10
    if all(isnan(maskImage(:,i))) 
        chooseIndex(i)=0;
        chooseIndex(i+14)=0;
    end
end
if isempty(trackingResult)
    chooseIndex=chooseIndex([1:10,25:27]);
    variableName=variableName([1:10,25:27]);
else
    for i=25:27
        if all(isnan(trackingResult(:,i-2)))
            chooseIndex(i)=0;
        end
    end
end
variableName=variableName(chooseIndex);
variableNum=numel(variableName);
line1=[];
line2=[];
line3=[];
for i=1:variableNum
    if i<=9
        line1=[line1,num2str(i),'. ',variableName{i},'  ;'];
    else
        if i<=18
            line2=[line2,num2str(i),'. ',variableName{i},'  ;'];
        else
            if i<=27
                line3=[line3,num2str(i),'. ',variableName{i},'  ;'];
            end
        end
    end
end
disp('now please input the variables and their range, if all have done, type("over")');
disp(line1)
disp(line2)
if ~isempty(line3)
    disp(line3)
end
end
