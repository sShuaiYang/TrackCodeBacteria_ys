function parameterAll=getFeatureParameter(time)
% now paraSheet has 22 column
% 1,2-meanSpeedHigh          3,4-meanSpeedLow
% 5,6-meanVelocityHigh       7,8-meanVelocityLow
% 9,10-jumpRatio              11 -orientationVelocity
% 12,13,14-msdSlope             15,16,17 SumDistance/allFrame
% 18-growthRate              19-averageLength
% 20-standPercent            21-time
% 22-frameNum                23-varOfLenthFit
%jumpRatio就是指jump的Frame总数除以Frame总数
paraNum=23;
disp('please choose the document where allData puts')
dirData=uigetdir();
% disp('please choose the save document')
dirSave=uigetdir();
% dirData='D:\细菌数据库，new 2012-10-13\数据allData\0010--2012-04-18';
nameList=dir(dirData);
allDataInfo=[];
allDataList=[];
threShold=0.15;
if numel(time)~=size(nameList,1)-2
    disp('Wrong time number')
    return
end
smallDataNumInfo=zeros(1,size(time,2));
featureParaNum=0;
serialNum=[];
for i=1:size(nameList,1)-2;
    fileName=strcat(dirData,'\','allData',num2str(i));
    load(fileName);
    smallDataNumInfo(1,i)=size(allData.dataInfo,2);
    for iCell=1:size(allData.dataInfo,2)
        featureParaNum=featureParaNum+1;
        featurePara{featureParaNum}.time=time(i);
        featurePara{featureParaNum}.cellInfo=[i,allData.dataList{iCell}];
        featurePara{featureParaNum}.frameNum=size(allData.dataInfo{iCell}.denosieData,1);
        serialNum=[serialNum;featurePara{featureParaNum}.cellInfo];
    end
    allDataInfo=[allDataInfo,allData.dataInfo];
    for iList=1:numel(allData.dataList)
        iDataList(iList,:)=allData.dataList{iList};
    end
    iDataList(:,3)=i;
    allDataList=[allDataList;iDataList];
    iDataList=[];
end
for iAllData=2:size(smallDataNumInfo,2)
    smallDataNumInfo(1,iAllData)=smallDataNumInfo(1,iAllData)+smallDataNumInfo(1,iAllData-1);
end
for iCell=1:size(allDataInfo,2)
    cellInfo=allDataInfo{iCell};
    [featurePara{iCell}.meanSpeedHigh,featurePara{iCell}.meanSpeedLow,featurePara{iCell}.meanVelocityHigh,featurePara{iCell}.meanVelocityLow,featurePara{iCell}.jumpRatio,featurePara{iCell}.omiga,featurePara{iCell}.msdSlope,featurePara{iCell}.averageVelocity]=getCaseParameter(cellInfo,threShold);
    [featurePara{iCell}.growthRate,featurePara{iCell}.averageLength,featurePara{iCell}.inclinationAngle,featurePara{iCell}.standPercent,featurePara{iCell}.standFrameNum,featurePara{iCell}.allFrameNum,featurePara{iCell}.varOfFit]=getLengthParameter(cellInfo.denosieData(:,9));
end
standFrameInfo=zeros(1,size(smallDataNumInfo,2));allFrameInfo=zeros(1,size(smallDataNumInfo,2));
paraSheet=zeros(size(allDataInfo,2),paraNum);
for iCell=1:size(featurePara,2)
    paraSheet(iCell,:)=[featurePara{iCell}.meanSpeedHigh,featurePara{iCell}.meanSpeedLow,featurePara{iCell}.meanVelocityHigh,featurePara{iCell}.meanVelocityLow,...
        featurePara{iCell}.jumpRatio,featurePara{iCell}.omiga,featurePara{iCell}.msdSlope,featurePara{iCell}.averageVelocity,featurePara{iCell}.growthRate,...
        featurePara{iCell}.averageLength,featurePara{iCell}.standPercent,featurePara{iCell}.time,featurePara{iCell}.frameNum,featurePara{iCell}.varOfFit];
    for iNumInfo=1:size(smallDataNumInfo,2)
        if iCell<=smallDataNumInfo(iNumInfo)
            standFrameInfo(iNumInfo)=standFrameInfo(iNumInfo)+featurePara{iCell}.standFrameNum;
            allFrameInfo(iNumInfo)=allFrameInfo(iNumInfo)+featurePara{iCell}.allFrameNum;
            break;
        end
    end
end
FieldMeanStandPercent=standFrameInfo./allFrameInfo;
parameterAll.FieldMeanStandPercent=FieldMeanStandPercent;
parameterAll.featurePara=featurePara;
parameterAll.paraSheet=paraSheet;
parameterAll.serialNum=serialNum;
parameterAll.time=time;
save(strcat(dirSave,'\','parameterResult','.mat'),'parameterAll');
end
function [meanSpeedHigh,meanSpeedLow,meanVelocityHigh,meanVelocityLow,jumpRatio,omiga,msdSlope,averageVelocity]=getCaseParameter(cellInfo,threShold)
continueMaxFrames=20000;
positionData1=cellInfo.denosieData(:,2:3);
velocityDataP1=cellInfo.velocityData(:,4);
velocityDataP1(end)=[];
[meanSpeedHigh(1,1),meanSpeedLow(1,1),meanVelocityHigh(1,1),meanVelocityLow(1,1),jumpRatio(1,1)]=P1CaseParameter(positionData1,velocityDataP1,threShold);
positionData2=cellInfo.denosieData(:,4:5);
velocityDataP2=cellInfo.velocityData(:,9);
velocityDataP2(end)=[];
positionDataCentroid=cellInfo.denosieData(:,6:7);
[meanSpeedHigh(1,2),meanSpeedLow(1,2),meanVelocityHigh(1,2),meanVelocityLow(1,2),jumpRatio(1,2)]=P1CaseParameter(positionData2,velocityDataP2,threShold);
orientationData=cellInfo.denosieData(:,8);
omiga=abs(diff(orientationData));
omiga=mean(omiga);
timeDelay=logspace(1,log10(continueMaxFrames),60);
[msd1,tao]=getMSD(positionData1,timeDelay);
[msd2,~]=getMSD(positionData2,timeDelay);
[msdCentroid,~]=getMSD(positionDataCentroid,timeDelay);
%here get 1,2 and centroid MSDslope by linear fitting start from msd=0.1
tao1=tao;tao2=tao;taoCentroid=tao;
if ~isempty(msd1(msd1>0.1))
tao1(msd1<0.1)=[];msd1(msd1<0.1)=[];
end
p1=polyfit(log(tao1),log(msd1),1);
msdSlope(1,1)=p1(1);
if ~isempty(msd2(msd2>0.1))
tao2(msd2<0.1)=[];msd2(msd2<0.1)=[];
end
p2=polyfit(log(tao2),log(msd2),1);
msdSlope(1,2)=p2(1);
if ~isempty(msdCentroid(msdCentroid>0.1))
taoCentroid(msdCentroid<0.1)=[];msdCentroid(msdCentroid<0.1)=[];
end
pCentroid=polyfit(log(taoCentroid),log(msdCentroid),1);
msdSlope(1,3)=pCentroid(1);
averageVelocity(1,1)=sqrt(sum((positionData1(1,:)-positionData1(end,:)).^2))/(size(positionData1,1)-1);
averageVelocity(1,2)=sqrt(sum((positionData2(1,:)-positionData1(end,:)).^2))/(size(positionData2,1)-1);
averageVelocity(1,3)=sqrt(sum((positionDataCentroid(1,:)-positionDataCentroid(end,:)).^2))/(size(positionDataCentroid,1)-1);
end
function [meanSpeedHigh,meanSpeedLow,meanVelocityHigh,meanVelocityLow,jumpRatio]=P1CaseParameter(positionData,velocityDataP1,threShold)
STATSJump=regionprops(velocityDataP1>=threShold,velocityDataP1,'pixelIdxList');
STATSwalk=regionprops(velocityDataP1<threShold,velocityDataP1,'pixelIdxList');
jumpCase=size(STATSJump,1);
walkCase=size(STATSwalk,1);
meanSpeedHigh=[];
meanSpeedLow=[];
meanVelocityHigh=[];
meanVelocityLow=[];
frameWeightFactor=[];
for i=1:jumpCase
    meanSpeedHigh=[meanSpeedHigh;mean(abs(velocityDataP1(STATSJump(i).PixelIdxList)))];
    traveledDistance=sqrt(sum((positionData(STATSJump(i).PixelIdxList(1),:)-positionData(STATSJump(i).PixelIdxList(end)+1,:)).^2));
    meanVelocityHigh=[meanVelocityHigh;traveledDistance/numel(STATSJump(i).PixelIdxList)];
end
for i=1:walkCase
    meanSpeedLow=[meanSpeedLow;mean(abs(velocityDataP1(STATSwalk(i).PixelIdxList)))];
    traveledDistance=sqrt(sum((positionData(STATSwalk(i).PixelIdxList(1),:)-positionData(STATSwalk(i).PixelIdxList(end)+1,:)).^2));
    meanVelocityLow=[meanVelocityLow;traveledDistance/numel(STATSwalk(i).PixelIdxList)];
    frameWeightFactor=[frameWeightFactor;numel(STATSwalk(i).PixelIdxList)/numel(velocityDataP1)];
end
meanSpeedHigh=mean(meanSpeedHigh);
meanSpeedLow=mean(meanSpeedLow);
meanVelocityHigh=mean(meanVelocityHigh);
meanVelocityLow=sum(meanVelocityLow.*frameWeightFactor); %here and a frame number weight to the calculation of mean velocity of every walk case.
% meanVelocityLow=mean(meanVelocityLow);
jumpRatio=numel(velocityDataP1(velocityDataP1>=threShold))/size(velocityDataP1,1);
end
function [msd,tao]=getMSD(position,timeDelay)
positionSize=size(position,1);
timeDelay(timeDelay>positionSize-1000)=[];
tao=[1,2,3,4,5,6,7,8,9,fix(timeDelay)];
msd=[];
for iTao=1:size(tao,2)
    pos_pre=position(1:end-tao(iTao),:);
    pos_next=position(1+tao(iTao):end,:);
    dataTemp=(pos_next-pos_pre).^2;
    msd=[msd,mean(dataTemp(:,1)+dataTemp(:,2))];
end
end
function [growthRate,averageLength,inclinationAngle,standPercent,standFrameNum,allFrameNum,varOfFit]=getLengthParameter(lengthInfo)
p=polyfit(1:numel(lengthInfo),lengthInfo',1);
growthRate=p(1);
possibleLength=polyval(p,1:numel(lengthInfo))';
varOfFit=sum((possibleLength-lengthInfo).^2)/numel(lengthInfo);
if (varOfFit>0.2&&numel(lengthInfo)<6000)||varOfFit>3||p(1)<0
    if max(lengthInfo)<30
        inclinationAngle=90*ones(numel(lengthInfo),1);
        standFrameNum=numel(lengthInfo);
    else
        lengthRatio=lengthInfo./max(lengthInfo);
        inclinationAngle=acos(lengthRatio)*180/pi;
        standFrameNum=numel(inclinationAngle(inclinationAngle>20));
    end
else
    if max(lengthInfo)<30
        inclinationAngle=90*ones(numel(lengthInfo),1);
        standFrameNum=numel(lengthInfo);
    else
        lengthRatio=lengthInfo./possibleLength;
        lengthRatio(lengthRatio>1)=1;
        inclinationAngle=acos(lengthRatio)*180/pi;
        standFrameNum=numel(inclinationAngle(inclinationAngle>20));
    end
end
allFrameNum=numel(lengthInfo); standPercent=standFrameNum./allFrameNum;
averageLength=mean(lengthInfo);
end