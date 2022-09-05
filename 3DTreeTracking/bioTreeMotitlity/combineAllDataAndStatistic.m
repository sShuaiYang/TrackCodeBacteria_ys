function combineAllDataAndStatistic(clusterInf)
dirSave=uigetdir();
cd(dirSave)
save(strcat(dirSave,'\clusterInf','.mat'),'clusterInf')
for i=1:numel(clusterInf)
    if ~isempty(clusterInf{i})
        NewAllData=makeNewAllData(clusterInf{i});
        mkdir(num2str(i));
        saveFile=strcat(dirSave,'\',num2str(i));
        %here caculate the histResult
        velocityList=allDataforHistVelocity(NewAllData);
        histResult=histVelocityAll(velocityList);
        %here caculate the timeResult
%         AllTimeResult=caculateAllDataforTimeSeries(NewAllData);
%         plotResultsAndSave(strcat(saveFile,'\Graphic'),saveFile,AllTimeResult,histResult,[0,0],i);
        result=JumpTimeHist(NewAllData,saveFile,0.15);
        save(strcat(saveFile,'velocityTimeStatistic'),'result')
        close all
    end
end
end

function NewAllData=makeNewAllData(oneCluster)
t=0;
for i=1:size(oneCluster,1)
    thisData=oneCluster(i,:);
    if numel(thisData(thisData>0))>0
        load(strcat('allData',num2str(i)))
        for j=1:numel(thisData(thisData>0))
            t=t+1;
            giveNum=oneCluster(i,j);
            giveNum2=[fix(giveNum) round((giveNum-fix(giveNum))*10)];
            for p=1:numel(allData.dataList)
                 equalorNot=allData.dataList{p}-giveNum2;
                 if numel(equalorNot(equalorNot==0))==2
                    NewAllData(t)=allData.dataInfo(p);
                    break
                end
            end
        end
        clear allData
    end
end
end
function velocityList=allDataforHistVelocity(NewAllData)
velocityList=[];
for i=1:numel(NewAllData)
    velocityList=[velocityList;NewAllData{i}.velocityData];
end
velocityList(:,1)=1:size(velocityList,1);
end

%% here caculate the timeSeries for all data
function AllTimeResult=caculateAllDataforTimeSeries(NewAllData)
continueMaxFrames=15000;                                   %here we have the Max Frame 
for i=1:numel(NewAllData)
    AllTimeResultOri{i}=timeSeriesForAll(NewAllData{i});
end
AllTimeResult.Msd=AllTimeResultOri{1}.Msd;
AllTimeResult.Correlation=AllTimeResultOri{1}.Correlation;
AllTimeResult.PowerSpectrum=AllTimeResultOri{1}.PowerSpectrum;
%here we get the hole Msd and MsdDiff
for i=2:numel(NewAllData)
    AllTimeResult.Msd=matrixPlus(AllTimeResult.Msd,AllTimeResultOri{i}.Msd);
end
AllTimeResult.Msd=AllTimeResult.Msd/numel(NewAllData);
AllTimeResult.MsdDiff=makeDiff(AllTimeResult.Msd(:,1)',AllTimeResult.Msd(:,2)',AllTimeResult.Msd(:,3)',AllTimeResult.Msd(:,4)');
    timeDelay=logspace(1,log10(continueMaxFrames),60);
    tao=[1,2,3,4,5,6,7,8,9,fix(timeDelay)]';
    tao(size(AllTimeResult.Msd,1)+1:end)=[];
    AllTimeResult.Msd(:,1)=tao;
    tao=[1,2,3,4,5,6,7,8,9,fix(timeDelay)]';
    tao(size(AllTimeResult.MsdDiff,1)+1:end)=[];
    AllTimeResult.MsdDiff(:,1)=tao;
%here we get the hole correlation and Power Spectrum
for i=2:numel(NewAllData)
    AllTimeResult.Correlation=matrixPlus(AllTimeResult.Correlation,AllTimeResultOri{i}.Correlation);
    AllTimeResult.PowerSpectrum=matrixPlus(AllTimeResult.PowerSpectrum,AllTimeResultOri{i}.PowerSpectrum);
end
AllTimeResult.Correlation(:,1)=0:size(AllTimeResult.Correlation,1)-1;
AllTimeResult.PowerSpectrum(:,1)=(0:(5000-1))/5000;
AllTimeResult.PowerSpectrum(:,2)=AllTimeResult.PowerSpectrum(:,2)/sum(AllTimeResult.PowerSpectrum(:,2));
AllTimeResult.PowerSpectrum(:,3)=AllTimeResult.PowerSpectrum(:,3)/sum(AllTimeResult.PowerSpectrum(:,3));
end
function MatrixA=matrixPlus(MatrixB,MatrixC)
a=size(MatrixB,1);
b=size(MatrixC,1);
if a>b
    MatrixB(b+1:a,:)=[];
else if a<b
        MatrixC(a+1:b,:)=[];
    end
end
MatrixA=MatrixB+MatrixC;
end

%% here we make timeSeries for one dateAll,the time is devided by 5000 and the continue Max frames is 15000
function timeResult=timeSeriesForAll(dataAll)
continueMaxFrames=15000;
timeDelay=logspace(1,log10(continueMaxFrames),60); %time series lengh need > 2000
originalData=dataAll.originalData;
denosieData=dataAll.denosieData;
velocityData=dataAll.velocityData;
[msdP1,tao1]=getMSD(denosieData(:,2:3),timeDelay);
[msdP2,~]=getMSD(denosieData(:,4:5),timeDelay);
[msdCentroid,~]=getMSD(denosieData(:,6:7),timeDelay);
timeResult.Msd=[tao1',msdP1',msdP2',msdCentroid'];
[~,tao1]=xcorr(velocityData(:,2));
corrVP1=xcorr(velocityData(:,2),velocityData(:,2))+xcorr(velocityData(:,3),velocityData(:,3));
corrVp2=xcorr(velocityData(:,7),velocityData(:,7))+xcorr(velocityData(:,8),velocityData(:,8));
corrVp1p2=xcorr(velocityData(:,2),velocityData(:,7))+xcorr(velocityData(:,3),velocityData(:,8));
n=numel(velocityData(:,2));
timeResult.Correlation=[tao1(end-(n-1):end)',corrVP1(end-(n-1):end),corrVp2(end-(n-1):end),corrVp1p2(end-(n-1):end)];

powerV1=caculatePower(velocityData(:,2),velocityData(:,3));
powerV2=caculatePower(velocityData(:,7),velocityData(:,8));
% tao=(0:(numel(powerV1)-1))/numel(powerV1);
tao=(0:(5000-1))/5000;
timeResult.PowerSpectrum=[tao',powerV1,powerV2];
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
function MsdDiff=makeDiff(tao1,msdP1,msdP2,msdCentroid)
msdP1Diff=diff(log(msdP1'));
msdP2Diff=diff(log(msdP2'));
msdCentroidDiff=diff(log(msdCentroid'));
taoDiff=diff(log(tao1'));
tao1Diff=tao1';
tao1Diff(end)=[];
MsdDiff=[tao1Diff,msdP1Diff./taoDiff,msdP2Diff./taoDiff,msdCentroidDiff./taoDiff];
end
function powerSpectrum=caculatePower(velocityX,velocityY)
powerX=abs(fft(velocityX,5000)).^2;
powerY=abs(fft(velocityY,5000)).^2;
powerSpectrum=powerX+powerY;
end

%% here is histVelocityAll   ( make more points)
function histResult=histVelocityAll(velocityList) % displacement weighted hisotgram
binSpace=logspace(-3,2,121); % the las number must be an odd, and the data nmber is (n-1)/2
angleSpace=linspace(-180,180,61); %same as above
[countP1,n1]=get1Dhist(velocityList(:,4),binSpace);
[countP2,n2]=get1Dhist(velocityList(:,9),binSpace);
p1threshold=getThreshold(n1,countP1);
p2threshold=getThreshold(n2,countP2);
[countMapP2sum,nMap]=get2Dhist(velocityList(:,4),velocityList(:,9),binSpace);
[countMapP1sum,~]=get2Dhist(velocityList(:,9),velocityList(:,4),binSpace);
histResult.p1=[n1',countP1'];
histResult.p2=[n2',countP2'];
histResult.p1p2MapP1sum=countMapP1sum;
histResult.p1p2MapP2sum=countMapP2sum';
histResult.p1p2Bin=nMap;
[p1VelocityOrientation,p2VelocityOrientation]=histOrientation(velocityList,angleSpace,p1threshold,p2threshold);
histResult.p1VelocityOrientation=p1VelocityOrientation;
histResult.p2VelocityOrientation=p2VelocityOrientation;
end
function threshold=getThreshold(n,countP)
threshold=[];
[peaks,location]=findpeaks(countP);
if size(peaks,2)==1
    return;
elseif size(peaks,2)>1
    for i=2:size(peaks,2)
        [~,index]=min(countP(location(i-1): location(i)));
        threshold=[threshold,n(location(i-1)+index-1)];
    end
    return;
end
end
function [velocityLow,velocityHigh,velocityOrientationHigh,velocityOrientationLow]=thresholdVelocity(velocity,velocityOrientaion,threshold)
velocityLow=velocity(velocity<=threshold);
velocityHigh=velocity(velocity>threshold);
velocityOrientationLow=velocityOrientaion(velocity<=threshold);
velocityOrientationHigh=velocityOrientaion(velocity>threshold);
end
function [count n]=get1Dhist(velocityData,binSpace)
nBin=binSpace;
count=[];
n=[];
for i=2:2:size(nBin,2)
    count=[count,sum(velocityData(velocityData>=nBin(i-1)&velocityData<nBin(i+1)))];
    n=[n,nBin(i)];
end
count=count./size(velocityData,1);
end
function [countMap,n]=get2Dhist(velocityP1,velocityP2,binSpace)
nBin=binSpace;
n=[];
countMap=zeros((size(nBin,2)-1)/2,(size(nBin,2)-1)/2);
for i=2:2:size(nBin,2)
    dataCol=velocityP2(velocityP1>=nBin(i-1)&velocityP1<nBin(i+1));
    for j=2:2:size(nBin,2)
        countMap(i/2,j/2)=sum(dataCol(dataCol>=nBin(j-1)&dataCol<nBin(j+1)));
    end
     n=[n,nBin(i)];
end
end
function [countP1angle,angle]=getAnglehist(velocityP1,velocityAngle,angleSpace)
countP1angle=[];
angle=[];
for i=2:2:size(angleSpace,2)
    countP1angle=[countP1angle,sum(velocityP1(velocityAngle>=angleSpace(i-1)&velocityAngle<angleSpace(i+1)))];
    angle=[angle,angleSpace(i)];
end
angleRed=angle.*(pi/180);
angle=[angleRed',angle'];
end
function [p1VelocityOrientation,p2VelocityOrientation]=histOrientation(velocityList,angleSpace,p1threshold,p2threshold)
if ~isempty(p1threshold)
    for iThreshold=1:size(p1threshold,2)
        [p1VelocityLow,p1VelocityHigh,p1VelocityOrientationHigh,p1VelocityOrientationLow]=thresholdVelocity(velocityList(:,4),velocityList(:,5),p1threshold(iThreshold));
        [countP1angle,angle]=getAnglehist(velocityList(:,4),velocityList(:,5),angleSpace);
        [countP1angleLow,~]=getAnglehist(p1VelocityLow,p1VelocityOrientationLow,angleSpace);
        [countP1angleHigh,~]=getAnglehist(p1VelocityHigh,p1VelocityOrientationHigh,angleSpace);
        p1VelocityOrientation{iThreshold}.threshold=p1threshold(iThreshold);
        p1VelocityOrientation{iThreshold}.histOrientation=[angle,countP1angle',countP1angleLow',countP1angleHigh'];
    end
elseif isempty(p1threshold)
    [countP1angle,angle]=getAnglehist(velocityList(:,4),velocityList(:,5),angleSpace);
    p1VelocityOrientation{1}.threshold=[];
    p1VelocityOrientation{1}.histOrientation=[angle,countP1angle'];
end
if ~isempty(p2threshold)
    for iThreshold=1:size(p2threshold,2)
        [p2VelocityLow,p2VelocityHigh,p2VelocityOrientationHigh,p2VelocityOrientationLow]=thresholdVelocity(velocityList(:,9),velocityList(:,10),p2threshold(iThreshold));
        [countP2angle,angle]=getAnglehist(velocityList(:,9),velocityList(:,10),angleSpace);
        [countP2angleLow,~]=getAnglehist(p2VelocityLow,p2VelocityOrientationLow,angleSpace);
        [countP2angleHigh,~]=getAnglehist(p2VelocityHigh,p2VelocityOrientationHigh,angleSpace);
        p2VelocityOrientation{iThreshold}.threshold=p2threshold(iThreshold);
        p2VelocityOrientation{iThreshold}.histOrientation=[angle,countP2angle',countP2angleLow',countP2angleHigh'];
    end
elseif isempty(p2threshold)
    [countP2angle,angle]=getAnglehist(velocityList(:,9),velocityList(:,10),angleSpace);
    p2VelocityOrientation{1}.threshold=[];
    p2VelocityOrientation{1}.histOrientation=[angle,countP2angle'];
end
end

    
