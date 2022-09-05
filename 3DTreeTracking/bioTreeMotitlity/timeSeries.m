function timeResult=timeSeries(dataAll)
continueMaxFrames=15000;
timeDelay=logspace(1,log10(continueMaxFrames),60); %time series lengh need > 2000
originalData=dataAll.originalData;
denosieData=dataAll.denosieData;
velocityData=dataAll.velocityData;
[msdP1,tao1]=getMSD(denosieData(:,2:3),timeDelay);
[msdP2,~]=getMSD(denosieData(:,4:5),timeDelay);
[msdCentroid,~]=getMSD(denosieData(:,6:7),timeDelay);
timeResult.Msd=[tao1',msdP1',msdP2',msdCentroid'];
timeResult.MsdDiff=makeDiff(tao1,msdP1,msdP2,msdCentroid);

timeResult.beforeWdencmpLen=dataAll.beforeWdencmpLen;
timeResult.afterWdencmpLen=dataAll.afterWdencmpLen;
timeResult.diffLen=abs(diff(timeResult.afterWdencmpLen));

OritationInfo=dataAll.denosieData(:,8);
timeResult.OritationInfo.Distribution=OritationInfo-mean(OritationInfo);
angleVelocity=diff(OritationInfo);
angleVelocity(angleVelocity>180)=angleVelocity(angleVelocity>180)-360;
angleVelocity(angleVelocity<-180)=angleVelocity(angleVelocity<-180)+360;
timeResult.OritationInfo.angleVelocity=angleVelocity;
[corrAngle,tao2]=xcorr(angleVelocity);
timeResult.OritationInfo.correlation=[tao2(end-(numel(angleVelocity)-1):end)',corrAngle(end-(numel(angleVelocity)-1):end)];
powerAngle=caculatePower(angleVelocity,0);
tao=(0:(numel(powerAngle)-1))/numel(powerAngle);
timeResult.OritationInfo.power=[tao',powerAngle];

% [corrVP1,tao1]=getCorrelation(velocityData(:,2:3),velocityData(:,2:3),timeDelay);
% [corrVp2,~]=getCorrelation(velocityData(:,7:8),velocityData(:,7:8),timeDelay);
% [corrVp1p2,~]=getCorrelation(velocityData(:,2:3),velocityData(:,7:8),timeDelay);
[~,tao1]=xcorr(velocityData(:,2));
corrVP1=xcorr(velocityData(:,2),velocityData(:,2))+xcorr(velocityData(:,3),velocityData(:,3));
corrVp2=xcorr(velocityData(:,7),velocityData(:,7))+xcorr(velocityData(:,8),velocityData(:,8));
corrVp1p2=xcorr(velocityData(:,2),velocityData(:,7))+xcorr(velocityData(:,3),velocityData(:,8));
n=numel(velocityData(:,2));
timeResult.Correlation=[tao1(end-(n-1):end)',corrVP1(end-(n-1):end),corrVp2(end-(n-1):end),corrVp1p2(end-(n-1):end)];

powerV1=caculatePower(velocityData(:,2),velocityData(:,3));
powerV2=caculatePower(velocityData(:,7),velocityData(:,8));
tao=(0:(numel(powerV1)-1))/numel(powerV1);
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
powerX=abs(fft(velocityX)).^2;
powerY=abs(fft(velocityY)).^2;
powerSpectrum=powerX+powerY;
end
% function [correlation,tao]=getCorrelation(velocity1,velocity2,timeDelay)
% tao=1:1000;
% correlation=[];
% % coeff=sum(mean(velocity1).*mean(velocity2));
% coeff=1;
% for iTao=1:size(tao,2)
%     velocity1_pre=velocity1(1:end-tao(iTao),:);
%     velocity2_next=velocity2(1+tao(iTao):end,:);
%      dataTemp=velocity1_pre.*velocity2_next;
%      correlation=[ correlation,mean(dataTemp(:,1)+dataTemp(:,2))];
% end
% correlation=correlation./coeff;
% end