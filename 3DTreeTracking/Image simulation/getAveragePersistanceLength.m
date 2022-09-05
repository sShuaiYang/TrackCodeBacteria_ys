% function getAllDataPlot()
% dirFile=uigetdir();
% nameList=dir(dirFile);
% for i=1:numel(nameList)-2
%     if i==2
%         continue
%     end
%     cd(strcat(dirFile,'\',nameList(i+2).name))
%     dirFrameInfo=strcat(dirFile,'\',nameList(i+2).name,'\bestBioTree\bioTree_5');
%     load(dirFrameInfo)
%     cosSitaAndLength=getAveragePersistanceLength(bioTree,'simulation');
%     save('cosSitaAndLength','cosSitaAndLength')
% end
% end
function cosSitaAndLength=getAveragePersistanceLength(bioTree,bioTreeType)
cosSitaAndLength=[];
% for iframe=1:size(bioTree,2)
for iframe=1:7000
    for iRoot=1:size(bioTree{iframe}.root,2)
        if strcmp(bioTreeType,'simulation')
            centroidInfo=bioTree{iframe}.root{iRoot}.traceInfo.traceCentroid;
        end
        if strcmp(bioTreeType,'real')
            measurment=bioTree{iframe}.root{iRoot}.traceInfo.measurment;
            centroidInfo=[];
            for iCen=1:numel(measurment)
                centroidInfo=[centroidInfo;measurment{iCen}(1).Centroid];
            end
        end
        cosSitaAndLength=[cosSitaAndLength;getPersistanceLength(centroidInfo)];
    end
    for iNode=1:size(bioTree{iframe}.node,2)
        for iOut=1:size(bioTree{iframe}.node{iNode}.Out,2)
            if strcmp(bioTreeType,'simulation')
                centroidInfo=bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.traceCentroid;
            end
            if strcmp(bioTreeType,'real')
                measurment=bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.measurment;
                centroidInfo=[];
                for iCen=1:numel(measurment)
                    centroidInfo=[centroidInfo;measurment{iCen}(1).Centroid];
                end
            end
            cosSitaAndLength=[cosSitaAndLength;getPersistanceLength(centroidInfo)];
        end
    end
end
cosSitaAndLength=getAverageCosSita(cosSitaAndLength);
end
function cosSitaAndLength=getPersistanceLength(centroidInfo)
if size(centroidInfo,1)==1
    cosSitaAndLength=[];
    return
end
positionInfo=centroidInfo(1:end-1,:);
positionInfo1=centroidInfo(2:end,:);
slopeInfo=(positionInfo1(:,2)-positionInfo(:,2))./(positionInfo1(:,1)-positionInfo(:,1));
slopeSita=atan(slopeInfo);
slopeSita=getRightSlope(slopeSita);
cosSitaAndLength=[];
timeSeries=linspace(1,2000,500);
timeSeries=round(timeSeries);
dataSize=size(centroidInfo,1);
timeSeries(timeSeries>=dataSize-50)=[];
for tao=1:numel(timeSeries)
    positionInfo1=positionInfo(1:end-timeSeries(tao),:);
    positionInfo2=positionInfo(1+timeSeries(tao):end,:);
    lengthInfo=(positionInfo1-positionInfo2).^2;
    slopeSita1=slopeSita(1:end-timeSeries(tao),:);
    slopeSita2=slopeSita(1+timeSeries(tao):end,:);
    deltaSlopeSita=slopeSita2-slopeSita1;
    cosSitaAndLength=[cosSitaAndLength;[(lengthInfo(:,1)+lengthInfo(:,2)).^0.5,cos(deltaSlopeSita)]];
end
end
function cosSitaAndLength=getAverageCosSita(cosSitaAndLength)
cosSitaAndLength(isnan(cosSitaAndLength(:,2)),:)=[];
maxLength=fix(max(cosSitaAndLength(:,1)))+1;
if maxLength>100;
    maxLength=maxLength-100;
end
lengthRange=2:2:maxLength;
resultNew=[];
sitaInfo=cosSitaAndLength(:,2);
lengthInfo=cosSitaAndLength(:,1);
for i=1:numel(lengthRange)
    meanCosSita=mean(sitaInfo(lengthInfo<lengthRange(i)+1 & lengthInfo>=lengthRange(i)-1));
    resultNew=[resultNew;lengthRange(i),meanCosSita];
end
cosSitaAndLength=resultNew;
cosSitaAndLength(cosSitaAndLength(:,2)<=0,:)=[];
end
function slopeSita=getRightSlope(slopeSita)
slopeSita=slopeSita*180/pi;
deltaSita=diff(slopeSita);
deltaSita(deltaSita>90)=180-deltaSita(deltaSita>90);
deltaSita(deltaSita<-90)=180+deltaSita(deltaSita<-90);
if slopeSita(1)<=-180
    slopeSita(1)=slopeSita(1)+360;
end
if slopeSita(1)>180
    slopeSita(1)=slopeSita-360;
end
for i=2:numel(slopeSita)
    slopeSita(i)=slopeSita(i-1)+deltaSita(i-1);
    if slopeSita(i)>180
        slopeSita(i)=slopeSita(i)-360;
    end
    if slopeSita(i)<=-180
        slopeSita(i)=slopeSita(i)+360;
    end
end
slopeSita=slopeSita*pi/180;
end