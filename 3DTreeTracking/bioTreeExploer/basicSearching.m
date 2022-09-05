function [resultCount,resultDivision,resultLength,resultDAevent]=basicSearching(bioTree,treeGraph)
stackSize=1999;
xSize=1024;
ySize=1024;

windowSize=120;
lengthCutoff=200;
lengthStep=5;
timeCutoff=240;
timeStep=5;

minLength=0;
stepLength=1;
maxLength=70;

stepDuration=1:2:9000;
resultCount=bioTreeCounter(bioTree,windowSize);
resultDivision=divisionExploer(bioTree,lengthCutoff,lengthStep,timeCutoff,timeStep,xSize,ySize);
[resultLength,lengthTime]=countLength(bioTree,windowSize,minLength,stepLength,maxLength);
resultDAevent=DAeventExploer(bioTree,treeGraph,lengthTime,windowSize,minLength,stepLength,maxLength,stepDuration,lengthCutoff,lengthStep,timeCutoff,timeStep,xSize,ySize,stackSize);

end
function result=bioTreeCounter(bioTree,windowSize)
endTime=size(bioTree,2);
timeSeries=1:1:endTime-windowSize;
attechingNumber=zeros(1,endTime-windowSize);
attechingType0Number=zeros(1,endTime-windowSize);
attechingType1Number=zeros(1,endTime-windowSize);
detechingType0Number=zeros(1,endTime-windowSize);
detechingType1Number=zeros(1,endTime-windowSize);
detechingNumber=zeros(1,endTime-windowSize);
nodeType1Number=zeros(1,endTime-windowSize);
nodeType2Number=zeros(1,endTime-windowSize);
nodeType3Number=zeros(1,endTime-windowSize);
countRoot=[];
countRootType0=[];
countRootType1=[];
countLeaf=[];
countLeafType0=[];
countLeafType1=[];
countNodeType1=[];
countNodeType2=[];
countNodeType3=[];
countNodeOutIn=[];
for iframe=1:endTime
    rootType0=0;
    rootType1=0;
    leafType0=0;
    leafType1=0;
    if ~isempty(bioTree{iframe}.root)
        countRoot=[countRoot,size(bioTree{iframe}.root,2)];
        for iroot=1:size(bioTree{iframe}.root,2)
            if bioTree{iframe}.root{iroot}.isCluster==false
                rootType0=rootType0+1;
            end
            if bioTree{iframe}.root{iroot}.isCluster==true
                rootType1=rootType1+1;
            end
        end
        countRootType0=[countRootType0,rootType0];
        countRootType1=[countRootType1,rootType1];
    else
        countRoot=[countRoot,0];
        countRootType0=[countRootType0,0];
        countRootType1=[countRootType1,0];
    end
    if ~isempty(bioTree{iframe}.leavies)
        countLeaf=[countLeaf,size(bioTree{iframe}.leavies,2)];
        for ileaf=1:size(bioTree{iframe}.leavies,2)
            if bioTree{iframe}.leavies{ileaf}.isCluster==false
                leafType0=leafType0+1;
            end
            if bioTree{iframe}.leavies{ileaf}.isCluster==true
                leafType1=leafType1+1;
            end
        end
        countLeafType0=[countLeafType0,leafType0];
        countLeafType1=[countLeafType1,leafType1];
    else
        countLeaf=[countLeaf,0];
        countLeafType0=[countLeafType0,0];
        countLeafType1=[countLeafType1,0];
    end
    if ~isempty(bioTree{iframe}.node)
        type1=0;
        type2=0;
        type3=0;
        Out=0;
        In=0;
        for iNode=1:size(bioTree{iframe}.node,2)
            if bioTree{iframe}.node{iNode}.type==1
                type1=type1+1;
            end
            if bioTree{iframe}.node{iNode}.type==2
                type2=type2+1;
            end
            if bioTree{iframe}.node{iNode}.type==3
                type3=type3+1;
            end
            Out=Out+size(bioTree{iframe}.node{iNode}.Out,2);
            In=In+size(bioTree{iframe}.node{iNode}.In,2);
        end
        countNodeType1=[countNodeType1,type1];
        countNodeType2=[countNodeType2,type2];
        countNodeType3=[countNodeType3,type3];
        countNodeOutIn=[countNodeOutIn,Out-In];
    else
        countNodeType1=[countNodeType1,0];
        countNodeType2=[countNodeType2,0];
        countNodeType3=[countNodeType3,0];
        countNodeOutIn=[countNodeOutIn,0];
    end
end

for i=1: endTime-windowSize
    attechingNumber(i)=sum(countRoot(i:i+windowSize));
    attechingType0Number(i)=sum(countRootType0(i:i+windowSize));
    attechingType1Number(i)=sum(countRootType1(i:i+windowSize));
    detechingNumber(i)=sum(countLeaf(i:i+windowSize-1));
    detechingType0Number(i)=sum(countLeafType0(i:i+windowSize));
    detechingType1Number(i)=sum(countLeafType1(i:i+windowSize));
    nodeType1Number(i)=sum(countNodeType1(i:i+windowSize));
    nodeType2Number(i)=sum(countNodeType2(i:i+windowSize));
    nodeType3Number(i)=sum(countNodeType3(i:i+windowSize));
end
traceNumber=countTrace(bioTree,windowSize);
traceNumber0_temp=countTrace(bioTree,0);
traceNumber0=traceNumber0_temp(1:end-windowSize);
result.Number=[timeSeries',traceNumber',attechingNumber',attechingType0Number',attechingType1Number',detechingNumber',detechingType0Number',detechingType1Number',nodeType1Number',nodeType2Number',nodeType3Number',traceNumber0',(attechingNumber-detechingNumber)'];
result.Rate=[timeSeries',((attechingNumber./traceNumber)./windowSize)',((detechingNumber./traceNumber)./windowSize)',((nodeType1Number./traceNumber)./windowSize)',((nodeType2Number./traceNumber)./windowSize)',((nodeType3Number./traceNumber)./windowSize)'];
end
function traceNumber=countTrace(bioTree,windowSize)
endTime=size(bioTree,2);
traceNumber=zeros(1,endTime-windowSize);
for iframe=1:endTime
    if ~isempty(bioTree{iframe}.root)
        for iroot=1:size(bioTree{iframe}.root,2)
            if iframe-windowSize<=0
                traceNumberStartPoint=1;
            else
                traceNumberStartPoint=iframe-windowSize;
            end
            if bioTree{iframe}.root{iroot}.is2Node==false
                leafInfo=bioTree{iframe}.root{iroot}.leafInfo;
                if leafInfo(1)>=endTime-windowSize
                    traceNumberEndPoint=endTime-windowSize;
                else
                    traceNumberEndPoint=leafInfo;
                end
            end
            if bioTree{iframe}.root{iroot}.is2Node==true
                nodeInfo=bioTree{iframe}.root{iroot}.nodeInfo;
                if nodeInfo(1)>=endTime-windowSize
                    traceNumberEndPoint=endTime-windowSize;
                else
                    traceNumberEndPoint=nodeInfo(1);
                end
            end  
            traceNumber(traceNumberStartPoint:traceNumberEndPoint)= traceNumber(traceNumberStartPoint:traceNumberEndPoint)+1;
        end
    end
    if ~isempty(bioTree{iframe}.node)
        for iNode=1:size(bioTree{iframe}.node,2)
            for iOut=1:size(bioTree{iframe}.node{iNode},2)
                if iframe-windowSize<=0
                    traceNumberStartPoint=1;
                else
                    traceNumberStartPoint=iframe-windowSize;
                end
                if bioTree{iframe}.node{iNode}.Out{iOut}.is2Node==false
                    leafInfo=bioTree{iframe}.node{iNode}.Out{iOut}.leafInfo;
                    if leafInfo(1)>=endTime-windowSize
                        traceNumberEndPoint=endTime-windowSize;
                    else
                        traceNumberEndPoint=leafInfo;
                    end
                end
                if bioTree{iframe}.node{iNode}.Out{iOut}.is2Node==true
                    nodeInfo=bioTree{iframe}.node{iNode}.Out{iOut}.nodeInfo;
                    if nodeInfo(1)>=endTime-windowSize
                        traceNumberEndPoint=endTime-windowSize;
                    else
                        traceNumberEndPoint=nodeInfo(1);
                    end
                end              
                traceNumber(traceNumberStartPoint:traceNumberEndPoint)= traceNumber(traceNumberStartPoint:traceNumberEndPoint)+1;           
            end
        end
    end
end
end
function resultDivision=divisionExploer(bioTree,lengthCutoff,lengthStep,timeCutoff,timeStep,xSize,ySize)
divisionEvent=[];
corSpaceResult=[];
corTimeResult=[];
for iframe=1:size(bioTree,2)
    if ~isempty(bioTree{iframe}.node)
        for iNode=1:size(bioTree{iframe}.node,2)
            if bioTree{iframe}.node{iNode}.type==3
                centroidOut=[];
                for iOut=1:size(bioTree{iframe}.node{iNode}.Out,2)
                    centroidOut=[centroidOut;bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.measurment{1}(1).Centroid];
                end
                divisionEvent=[divisionEvent;[iframe,mean(centroidOut)]];
            end
        end
    end
end
traceNumber=countTrace(bioTree,0);
divisionEvent(:,4)=traceNumber(divisionEvent(:,1)')';
for itimeCutoff=[50,200,800,1600,3200]
    tempDistResult=spaceExploer(divisionEvent,lengthCutoff,lengthStep,itimeCutoff,xSize,ySize);
    corSpaceResult=[corSpaceResult,tempDistResult'];
end
for iDistCutoff=[50,128,256,512,1024]
    tempTimeResult=timeExploer(divisionEvent,timeCutoff,timeStep,iDistCutoff);
    corTimeResult=[corTimeResult,(tempTimeResult)'];
end
resultDivision{1}=[(0:lengthStep:lengthCutoff-lengthStep)',corSpaceResult];
resultDivision{2}=[(0:timeStep:timeCutoff-timeStep)',corTimeResult];
end
function corSpaceResult=spaceExploer(divisionEvent,lengthCutoff,lengthStep,timeThreshold,xSize,ySize)
centroidInWindowtemp=divisionEvent((divisionEvent(:,2)>=lengthCutoff & divisionEvent(:,2)<=xSize-lengthCutoff),:);
centroidInWindow=centroidInWindowtemp((centroidInWindowtemp(:,3)>=lengthCutoff & centroidInWindowtemp(:,3)<=ySize-lengthCutoff),:);
lengthBin=0:lengthStep:lengthCutoff-lengthStep;
areaBin=((lengthStep:lengthStep:lengthCutoff).^2-lengthBin.^2).*3.1416;
corSpaceResult=zeros(1,size(lengthBin,2));
for iEvent=1:size(centroidInWindow,1)
    distSquare=(divisionEvent(:,2)-centroidInWindow(iEvent,2)).^2+(divisionEvent(:,3)-centroidInWindow(iEvent,3)).^2;
    distSquareInCutOff=distSquare<=(lengthCutoff-lengthStep)*(lengthCutoff-lengthStep);
    centroidInWindowFinal=divisionEvent(distSquareInCutOff,:);
    distHist(:,2)=centroidInWindowFinal(:,4);
    distHist(:,1)=sqrt(distSquare(distSquareInCutOff));
    distHist(:,3)=abs(centroidInWindowFinal(:,1)-centroidInWindow(iEvent,1));
    distHistTime=distHist(distHist(:,3)<=timeThreshold,:);
    distResult=distCor(distHistTime,lengthCutoff,lengthStep);
    corSpaceResult=corSpaceResult+distResult;
    distHist=[];
end
corSpaceResult=((corSpaceResult./size(centroidInWindow,1))./areaBin)./timeThreshold;
end
function distResult=distCor(distHist,lengthCutoff,lengthStep)
lengthBin=0:lengthStep:lengthCutoff-lengthStep;
corSpaceResult=zeros(1,size(lengthBin,2));
count=0;
for iDist=1:size(distHist,1)
    distIndex=lengthBin<=distHist(iDist,1) & lengthBin>distHist(iDist,1)-lengthStep;
    corSpaceResult(distIndex)=corSpaceResult(distIndex)+1/distHist(iDist,2);
    if distHist(iDist,1)==0 && count==0
        corSpaceResult(1)=corSpaceResult(1)-1/distHist(iDist,2);
        count=count+1;
    end
end
distResult=corSpaceResult;
end
function corResult=timeExploer(divisionEvent,timeCutoff,timeStep,lengthThresohold)
eventInWindow=findEvent(divisionEvent,timeCutoff);
timeBin=0:timeStep:timeCutoff-timeStep;
corResult=zeros(1,size(timeBin,2));
for iEvent=1:size(eventInWindow,1)
    timeHistTemp=abs(divisionEvent(:,1)-eventInWindow(iEvent,1))';
    histIndex=timeHistTemp<=timeCutoff;
    timeHist=divisionEvent(histIndex,:);
    timeHist(:,1)=timeHist(:,1)-eventInWindow(iEvent,1);
    timeHist(:,5)=sqrt((timeHist(:,2)-eventInWindow(iEvent,2)).^2+(timeHist(:,3)-eventInWindow(iEvent,3)).^2);
    distIndex= timeHist(:,5)<=lengthThresohold;
    timeHistDist=timeHist(distIndex,:);
    tempResult=timeCor(timeHistDist,timeCutoff,timeStep);
    corResult=corResult+tempResult;
end
corResult=((corResult./size(eventInWindow,1))./timeStep)./(lengthThresohold*lengthThresohold);
end
function eventInWindow=findEvent(divisionEvent,timeCutoff)
firstEvent=divisionEvent(1,1);
lastEvent=divisionEvent(end,1);
timeEvent=divisionEvent(:,1)';
indexEvent= timeEvent>=firstEvent+timeCutoff & timeEvent<=lastEvent-timeCutoff;
eventInWindow=divisionEvent(indexEvent,:);
end
function corResult=timeCor(timeHist,timeCutoff,timeStep)
timeBin=0:timeStep:timeCutoff-timeStep;
corResult=zeros(1,size(timeBin,2));
count=0;
for iTime=1:size(timeHist,1)
    timeIndex=timeBin<=timeHist(iTime,1) & timeBin>timeHist(iTime,1)-timeStep;
    corResult(timeIndex)=corResult(timeIndex)+1/timeHist(iTime,4);
    if timeHist(iTime,1)==0 && count==0
        corResult(1)=corResult(1)-1/timeHist(iTime,4);
        count=count+1;
    end
end
end
function resultDAevent=DAeventExploer(bioTree,treeGraph,lengthTime,windowSize,minLength,stepLength,maxLength,stepDuration,lengthCutoff,lengthStep,timeCutoff,timeStep,xSize,ySize,stackSize)
attaching0Event=[];
detaching0Event=[];
attaching1Event=[];
detaching1Event=[];
traceNumber=countTrace(bioTree,0);
for iframe=1:size(bioTree,2);
%     if rem(iframe,stackSize)~=0
    if ~isempty(bioTree{iframe}.root)       
        for iroot=1:size(bioTree{iframe}.root,2)
            if bioTree{iframe}.root{iroot}.isCluster==false
                leafInfo=bioTree{iframe}.root{iroot}.leafInfo;
                attachingLength=bioTree{iframe}.root{iroot}.rootMeasurment.MajorAxisLength;
                attachingCentroid=bioTree{iframe}.root{iroot}.rootMeasurment.Centroid;
                attaching0Event=[attaching0Event;[iframe,leafInfo(1)-iframe,attachingLength,attachingCentroid,traceNumber(iframe)]];
            end
            if bioTree{iframe}.root{iroot}.isCluster==true                
                nodeIndex=bioTree{iframe}.root{iroot}.nodeIndex;
                desNode=getdescendants(treeGraph.nodes(nodeIndex),100);
                isLeaf=get(desNode,'Label');
                timeFrame=get(desNode,'UserData');
                attachingLength=bioTree{iframe}.root{iroot}.rootMeasurment.MajorAxisLength;
                attachingCentroid=bioTree{iframe}.root{iroot}.rootMeasurment.Centroid;
                for iList=1:size(isLeaf,1)
                    if strcmp(isLeaf{iList},'leaf') || strcmp(isLeaf{iList},'endLeaf')
                        duratioTime=timeFrame{iList}-iframe;
                        attaching1Event=[attaching1Event;[iframe,duratioTime,attachingLength,attachingCentroid,traceNumber(iframe)]];
                    end
                end                
            end
        end  
    end
    if ~isempty(bioTree{iframe}.leavies)      
        for ileaf=1:size(bioTree{iframe}.leavies,2)
            if bioTree{iframe}.leavies{ileaf}.isCluster==false
                rootInfo=bioTree{iframe}.leavies{ileaf}.rootInfo;
                duringTime=iframe-rootInfo(1);
                detachingLength=bioTree{iframe}.leavies{ileaf}.leafMeasurment.MajorAxisLength;
                detachingCentroid=bioTree{iframe}.leavies{ileaf}.leafMeasurment.Centroid;
                detaching0Event=[detaching0Event;[iframe,duringTime,detachingLength,detachingCentroid,traceNumber(iframe)]];
            end  
            if bioTree{iframe}.leavies{ileaf}.isCluster==true
                nodeIndex=bioTree{iframe}.leavies{ileaf}.nodeIndex;
                ancNode=getancestors(treeGraph.nodes(nodeIndex),100);
                isRoot=get(desNode,'Label');
                timeFrame=get(ancNode,'UserData');
                detachingLength=bioTree{iframe}.leavies{ileaf}.leafMeasurment.MajorAxisLength;
                detachingCentroid=bioTree{iframe}.leavies{ileaf}.leafMeasurment.Centroid;
                for iList=1:size(isRoot,1)
                    if strcmp(isRoot{iList},'root')
                        duratioTime=iframe-timeFrame{iList};
                        detaching1Event=[detaching1Event;[iframe,duratioTime,detachingLength,detachingCentroid,traceNumber(iframe)]];
                    end
                end
            end
        end
    end
%     end
end
resultDAevent0=surriveEvent(attaching0Event,detaching0Event);
resultDAevent1=surriveEvent(attaching1Event,detaching1Event);
resultDAevent.type0=resultDAevent0;
resultDAevent.type1=resultDAevent1;
end
function resultDAevent=surriveEvent(attachingEvent,detachingEvent)
stepTime=2;
endTime=10000;
stepDuration=1:stepTime:endTime;
surriveAttachingCount=zeros(1,size(stepDuration,2));
surriveDetachingCount=zeros(1,size(stepDuration,2));
for iTime=1:size(stepDuration,2)
    attachingSurrive=attachingEvent(attachingEvent(:,2)>=stepDuration(iTime),:);
    detachingSurrive=detachingEvent(detachingEvent(:,2)>=stepDuration(iTime),:);
    surriveAttachingCount(iTime)=size(attachingSurrive,1);
    surriveDetachingCount(iTime)=size(detachingSurrive,1);
end
resultDAevent.attachingSurrive=[stepDuration',surriveAttachingCount'];
resultDAevent.detachingSurrive=[stepDuration',surriveDetachingCount'];
end
function resultDAevent=satDAevent(bioTree,attachingEvent,detachingEvent,lengthTime,windowSize,minLength,stepLength,maxLength,stepDuration,lengthCutoff,lengthStep,timeCutoff,timeStep,xSize,ySize)
% corSpaceAResult=[];
% corTimeAResult=[];
% corSpaceDResult=[];
% corTimeDResult=[];
% [attachingHist,detachingHist]=histLength(bioTree,attachingEvent(:,1:5),detachingEvent(:,1:5),lengthTime,windowSize,minLength,stepLength,maxLength);
[timeDurationCount,durationCount,durationHist,durationLength]=durationSat(bioTree,detachingEvent(:,1:5),windowSize,stepDuration,minLength,stepLength,maxLength);
% corAttachingEvent=[attachingEvent(:,1),attachingEvent(:,4:6)];
% corDetachingEvent=[detachingEvent(:,1),detachingEvent(:,4:6)];
% for itimeCutoff=[5,20,80,320,1000]
%     tempDistResult=spaceExploer(corAttachingEvent,lengthCutoff,lengthStep,itimeCutoff,xSize,ySize);
%     corSpaceAResult=[corSpaceAResult,tempDistResult'];
%     tempDistResult=spaceExploer(corDetachingEvent,lengthCutoff,lengthStep,itimeCutoff,xSize,ySize);
%     corSpaceDResult=[corSpaceDResult,tempDistResult'];
% end
% for iDistCutoff=[30,60,120,240,480]
%     tempTimeResult=timeExploer(corAttachingEvent,timeCutoff,timeStep,iDistCutoff);
%     corTimeAResult=[corTimeAResult,(tempTimeResult)'];
%     tempTimeResult=timeExploer(corDetachingEvent,timeCutoff,timeStep,iDistCutoff);
%     corTimeDResult=[corTimeDResult,(tempTimeResult)'];
% end
% resultDAevent.spaceCorrelation.attaching=[(0:lengthStep:lengthCutoff-lengthStep)',corSpaceAResult];
% resultDAevent.timeCorrelation.attaching=[(0:timeStep:timeCutoff-timeStep)',corTimeAResult];
% resultDAevent.spaceCorrelation.detaching=[(0:lengthStep:lengthCutoff-lengthStep)',corSpaceDResult];
% resultDAevent.timeCorrelation.detaching=[(0:timeStep:timeCutoff-timeStep)',corTimeDResult];
% resultDAevent.length.attaching=attachingHist;
% resultDAevent.length.detaching=detachingHist;
% resultDAevent.duration.Hist2D=durationHist;
resultDAevent.duration.length=durationLength;
resultDAevent.duration.count=durationCount;
% resultDAevent.duration.timeCount=timeDurationCount;
end
function [attachingHist,detachingHist]=histLength(bioTree,attachingEvent,detachingEvent,lengthTime,windowSize,minLength,stepLength,maxLength)
histLengthBin=minLength:stepLength:maxLength-stepLength;
histTimeBin=0:windowSize:size(bioTree,2)-windowSize;
attachingHist2D=zeros(size(histLengthBin,2),size(histTimeBin,2));
detachingHist2D=zeros(size(histLengthBin,2),size(histTimeBin,2));
for iTime=1:size(histTimeBin,2)
    attachingNorm=[];
    detachingNorm=[];
    attachingLengthHisinTimeBin=zeros(1,size(histLengthBin,2));
    detachingLengthHisinTimeBin=zeros(1,size(histLengthBin,2));
    attachingEventinTimeBin=attachingEvent((attachingEvent(:,1)>=(iTime-1)*windowSize & attachingEvent(:,1)<iTime*windowSize),:);
    detachingEventinTimeBin=detachingEvent((detachingEvent(:,1)>=(iTime-1)*windowSize & detachingEvent(:,1)<iTime*windowSize),:);
    for iAttaching=1:size(attachingEventinTimeBin,1)
        lengthCountAttaching=size(find(lengthTime{attachingEventinTimeBin(iAttaching,1)}<=attachingEventinTimeBin(iAttaching,3)+stepLength/2 & lengthTime{attachingEventinTimeBin(iAttaching,1)}>=attachingEventinTimeBin(iAttaching,3)-stepLength/2),2);
        attachingNorm=[attachingNorm,lengthCountAttaching];
    end
    for iDetaching=1:size(detachingEventinTimeBin,1)
        lengthCountdetaching=size(find(lengthTime{detachingEventinTimeBin(iDetaching,1)}<=detachingEventinTimeBin(iDetaching,3)+stepLength/2 & lengthTime{detachingEventinTimeBin(iDetaching,1)}>=detachingEventinTimeBin(iDetaching,3)-stepLength/2),2);
        detachingNorm=[detachingNorm,lengthCountdetaching];
    end
    attachingEventinTimeBin=[attachingEventinTimeBin,attachingNorm'];
    detachingEventinTimeBin=[detachingEventinTimeBin,detachingNorm'];
    if ~isempty(attachingEventinTimeBin)
        for iLength=1:size(histLengthBin,2)
            attachingEventinLengthBin=attachingEventinTimeBin(attachingEventinTimeBin(:,3)>=(iLength-1)*stepLength & attachingEventinTimeBin(:,3)<iLength*stepLength,:);
            if ~isempty(attachingEventinLengthBin)
                attachingLengthHisinTimeBin(iLength)=mean(1./attachingEventinLengthBin(:,6));
            else
                attachingLengthHisinTimeBin(iLength)=0;
            end
        end
    end
    
    if ~isempty(detachingEventinTimeBin)
        for iLength=1:size(histLengthBin,2)
            detachingEventinLengthBin=detachingEventinTimeBin(detachingEventinTimeBin(:,3)>=(iLength-1)*stepLength & detachingEventinTimeBin(:,3)<iLength*stepLength,:);
            if ~isempty(detachingEventinLengthBin)
                detachingLengthHisinTimeBin(iLength)=mean(1./detachingEventinLengthBin(:,6));
            else
                detachingLengthHisinTimeBin(iLength)=0;
            end
        end
    end
    
    attachingHist2D(:,iTime)=attachingLengthHisinTimeBin;
    detachingHist2D(:,iTime)=detachingLengthHisinTimeBin;
end

attachingHist=attachingHist2D;
detachingHist=detachingHist2D;
end
function [resultLength,lengthTime]=countLength(bioTree,windowSize,minLength,stepLength,maxLength)
endTime=size(bioTree,2);
lengthTime{endTime}=[];
beforeDivisionLength{endTime}=[];
afterDivisionLength{endTime}=[];
timeBin=0:windowSize:endTime-windowSize;
lengthBin=minLength:stepLength:maxLength;
for iframe=1:endTime
    if ~isempty(bioTree{iframe}.root)
        for iroot=1:size(bioTree{iframe}.root,2)
            for iTrace=1:size(bioTree{iframe}.root{iroot}.traceInfo.measurment,2)
                lengthTime{iframe+iTrace-1}=[lengthTime{iframe+iTrace-1},bioTree{iframe}.root{iroot}.traceInfo.measurment{iTrace}.MajorAxisLength];
            end
        end
    end
    if ~isempty(bioTree{iframe}.node)
        for iNode=1:size(bioTree{iframe}.node,2)
            for iOut=1:size(bioTree{iframe}.node{iNode}.Out,2)
                for iTrace=1:size(bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.measurment,2)
                    lengthTime{iframe+iTrace-1}=[lengthTime{iframe+iTrace-1},bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.measurment{iTrace}.MajorAxisLength];
                end
            end
            if bioTree{iframe}.node{iNode}.type==3
                if size(bioTree{iframe}.node{iNode}.In,2)==1 && size(bioTree{iframe}.node{iNode}.Out,2)==2
                    if bioTree{iframe}.node{iNode}.In{1}.isNode==false
                        rootInfo=bioTree{iframe}.node{iNode}.In{1}.rootInfo;
                        beforeDivisionLength{iframe}=[beforeDivisionLength{iframe},sum([bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.measurment{end}.MajorAxisLength])];
                        afterDivisionLength{iframe}=[afterDivisionLength{iframe},mean([bioTree{iframe}.node{iNode}.Out{1}.traceInfo.measurment{1}.MajorAxisLength,bioTree{iframe}.node{iNode}.Out{2}.traceInfo.measurment{1}.MajorAxisLength])];
                    end
                    if bioTree{iframe}.node{iNode}.In{1}.isNode==true
                        nodeInfo=bioTree{iframe}.node{iNode}.In{1}.nodeInfo;
                        beforeDivisionLength{iframe}=[beforeDivisionLength{iframe+iTrace-1},sum([bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.measurment{end}.MajorAxisLength])];
                        afterDivisionLength{iframe}=[afterDivisionLength{iframe},mean([bioTree{iframe}.node{iNode}.Out{1}.traceInfo.measurment{1}.MajorAxisLength,bioTree{iframe}.node{iNode}.Out{2}.traceInfo.measurment{1}.MajorAxisLength])];
                    end
                end
            end
        end
    end
end
averageLength=zeros(1,size(timeBin,2));
averageLengthBeforDivision=zeros(1,size(timeBin,2));
averageLengthAfterDivision=zeros(1,size(timeBin,2));
hist2DLength=zeros(size(lengthBin,2),size(timeBin,2));
for iTime=1:size(timeBin,2)
    timeLenfthTemp=[];
    lengthBeforeDivision=[];
    lengthAfterDivision=[];
    countLength=[];
    for iframe=(iTime-1)*windowSize+1:iTime*windowSize;
        timeLenfthTemp=[timeLenfthTemp,lengthTime{iframe}];
        if ~isempty(beforeDivisionLength{iframe})
            lengthBeforeDivision=[lengthBeforeDivision,beforeDivisionLength{iframe}];
        end
        if ~isempty(afterDivisionLength{iframe})
            lengthAfterDivision=[lengthAfterDivision,afterDivisionLength{iframe}];
        end
    end
    averageLength(iTime)=mean(timeLenfthTemp);
    countLength=hist(timeLenfthTemp,lengthBin);
    hist2DLength(:,iTime)=countLength./sum(countLength);
    if ~isempty(lengthBeforeDivision)
        averageLengthBeforDivision(iTime)=mean(lengthBeforeDivision);
    end
    if ~isempty(lengthAfterDivision)
        averageLengthAfterDivision(iTime)=mean(lengthAfterDivision);
    end
end
timeDivsion=timeBin(averageLengthBeforDivision>0);
averageLengthBeforDivision=averageLengthBeforDivision(averageLengthBeforDivision>0);
averageLengthAfterDivision=averageLengthAfterDivision(averageLengthAfterDivision>0);
resultLength.averageLength(:,1)=(timeBin./windowSize)+1;
resultLength.averageLength(:,2)=averageLength;
resultLength.hist2D=hist2DLength;
resultLength.divisionLength(:,1)=(timeDivsion./windowSize)+1;
resultLength.divisionLength(:,2)=averageLengthBeforDivision;
resultLength.divisionLength(:,3)=averageLengthAfterDivision;
end
function [timeDurationCount,durationCount,durationHist,durationLength]=durationSat(bioTree,detachingEvent,windowSize,stepDuration,minLength,stepLength,maxLength)
histTimeBin=0:windowSize:size(bioTree,2)-windowSize;
dutaionTimeBin=stepDuration(1:end);
timeDurationCount=zeros(size(dutaionTimeBin,2),size(histTimeBin,2));
for iTime=1:size(histTimeBin,2)
    detachingEnventinWindow=detachingEvent((detachingEvent(:,1)>=(iTime-1)*windowSize & detachingEvent(:,1)<(iTime)*windowSize),:);
    [durationCount,~,~]=durationHistCount(detachingEnventinWindow,stepDuration,minLength,stepLength,maxLength);
    timeDurationCount(:,iTime)=durationCount./sum(durationCount);
end
[durationCount,durationHist,durationLength]=durationHistCount(detachingEvent,stepDuration,minLength,stepLength,maxLength);
durationCount=[stepDuration',durationCount'];
durationLength=[dutaionTimeBin',durationLength'];
end
function [durationCount,durationHist,durationLength]=durationHistCount(detachingEvent,stepDuration,minLength,stepLength,maxLength)
histLengthBin=minLength:stepLength:maxLength-stepLength;
dutaionTimeBin=stepDuration(2:end)-stepDuration(1:end-1);
durationHist=zeros(size(histLengthBin,2),size(dutaionTimeBin,2));
durationLength=zeros(1,size(dutaionTimeBin,2));
durationCount=zeros(1,size(dutaionTimeBin,2));
for iDuration=1:size(stepDuration,2)
    eventInDuration=detachingEvent(detachingEvent(:,2)>=stepDuration(iDuration),:);
    [countTemp,~]=hist(eventInDuration(:,3),histLengthBin);
    durationCount(iDuration)=size(eventInDuration(:,3),1);
    durationHist(:,iDuration)=countTemp./sum(countTemp);
    durationLength(iDuration)=sum((durationHist(:,iDuration))'.*histLengthBin);
end
end








