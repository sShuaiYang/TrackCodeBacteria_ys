function [velocityResult,lengthTime]=lengthAndVelocityVsTime(bioTree,bacteriaFrameInfo)
%% calculate ave length vs time
step=20;
timeInterval=bioTree{1}.scaleInfo.timeInterval;
scaleInfo=bioTree{1}.scaleInfo.scaleInfo;
n=0;
for i=step:step:numel(bioTree)
    n=n+1;
    lengthTime(n,1)=i*timeInterval/60;
    lengthInfo=bacteriaFrameInfo{i}.lengthInfo;
    lengthTime(n,2)=mean(lengthInfo)*scaleInfo;
end

%% calculate mean velocity vs time
stepV=5;
vInfo=[];
for iframe=1:numel(bioTree)
    for iRoot=1:size(bioTree{iframe}.root,2)
        traceNum=numel(bioTree{iframe}.root{iRoot}.traceInfo.measurment);
        if traceNum>=n+1;
            positionInfo=[];
            for iTrace=1:stepV:traceNum
                positionInfo((iTrace-1)/stepV+1,:)=[(iTrace+iframe-1)*timeInterval/60,bioTree{iframe}.root{iRoot}.traceInfo.measurment{iTrace}(1).Centroid(1,:)];
            end
            velocityInfo=sum((positionInfo(2:end,2:3)-positionInfo(1:end-1,2:3)).^2,2).^(1/2);
            velocityInfo=velocityInfo*scaleInfo/(stepV*timeInterval);   % um/s
            velocityInfo=cat(2,positionInfo(1:end-1,1),velocityInfo);
            vInfo=[vInfo;velocityInfo];
        end
    end
    for iNode=1:size(bioTree{iframe}.node,2)
        for iOut=1:size(bioTree{iframe}.node{iNode}.Out,2)
            traceNum=numel(bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.measurment);
            if traceNum>=n+1;
                positionInfo=[];
                for iTrace=1:stepV:traceNum
                    positionInfo((iTrace-1)/stepV+1,:)=[(iTrace+iframe-1)*timeInterval/60,bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.measurment{iTrace}(1).Centroid(1,:)];
                end
                velocityInfo=sum((positionInfo(2:end,2:3)-positionInfo(1:end-1,2:3)).^2,2).^(1/2);
                velocityInfo=velocityInfo*scaleInfo/(stepV*timeInterval);   % um/s
                velocityInfo=cat(2,positionInfo(1:end-1,1),velocityInfo);
                vInfo=[vInfo;velocityInfo];
            end
        end
    end
end
maxTime=max(vInfo(:,1));
timeSeries=linspace(0,maxTime,600);
velocityAll=vInfo(:,2);
timeAll=vInfo(:,1);
for iTime=1:numel(timeSeries)-1
    velocityResult(iTime,1)=(timeSeries(iTime)+timeSeries(iTime+1))/2;
    velocityResult(iTime,2)=mean(velocityAll(timeAll<timeSeries(iTime+1) & timeAll>=timeSeries(iTime)));
end