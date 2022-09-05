function histResult=histVelocity(velocityList) % displacement weighted hisotgram
binSpace=logspace(-3,2,81); % the las number must be an odd, and the data nmber is (n-1)/2
angleSpace=linspace(-180,180,37); %same as above
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