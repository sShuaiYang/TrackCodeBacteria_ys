function degreeDistributionVsTime(clusterTree,step,beginFrame)
beginOne=beginFrame;
timeInterval=[beginOne:step:size(clusterTree,2),size(clusterTree,2)];
for i=timeInterval
    [degreeDistribution]=getDistDistribution(clusterTree{i}.weightedMatrix);
    degreeDistribution(2,:)=degreeDistribution(2,:)/max(degreeDistribution(2,:));
    if size(degreeDistribution,2)>=100
        p=[];
        newY=[];
    else
        newY=[];
        p=[];
    end
    if ~isempty(degreeDistribution)
        hold on; plot(degreeDistribution(1,:), degreeDistribution(2,:),'r');
    end
end
end
function regionHist=getDistDistribution(distMatrix)
for i=1:size(distMatrix,2)
    iLine=distMatrix(i,:);
    degreeNum(i)=sum(iLine(iLine~=0));
end
maxDegree=ceil(max(degreeNum));
regionHist(1,:)=1:maxDegree;
n=0;
for i=1:maxDegree
    n=n+1;
    regionHist(2,n)=numel(degreeNum(degreeNum>=i));
end
end