number=67;
branchMatrix=result(number).branchMatrix;
branchMatrix(branchMatrix>=result(number).frameNum)=0;
branchMatrix(branchMatrix==-10000)=0;
branchMatrix(branchMatrix~=0)=1;

degreeNum=[];
regionHist=[];
distMatrix=branchMatrix;
for i=1:size(distMatrix,1)
    iLine=distMatrix(i,:);
    degreeNum(i)=sum(iLine);
end
maxDegree=max(degreeNum);
regionHist(1,:)=1:maxDegree;
for i=1:maxDegree
    regionHist(2,i)=numel(degreeNum(degreeNum>=i));
end

degreeNum=[];
for i=1:size(distMatrix,1)
iLine=distMatrix(i,:);
iLine=iLine/160000*frameStep/10/transRate*10^(-4);
iLine=1-iLine;
degreeNum(i)=1-prod(iLine);
end
regionHist=[];
regionHist(1,:)=0:0.0001:1;
for i=1:numel(regionHist)
regionHist(2,i)=numel(degreeNum(degreeNum>=regionHist(1,i)));
end
plot(regionHist(1,:),regionHist(2,:))