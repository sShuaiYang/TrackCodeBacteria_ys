function allAve=timeLapsePlotForAveNewBranchRatio(bacteriaFrameInfo)
t=0;
frameNeed=100:100:14900;
for aimFrame=frameNeed
    finalInfo=bacteriaFrameInfo{aimFrame};
    aveLen=mean(finalInfo.lengthInfo(finalInfo.lengthInfo<=80));
    distMatrix=pdist2(finalInfo.centroidInfo,finalInfo.centroidInfo);
    distMatrix(distMatrix<=aveLen*3 & distMatrix~=0)=1;
    distMatrix(distMatrix~=1)=0;
    distMatrix=logical(distMatrix);
    degree=[];
    attachNum=[];
    branchIndex=bacteriaFrameInfo{aimFrame}.bacteriaInfo(:,5);
    for i=1:size(distMatrix,2)
        iLine=distMatrix(i,:);
        degree(i,1)=numel(iLine(iLine==1));
        attachBranch=branchIndex(iLine);
        bacBranch=branchIndex(i);
        attachBranch=unique(attachBranch);
        attachBranch(attachBranch==bacBranch)=[];
        attachNum(i,1)=numel(attachBranch);
    end
    attachNum=attachNum./degree;
    attachNum(degree==0)=0;
    t=t+1;
    if ~isempty(attachNum)
        allAve(t)=mean(attachNum);
    end
    allAve(isnan(allAve))=0;
end
allAve=cat(1,frameNeed/1200,allAve);
end