function [degreeResult] = mixNumPlot(bacteriaFrameInfo)
% generate direct link matrix
multi=3;
aimFrame=numel(bacteriaFrameInfo);
finalInfo=bacteriaFrameInfo{aimFrame};
aveLen=mean(finalInfo.lengthInfo(finalInfo.lengthInfo<=80));
distMatrix=pdist2(finalInfo.centroidInfo,finalInfo.centroidInfo);
distMatrix(distMatrix<=aveLen*multi & distMatrix~=0)=1;
distMatrix(distMatrix~=1)=0;
distMatrix=logical(distMatrix);

% branch mix
degree=[];
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
degreeResult=[];
for i=1:max(degree)
    degreeResult=[degreeResult;[i,mean(attachNum(degree==i))]];
end
end