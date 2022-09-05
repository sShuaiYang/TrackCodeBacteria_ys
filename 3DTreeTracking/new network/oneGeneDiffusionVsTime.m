function [resultMatrix,bacMatrix]=oneGeneDiffusionVsTime(clusterTree,resultAll)
branchNum=size(resultAll(end).finalMatrix,2);
timeAll=size(resultAll,2);
resultMatrix=zeros(timeAll,branchNum);
bacMatrix=zeros(timeAll,branchNum);
timeStep=200*3;   %(second);
for iTime=1:timeAll
    finalMatrix=resultAll(iTime).finalMatrix;
    branchMatrix=clusterTree(iTime*timeStep/3).branchMatrix;
%     finalMatrix=finalMatrix/160000;
    finalMatrix=finalMatrix+branchMatrix;
    nodeInfo=clusterTree(iTime*timeStep/3).nodeInfo;
    for iBranch=1:branchNum
        linkInfo=finalMatrix(:,iBranch);
        if all(linkInfo==0)
            resultMatrix(iTime,iBranch)=0;
            bacMatrix(iTime,iBranch)=0;
        else
            centroidList=nodeInfo(linkInfo~=0,6:7);
            allCen=mean(centroidList,1);
            bacBranch=nodeInfo(:,5);
            bacBranchNum=numel(linkInfo(linkInfo==1));
            branchCentroid=nodeInfo(bacBranch==iBranch,6:7);
            allCen2=mean(branchCentroid,1);
            allBac=nodeInfo(:,3);
            allBacNum=numel(allBac(allBac~=-1));
            allBacNum=1;
            bacBranchNum=1;
            if bacBranchNum==0
                resultMatrix(iTime,iBranch)=mean(pdist2(centroidList,allCen).^2)^(0.5)/allBacNum;
            else
                resultMatrix(iTime,iBranch)=mean(pdist2(centroidList,allCen).^2)^(0.5)/allBacNum/bacBranchNum;
                bacMatrix(iTime,iBranch)=mean(pdist2(branchCentroid,allCen2).^2)^(0.5)/allBacNum/bacBranchNum;
            end
        end
    end
end
end
       