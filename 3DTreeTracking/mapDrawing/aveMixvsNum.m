function aveMixvsNum(bacteriaFrameInfo,color)
n=0;
for iframe=1:200:numel(bacteriaFrameInfo)
n=n+1;
multi=3;
finalInfo=bacteriaFrameInfo{iframe};
aveLen=mean(finalInfo.lengthInfo(finalInfo.lengthInfo<=80));
distMatrix=pdist2(finalInfo.centroidInfo,finalInfo.centroidInfo);
distMatrix(distMatrix<=aveLen*multi & distMatrix~=0)=1;
distMatrix(distMatrix~=1)=0;
distMatrix=logical(distMatrix);
degree=[];
branchIndex=bacteriaFrameInfo{iframe}.bacteriaInfo(:,5);
t=0;
attachNum=[];
for i=1:size(distMatrix,2)
    iLine=distMatrix(i,:);
    bacNum=numel(iLine(iLine==1));
    if bacNum>=10
        t=t+1;
        attachBranch=branchIndex(iLine);
        attachBranch=unique(attachBranch);
        attachNum(t,1)=numel(attachBranch);
    end
end
bacteriaNum(n)=numel(attachNum);
aveMatchNum(n)=mean(attachNum);
end
hold on; plot(bacteriaNum,aveMatchNum,color)
end

% load('D:\12-1数据总结\2013-04-13 F25 longTime jzy（WD 9）\newCluster\bacteriaFrameInfo.mat')
% aveMixvsNum(bacteriaFrameInfo,'r')
% load('D:\12-1数据总结\2013-05-05 jzy F25 longtime（WD 9）\newCluster\bacteriaFrameInfo.mat')
% aveMixvsNum(bacteriaFrameInfo,'r')
% load('D:\12-1数据总结\2013-05-07 jzy F15 longTime（WD 9）\newCluster\bacteriaFrameInfo.mat')
% aveMixvsNum(bacteriaFrameInfo,'b')
% load('D:\12-1数据总结\2013-03-31 F25 jzy longTIme（gp D）\newCluster\bacteriaFrameInfo.mat')
% aveMixvsNum(bacteriaFrameInfo,'r')
% load('D:\12-1数据总结\2012-12-15 jzy F1（WD 5）\newCluster\bacteriaFrameInfo.mat')
% aveMixvsNum(bacteriaFrameInfo,'g')
% load('F:\2013-05-18 wyj F1 longtime\newCluster\bacteriaFrameInfo.mat')
% aveMixvsNum(bacteriaFrameInfo,'g')
% load('F:\2013-05-31 jzy F1 longtime\newCluster\bacteriaFrameInfo.mat')
% aveMixvsNum(bacteriaFrameInfo,'g')